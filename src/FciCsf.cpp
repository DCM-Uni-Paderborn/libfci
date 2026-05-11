/* This program is free software. It comes without any warranty, to the extent
 * permitted by applicable law. You may use it, redistribute it and/or modify
 * it, in whole or in part, provided that you do so at your own risk and do not
 * hold the developers or copyright holders liable for any claim, damages, or
 * other liabilities arising in connection with the software.
 *
 * Developed by Gerald Knizia, 2010--2012.
 */

#include "Fci.h"

// recursively add patterns which add up to nElecLeft electons in the nOrbLeft
// last orbitals.
static void ListBitPatternsR(FOrbPat *&pOut, FOrbPat OldPat, uint nOrb,
                             int nElecLeft, uint iFirstOrb) {
  if (nElecLeft == 0 && iFirstOrb <= nOrb) {
    *pOut = OldPat;
    ++pOut;
    return;
  };

  for (uint iOrb = iFirstOrb; iOrb < nOrb; ++iOrb) {
    FOrbPat NewPat = OldPat | (FOrbPat(1) << iOrb);
    ListBitPatternsR(pOut, NewPat, nOrb, nElecLeft - 1, iOrb + 1);
  };
}

uint64_t binomial_coefficient(uint N, uint k);

// make an explicit list of all bit patterns with nElec bits != 0 and all
// bits > nOrb zero.
void ListBitPatterns(FOrbPat *&pBeg, FOrbPat *&pEnd, uint nElec, uint nOrb,
                     FMemoryStack &Mem) {
  std::size_t nPatterns = binomial_coefficient(nOrb, nElec);
  Mem.Alloc(pBeg, nPatterns);
  pEnd = pBeg + nPatterns;
  FOrbPat *pOut = pBeg;
  ListBitPatternsR(pOut, 0, nOrb, nElec, 0);
  assert(pOut == pEnd);
};

// S,M: total spin *after* s/m have been added.
// (note: we might want to tabulate these if they become a problem)
static double ClebschGordan12(uint S_, int M_, bool sPlus, bool mPlus) {
  if ((uint)std::abs(M_) > S_)
    return 0.0;

  // all spins given in units of HALF integers.
  double S = S_ / 2.0, M = M_ / 2.0, m = mPlus ? +.5 : -.5;

  // purple book, p. 56, eq. 2.6.5, 2.6.6
  if (sPlus)
    return std::sqrt((S + 2 * m * M) / (2. * S));
  else
    return (-2 * m) * std::sqrt((S + 1 - 2 * m * M) / (2. * (S + 1)));
};

// calculate coefficient of determimant Det of CSF identified by the
// genealogical coupling vector Coup.
static double GetCsfDetCoeff(FOrbPat Det, FOrbPat Coup, uint S, int M,
                             uint nOrb) {
  // intermediate spin/m to which the current determinant is already coupled.
  uint S_ = 0;
  int M_ = 0;
  // coefficient of the current intermediate determinant.
  double Coeff = 1.;
  for (uint iOrb = 0; iOrb < nOrb; ++iOrb) {
    bool sPlus = (Coup >> iOrb) & 1, mPlus = (Det >> iOrb) & 1;
    S_ += 2 * (int)sPlus - 1;
    M_ += 2 * (int)mPlus - 1;
    Coeff *= ClebschGordan12(S_, M_, sPlus, mPlus);
    if (Coeff == 0.)
      return 0.;
  };
  return Coeff;
};

// Absorb the relative sign of permuting from the determinant order as
// alpha/beta bettern with respect to a fixed order of open-shell
// orbitals, into the string order, in which we have first all alpha
// orbitals (with their respective open orbitals), followed by all beta
// orbitals (with their respective open orbitals)
//
// pDetCoeffs: a nDets x nVecs matrix (need not be CSFs; just treated as a list
// of determinant coefficients).
static void AbsorbSignsForStringOrder1(double *pDetCoeffs, FOrbPat *pDetPat,
                                       uint nDets, uint nVecs, uint nAlpha,
                                       uint nOpen) {
  for (uint iDet = 0; iDet < nDets; ++iDet) {
    int RelativeSign = 1, nAlphaLeft = nAlpha;
    FOrbPat Det = pDetPat[iDet];
    for (uint iOrb = 0; iOrb < nOpen; ++iOrb) {
      if (Det & (FOrbPat(1) << iOrb)) {
        nAlphaLeft -= 1;
      } else {
        if ((nAlphaLeft & 1) == 1)
          RelativeSign *= -1;
      }
    }

    if (RelativeSign == -1)
      for (uint iCsf = 0; iCsf < nVecs; ++iCsf)
        pDetCoeffs[iDet + nDets * iCsf] *= -1;
  }
}

static void RemoveIllicitCouplings(FOrbPat *pCoupBeg, FOrbPat *&pCoupEnd,
                                   uint nOpen) {
  // remove illicit couplings: that is all couplings for which the
  // partial sum of the first N terms is < 0. These are illegal because
  // the sum of the first N terms represents the total spin after N steps.
  FOrbPat *pCoupEndOrig = pCoupEnd;

  pCoupEnd = pCoupBeg;
  for (FOrbPat *pCoup = pCoupBeg; pCoup != pCoupEndOrig; ++pCoup) {
    int iPartialSum = 0;
    for (uint iOrb = 0; iOrb < nOpen; ++iOrb) {
      iPartialSum += (2 * ((*pCoup >> iOrb) & 1) - 1);
      if (iPartialSum < 0)
        // illicit. delete configuration and go on.
        break;
    };
    if (iPartialSum >= 0) {
      *pCoupEnd = *pCoup;
      ++pCoupEnd;
    }
  };
}

// This implements the genealogical coupling scheme for generating the expansion
// of spin eigenfunctions (for configuration state functions) in terms of
// determinants. See purple book p. 54ff.
//
// Return values:
//    - pCsfCoeffs: nDets x nCsf matrix.
//    - nDets, nCsfs: number of determinants and csfs (respectively) for nOpen
//      orbitals and given S/M quantum numbers.
//    - pDetPat: determinant patterns nDets refers to. Each one is understood
//      as a vector of length nOpen, with bit #i indicating that open-shell
//      orbital #i is alpha. All open orbitals which are not alpha are beta.
// Note on inputs:
//    - If AbsorbSignsForStringOrder is given, additional signs are absorbed
//      into the CSFs with the assumption that instead of the default ordering,
//      where the spatial orbital order is fixed and alpha/beta varies,
//      we first have all alpha orbitals and then all beta orbitals, with the
//      spatial orbital order varying. The latter scheme is convenient in string
//      determinant methods.
//
// Note: Output data is allocated on Mem. pCsfCoeffs may not be the first
//       allocation on Mem! So keep base pointers separately.
void MakeCsfExpansion(double *&pCsfCoeffs, FOrbPat *&pDetPat, uint &nDets,
                      uint &nCsfs, uint S, uint M, uint nOpen,
                      bool AbsorbSignsForStringOrder, FMemoryStack &Mem) {
  if (nOpen < S || nOpen < M) {
    // not enough orbitals to produce the desired spin. no CSFs.
    pCsfCoeffs = 0;
    nCsfs = 0;
    assert(0); // should probably not be called in this case.
    return;
  }

  // nOpen = nAlpha + nBeta
  //     M = nAlpha - nBeta
  assert((nOpen + M) % 2 == 0 && (nOpen - M) % 2 == 0);
  assert(M <= nOpen);
  uint nAlpha = (nOpen + M) / 2;
  // nBeta = (nOpen - M)/2;

  FOrbPat *pDetEnd, *pCoupBeg, *pCoupEnd;
  ListBitPatterns(pDetPat, pDetEnd, nAlpha, nOpen, Mem);
  nDets = pDetEnd - pDetPat;

  // nPlus - nMinus = S;
  // nPlus + nMinus = nOpen;
  assert((nOpen + S) % 2 == 0);
  ListBitPatterns(pCoupBeg, pCoupEnd, (nOpen + S) / 2, nOpen, Mem);
  // ^- this lists all the bit patterns with nPlus bits positive;
  //    bit #i set means: |+>, bit #i not set means |->. E.g.,
  //    1101 denotes the |++-+> configuration (which has total spin
  //    1/2 + 1/2 - 1/2 + 1/2 = 2/2).
  //    In our current scheme we will use them as genealogical coupling
  //    vectors, each of them identifying a CSF with total spin s.

  // remove illicit couplings: that is all couplings for which the
  // partial sum of the first N terms is < 0. These are illegal because
  // the sum of the first N terms represents the total spin after N steps.
  RemoveIllicitCouplings(pCoupBeg, pCoupEnd, nOpen);
  nCsfs = pCoupEnd - pCoupBeg;

  Mem.Alloc(pCsfCoeffs, nDets * nCsfs);
  for (uint iCoup = 0; iCoup < pCoupEnd - pCoupBeg; ++iCoup)
    for (uint iDet = 0; iDet < pDetEnd - pDetPat; ++iDet) {
      pCsfCoeffs[iDet + nDets * iCoup] =
          GetCsfDetCoeff(pDetPat[iDet], pCoupBeg[iCoup], S, M, nOpen);
    };

  if (AbsorbSignsForStringOrder && !g_BosonicSigns)
    AbsorbSignsForStringOrder1(pCsfCoeffs, pDetPat, nDets, nCsfs, nAlpha,
                               nOpen);
}

// kate: indent-width 3
