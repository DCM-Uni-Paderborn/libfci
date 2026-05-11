/* This program is free software. It comes without any warranty, to the extent
 * permitted by applicable law. You may use it, redistribute it and/or modify
 * it, in whole or in part, provided that you do so at your own risk and do not
 * hold the developers or copyright holders liable for any claim, damages, or
 * other liabilities arising in connection with the software.
 *
 * Developed by Gerald Knizia, 2010--2012.
 */

#include "Fci.h"

// references:
//   [1]: Knowles & Handy, CPL 111 315 (1984):
//        "A new determinant based full configuration interaction algorithm"

// algorithm 1: for R := H|Psi> contraction, with
//      H = \sum_ijkl <ij|kl> c^i c^j c_l c_k
// (only two electron part since 1e^- can be absorbed).
// We use the RI in the N-2 electron space. Introduce K:
//      r[K] := \sum_ij \sum_K <R|c^i c^j|K> I[ij,K]
//      I[ij,K] = \sum_kl <ij|kl> <K|c_l c_k|I> c_I
//
// 1) form intermediate space determinants connecting to K
//      I[kl,K] = c_l c_k |I> c_I
// 2) contract with two-electron integrals
//      I'[ij,K] = <ij|kl> I[kl,K]
// 3) form contributions to residual. Set r := 0, then go
//      R[J] = c^i c^j I'[ij,K]
//
// note: may need to keep track of sign changes due to the fermion stuff.

// algorithm 2: re-cast into form of single-replacement operators.
//      H  = \sum_ijkl <ij|kl> c^i c^j c_l c_k
//         = \sum_ijkl <ij|kl> c^i c_k c^j c_l + delta[kj] ..
//         = \sum_ijkl (ik|jl) c^i c_k c^j c_l + delta[kj] ..
//         = \sum_ijkl (ij|kl) c^i c_j c^k c_l + delta[kj] ..
// then form RI in the E_ij space:
//      <R|H2e|I> c_I =  <R|c^i_j|K> (ij|kl) <K|c^k_l|I> C_I
// this has the advantage that the replacements are spin-pure,
// as explained by Knowles (when c^k_l = c^kA_lA + c^kB_lB, then
// each acts on only one component of the spin space).
//
// Note that the 2-destruction (or even the 2-creation form) may
// still be better, because the N-2 (or N+2) electron intermediate space
// may have less configurations. However, the current form is what
// we actually need for the Hubbard model, so we go with that.

FHamiltonianTerm::~FHamiltonianTerm() {}

std::ostream &operator<<(std::ostream &out, FSubstResult const &o) {
  out << format("-> c^%i c_%i -> %i%s") % (int)o.k % (int)o.l % (int)o.sign %
             FmtPat(o.Str, 9);
  return out;
};

// form a sparse list of all |K> which can be reached by applying
// c^k_l on string |I>. These are supposed to have the same spin
// projection (e.g., c^{k\alpha}_{l\alpha} |I\alpha>.
//
// data is created for all k,l in 0 .. Adr.nOrb-1;
// The nOrb x nOrb matrix p1OpMatrix, if supplied, is used to
// decide whether a specific kl entry is to be processed (if
// the matrix element is 0, it is not)
void FormStringSubstsForSpin(FSubstResult *&pResult, uint &nEntries,
                             FScalar const *p1OpMatrix,
                             FOrbStringAdrTable const &Adr, FOrbPat I,
                             FMemoryStack &Mem) {
  uint nOrb = Adr.nOrb();
  Mem.Alloc(pResult, nOrb * nOrb);
  nEntries = 0;

  for (uint l = 0; l < nOrb; ++l) {
    FOrbPat J = I, mask_l;
    mask_l = (FOrbPat(1) << l);
    if ((I & mask_l) == 0)
      continue; // c_l annihilates |I>.
    J &= ~mask_l;
    char sign1 = StringParityBeforePos(J, l);

    for (uint k = 0; k < nOrb; ++k) {
      FOrbPat K = J, mask_k;
      mask_k = (FOrbPat(1) << k);
      K |= mask_k;
      if ((J & mask_k) != 0 ||
          (p1OpMatrix && std::abs(p1OpMatrix[k + nOrb * l]) < ThrNegl))
        continue; // c^k annihilates c_l|I> or h_kl is zero.

      FSubstResult &r = pResult[nEntries];
      r.k = k;
      r.l = l;
      r.sign = 1 - 2 * (sign1 ^ StringParityBeforePos(K, k));
      r.Str = K;
      r.iStr = Adr(K);
      ++nEntries;
    }
  }
  if (g_BosonicSigns)
    // remove the signs.
    for (uint i = 0; i < nEntries; ++i)
      pResult[i].sign = 1;
};

// apply an orbital creation or destruction operator of spin X on a determinant
// string of spin X. Returns resulting determinant in Out (which will have N-1
// or N+1 electrons). Return values:
//     +1: operation has positive parity,
//     -1: operation negative partiy,
//      0: operation annihilates In
signed char ApplyCop1(FOrbPat &Out, FOrbPat In, uint iOrb,
                      int iCreateOrDestroy) {
  FOrbPat mask_l = (FOrbPat(1) << iOrb);
  Out = 0;
  if (iCreateOrDestroy == -1) {
    if ((In & mask_l) == 0)
      return 0; // c_l annihilates |I>.
    Out = In & ~mask_l;
  } else {
    if ((In & mask_l) != 0)
      return 0; // c^l annihilates |I>.
    Out = In | mask_l;
  }
  return 1 - 2 * StringParityBeforePos(Out, iOrb);
};

// applies the 1-electron operator p1Op, given by a nOrb x nOrb matrix,
// on spin branch 1 of the vector given in
//    pCoeffs[iAdr1 * St1 + iAdr2 * St2] where
// iAdr1/2 are string addresses as expressed by the Adr1/Adr2 objects,
// and St1/St2 are the coefficient strides.
//
// Results added to pResult[iAdr1 * St1 + iAdr2 * St2].
void Apply1eOp(FScalar *pResult, FScalar *pCoeff, double Prefactor,
               FScalar const *p1OpMatrix, FOrbStringAdrTable const &Adr1,
               FOrbStringAdrTable const &Adr2, uint St1, uint St2,
               FMemoryStack &Mem_) {
  void *pBeginOfMemory = Mem_.Alloc(0);
  uint nOrb = Adr1.nOrb();

  if (Adr1.nElec() == 0)
    return;
  FOrbStringAdrTable
      // address table for intermediate states with N-1 electrons:
      //    r[J] += \sum_(K,I) <J|c^k1|K> h[k,l] <K|c_l1|I> c[I]
      // note that contrary to the 2e operator version, we only need one
      // spin state here.
      AdrK1(Adr1.nElec() - 1, nOrb);
  //     cout << format("AdrK1: nElec = %i  nOrb = %i  nStr = %i") %
  //     AdrK1.nElec() % AdrK1.nOrb() % AdrK1.nStr() << std::endl;
  FAdr *piAdr1;
  signed char *pSigns;
  Mem_.Alloc(piAdr1, nOrb * AdrK1.nStr());
  Mem_.Alloc(pSigns, nOrb * AdrK1.nStr());
  Mem_.Align(64);

  {
    FMemoryStackArray MemStacks(Mem_);

#pragma omp parallel for schedule(dynamic)
    for (FAdr iStr1 = 0; iStr1 < AdrK1.nStr(); ++iStr1) {
      FOrbPat K1 = AdrK1.MakePattern(iStr1);
      for (uint k = 0; k < nOrb; ++k) {
        FAdr iKk = k + nOrb * iStr1;
        FOrbPat I1;
        pSigns[iKk] = ApplyCop1(I1, K1, k, +1);
        if (pSigns[iKk] != 0)
          piAdr1[iKk] = Adr1(I1);
        else
          piAdr1[iKk] = 0; // that's a valid index even for sign == 0.
      }
    }

// iterate through |K2> configurations (2 = alpha or beta)
#pragma omp parallel for schedule(dynamic)
    for (FAdr iStr2 = 0; iStr2 < Adr2.nStr(); ++iStr2) {
      FMemoryStack &Mem = MemStacks.GetStackOfThread();
      uint const nTgtBlkK = 64;
      double *pCoeff_ = pCoeff + St2 * iStr2, *pResult_ = pResult + St2 * iStr2,
             *pInpK, *pOutK;
      Mem.Alloc(pInpK, nOrb * nTgtBlkK);
      Mem.Alloc(pOutK, nOrb * nTgtBlkK);

      for (FAdr iBlockBegK1 = 0; iBlockBegK1 < AdrK1.nStr();
           iBlockBegK1 += nTgtBlkK) {
        uint nBlkK =
            std::min(iBlockBegK1 + nTgtBlkK, AdrK1.nStr()) - iBlockBegK1;

        FAdr *piAdr1_ = piAdr1 + nOrb * iBlockBegK1;
        signed char *pSigns_ = pSigns + nOrb * iBlockBegK1;

        for (uint iKk = 0; iKk < nOrb * nBlkK; ++iKk)
          pInpK[iKk] = pSigns_[iKk] * pCoeff_[St1 * piAdr1_[iKk]];

        Mxm(pOutK, 1, nOrb, p1OpMatrix, 1, nOrb, pInpK, 1, nOrb, nOrb, nOrb,
            nBlkK, false, Prefactor);

        for (uint iKk = 0; iKk < nOrb * nBlkK; ++iKk)
          pResult_[St1 * piAdr1_[iKk]] += pSigns_[iKk] * pOutK[iKk];
      }
      Mem.Free(pInpK);
    };
  }

  Mem_.Free(pBeginOfMemory);
};

// void Apply1eOp( FScalar *pResult, FScalar *pCoeff, double Prefactor,
//     FScalar const *p1OpMatrix,
//     FOrbStringAdrTable const &Adr1, FOrbStringAdrTable const &Adr2,
//     uint St1, uint St2, FMemoryStack &Mem_)
// {
//     uint
//         nOrb = Adr1.nOrb();
//     // note: I configurations are |I1,I2>, with one of 1,2
//     // being alpha and the other being beta spin.
//
//     FMemoryStackArray MemStacks(Mem_);
//
//     // iterate through |I1> configurations (1 = alpha or beta)
//     #pragma omp parallel for schedule(dynamic)
//     for ( FAdr iStr1 = 0; iStr1 < Adr1.nStr(); ++ iStr1 ) {
//         FMemoryStack &Mem = MemStacks.GetStackOfThread();
//         FOrbPat
//             I1 = Adr1.MakePattern(iStr1);
//         FSubstResult
//             *pSubst;
//         uint
//             nSubst;
//         FormStringSubstsForSpin( pSubst, nSubst,
//             p1OpMatrix, Adr1, I1, Mem );
//         for ( uint iSubst = 0; iSubst < nSubst; ++ iSubst ) {
//             FSubstResult
//                 &s = pSubst[iSubst];
//             double
//                 tkl = p1OpMatrix[s.k + nOrb * s.l] * s.sign * Prefactor;
//             if ( tkl == 0 )
//                 continue;
//             FAdr
//                 iAdrC1 = St1 * s.iStr,
//                 iAdrR1 = St1 * iStr1;
//
//             // iterate through complementary I2 spin part of I configuration.
//             // The I2 part is unaffected by the operator we applied; it
//             // is thus just passed through from |c> to |r>.
//             for ( FAdr iStr2 = 0; iStr2 < Adr2.nStr(); ++ iStr2 )
//                 pResult[iAdrR1 + St2 * iStr2] += tkl * pCoeff[iAdrC1 + St2 *
//                 iStr2];
//                 // ^- beware of lots of memory jumping.
//         };
//         Mem.Free(pSubst);
//     };
// };

void FHamiltonianTerm1e::Contract(FFciVector &r, FFciVector &c,
                                  double Prefactor, FMemoryStack &Mem) {
  assert(compatible(r, c));

  // perform replacement on alpha strings, leaving beta invariant
  Apply1eOp(&r[0], &c[0], BaseFactor * Prefactor, pOp1MatrixA, r.AdrA, r.AdrB,
            1, r.nStrA, Mem);
  // perform replacement on beta strings, leaving alpha invariant
  Apply1eOp(&r[0], &c[0], BaseFactor * Prefactor, pOp1MatrixB, r.AdrB, r.AdrA,
            r.nStrA, 1, Mem);
};

void Add1RdmForSpin(FScalar *pRdm, FScalar *pCoeffL, FScalar *pCoeffR,
                    FOrbStringAdrTable const &Adr1,
                    FOrbStringAdrTable const &Adr2, uint St1, uint St2,
                    FMemoryStack &Mem) {
  // Operator acts on Adr1, Adr2 is left invariant. see Apply1eOp().
  uint nOrb = Adr1.nOrb();

  // iterate through |I1> configurations (1 = alpha or beta)
  for (FAdr iStr1 = 0; iStr1 < Adr1.nStr(); ++iStr1) {
    FOrbPat I1 = Adr1.MakePattern(iStr1);
    FSubstResult *pSubst;
    uint nSubst;
    FormStringSubstsForSpin(pSubst, nSubst, 0, Adr1, I1, Mem);
    for (uint iSubst = 0; iSubst < nSubst; ++iSubst) {
      FSubstResult &s = pSubst[iSubst];
      FAdr iAdrR1 = St1 * s.iStr, iAdrC1 = St1 * iStr1;
      double tkl = 0;

      for (FAdr iStr2 = 0; iStr2 < Adr2.nStr(); ++iStr2)
        tkl += pCoeffL[iAdrR1 + St2 * iStr2] * pCoeffR[iAdrC1 + St2 * iStr2];

      pRdm[s.k + nOrb * s.l] += s.sign * tkl;
    };
    Mem.Free(pSubst);
  };
};

// pDataK:  nPairs x [nBlk1 x nBlk2] array.
// performs either
//      DataK[kl, (K1,K2)] += <K1|c^k\sigma c_l\sigma|J1> Coeff[J1,K2]
//      (+symmetrization over kl on rhs)
// or
//      Resid[I1,K2] += <K1|c^l\sigma c_k\sigma|I1> DataK[kl, (K1,K2)]
// for a nBlk1 x nBlk2 block of strings in K, depending on 'Direction' (c or R).
void BlockContractCc1(double *pDataK, FStrInfo *pInfo1, uint nBlk1, uint nBlk2,
                      uint64_t nStK1, uint64_t nStK2,
                      FOrbStringAdrTable const &Adr1,
                      FOrbStringAdrTable const &Adr2, double *pCoeffs,
                      uint64_t nStC1, uint64_t nStC2, char Direction,
                      double *pCSum, double Prefactor) {
  assert(Direction == 'c' || Direction == 'R');

  for (uint iBlk1 = 0; iBlk1 < nBlk1; ++iBlk1) {
    FStrInfo &Info1 = pInfo1[iBlk1];

    // loop over k,l
    for (uint iSubst = 0; iSubst < Info1.nSubst; ++iSubst) {
      FSubstResult &s = Info1.pSubst[iSubst];
      double pf = s.sign * Prefactor;
      uint k_ = std::max(s.k, s.l), l_ = std::min(s.k, s.l),
           kl = (k_ * (k_ + 1)) / 2 + l_;
      double *pC1 = &pCoeffs[s.iStr * nStC1],
             *pK1 = &pDataK[kl + nStK1 * iBlk1];

      if (Direction == 'c') {
        for (uint iBlk2 = 0; iBlk2 < nBlk2; ++iBlk2) {
          double t = pf * pC1[iBlk2 * nStC2];
          pK1[iBlk2 * nStK2] += t;
          *pCSum += t * t;
        }
      } else {
        for (uint iBlk2 = 0; iBlk2 < nBlk2; ++iBlk2)
          pC1[iBlk2 * nStC2] += pf * pK1[iBlk2 * nStK2];
      }
    };
  }
};

// pDataK:  (nOrb x nOrb) x [nBlk1 x nBlk2] array.
// performs
//      DataK[kl, (K1,K2)] += <K1|c^k\sigma c_l\sigma|J1> Coeff[J1,K2]
//      (*without* symmetrization over kl on rhs)
// for a nBlk1 x nBlk2 block of strings in K..
void BlockContractCc1_NoSym(double *pDataK, uint iStOrbK, uint iStOrbL,
                            FStrInfo *pInfo1, uint nBlk1, uint nBlk2,
                            uint64_t nStK1, uint64_t nStK2,
                            FOrbStringAdrTable const &Adr1,
                            FOrbStringAdrTable const &Adr2, double *pCoeffs,
                            uint64_t nStC1, uint64_t nStC2, double *pCSum,
                            double Prefactor) {
  for (uint iBlk1 = 0; iBlk1 < nBlk1; ++iBlk1) {
    FStrInfo &Info1 = pInfo1[iBlk1];

    // loop over k,l
    for (uint iSubst = 0; iSubst < Info1.nSubst; ++iSubst) {
      FSubstResult &s = Info1.pSubst[iSubst];
      double pf = s.sign * Prefactor;
      double *pC1 = &pCoeffs[s.iStr * nStC1],
             *pK1 = &pDataK[iStOrbK * s.k + iStOrbL * s.l + nStK1 * iBlk1];

      for (uint iBlk2 = 0; iBlk2 < nBlk2; ++iBlk2) {
        double t = pf * pC1[iBlk2 * nStC2];
        pK1[iBlk2 * nStK2] += t;
        *pCSum += t * t;
      }
    };
  }
};

void SymmetrizeCiVector(FFciVector &x) {
  assert(x.nElecA == x.nElecB);
  //     double pf = (x.nElecA%2 == 0)? 1.0 : -1.0;
  double pf = 1.;
  for (FAdr iStrB = 0; iStrB < x.AdrB.nStr(); ++iStrB)
    for (FAdr iStrA = 0; iStrA < x.AdrA.nStr(); ++iStrA) {
      double f = (x(iStrA, iStrB) + pf * x(iStrB, iStrA)) / 2;
      x(iStrA, iStrB) = f;
      x(iStrB, iStrA) = pf * f;
    }
};

void FHamiltonianTerm2e::Contract(FFciVector &r, FFciVector &c,
                                  double Prefactor, FMemoryStack &Mem_) {
  assert(compatible(r, c));
  // simplest way of getting in the Prefactor at this point...
  Scale(r.pData, 1.0 / (BaseFactor * Prefactor), r.nData());

  // we do:
  //   r_I += [<I|c^r\alpha c_s\alpha|K> + <I|c^r\beta c_s\beta|K>] *
  //           (rs|tu) *
  //          [<K|c^t\alpha c_u\alpha|J> + <K|c^t\beta c_u\beta|J>] c_J *
  //
  FAdr nTgtBlkK = 32, nTgtBlkKB = 32;
  uint nOrb = r.nOrb, nPairs = nOrb * (nOrb + 1) / 2;
  FOrbStringAdrTable &AdrA = r.AdrA, &AdrB = r.AdrB;

  bool
      // singlet case: alpha/beta strings are equal.
      UseSymmetryAB = Spatial && (r.nSpin() == 0) && (r.IsSpinProjected);
  if (UseSymmetryAB)
    SymmetrizeCiVector(c);
  // ^- might be slightly non-symmetric (p-space and stuff).
  //    in that case numerical problems!!

  FMemoryStackArray MemStacks(Mem_);

// loop over blocks of alpha/beta strings of K.
#pragma omp parallel for schedule(dynamic)
  for (FAdr iBlockBegKB = 0; iBlockBegKB < AdrB.nStr();
       iBlockBegKB += nTgtBlkKB) {
    FMemoryStack &Mem = MemStacks.GetStackOfThread();

    FAdr iBlockEndKB = std::min(iBlockBegKB + nTgtBlkKB, AdrB.nStr());
    uint nBlkB = iBlockEndKB - iBlockBegKB; // number of beta-strings in block
    FStrInfo *pInfoA, *pInfoB;

    // find beta determinant strings on |I> and |J> side
    // (both are equal) which connect to Ks within the current beta block.
    Mem.Alloc(pInfoB, nBlkB);
    for (uint i = 0; i < nBlkB; ++i) {
      pInfoB[i].Str = AdrB.MakePattern(iBlockBegKB + i);
      FormStringSubstsForSpin(pInfoB[i].pSubst, pInfoB[i].nSubst, 0, AdrB,
                              pInfoB[i].Str, Mem);
    }

    for (FAdr iBlockBegKA = 0; iBlockBegKA < AdrA.nStr();
         iBlockBegKA += nTgtBlkK) {
      double Scale = BaseFactor * Prefactor;
      if (UseSymmetryAB && iBlockBegKA > iBlockBegKB)
        continue;
      if (UseSymmetryAB && iBlockBegKA < iBlockBegKB)
        Scale *= 2.0;

      FAdr iBlockEndKA = std::min(iBlockBegKA + nTgtBlkK, AdrA.nStr());
      uint nBlkA =
          iBlockEndKA - iBlockBegKA; // number of alpha-strings in block

      // find alpha determinant strings on |I>/|J> connecting to K.
      Mem.Alloc(pInfoA, nBlkA);
      for (uint i = 0; i < nBlkA; ++i) {
        pInfoA[i].Str = AdrA.MakePattern(iBlockBegKA + i);
        FormStringSubstsForSpin(pInfoA[i].pSubst, pInfoA[i].nSubst, 0, AdrA,
                                pInfoA[i].Str, Mem);
      }

      // get memory for storage of intermediate Dat[kl,K] block
      // and it's mxm result Dat[mn,K] = (mn|kl) Dat[kl,K].
      double *pInpA, *pOutA, *pInpB, *pOutB;
      uint const nIntKmn = nPairs * nBlkA * nBlkB * (Spatial ? 1 : 2);
      Mem.ClearAlloc(pInpA, nIntKmn);
      Mem.Alloc(pOutA, nIntKmn);
      if (Spatial) {
        pInpB = pInpA;
        pOutB = pOutA;
      } else {
        pInpB = pInpA + nIntKmn / 2;
        pOutB = pOutA + nIntKmn / 2;
      }

      double csum = 0;

      // form Inp[kl,K\alpha K\beta] = <K\alpha|c^k\alpha c_l\alpha|J\alpha> *
      // c[J\alpha,K\beta] (note: includes symmetrization)
      BlockContractCc1(pInpA, pInfoA, nBlkA, nBlkB, 1 * nPairs, nPairs * nBlkA,
                       AdrA, AdrB, &c(0, iBlockBegKB), 1, c.nStrA, 'c', &csum);
      // add Inp[kl,K\alpha K\beta] += <K\beta|c^k\beta c_l\beta|J\beta> *
      // c[K\alpha,J\beta]
      BlockContractCc1(pInpB, pInfoB, nBlkB, nBlkA, nBlkA * nPairs, 1 * nPairs,
                       AdrB, AdrA, &c(iBlockBegKA, 0), c.nStrA, 1, 'c', &csum);

      if (csum > ThrNegl) { // <- mainly helpful for first iterations
                            //      where c is nearly empty.
        // contract with integrals:
        //     Dat[mn,K] = (mn|kl) Dat[kl,K].
        Mxm(pOutA, 1, nPairs, pOp2MatrixAA, 1, nPairs, pInpA, 1, nPairs, nPairs,
            nPairs, nBlkA * nBlkB);
        if (!Spatial) {
          Mxm(pOutB, 1, nPairs, pOp2MatrixBB, 1, nPairs, pInpB, 1, nPairs,
              nPairs, nPairs, nBlkA * nBlkB);
          Mxm(pOutA, 1, nPairs, pOp2MatrixAB, 1, nPairs, pInpB, 1, nPairs,
              nPairs, nPairs, nBlkA * nBlkB, true);
          Mxm(pOutB, 1, nPairs, pOp2MatrixAB, nPairs, 1, pInpA, 1, nPairs,
              nPairs, nPairs, nBlkA * nBlkB, true);
        }

#pragma omp critical
        {
          // contract with coupling coefficients for writing
          // back to residuals.
          BlockContractCc1(pOutA, pInfoA, nBlkA, nBlkB, 1 * nPairs,
                           nPairs * nBlkA, AdrA, AdrB, &r(0, iBlockBegKB), 1,
                           c.nStrA, 'R', 0, Scale);
          // add Inp[kl,K\alpha K\beta] += <K\beta|c^k\beta c_l\beta|J\beta> *
          // c[K\alpha,J\beta]
          BlockContractCc1(pOutB, pInfoB, nBlkB, nBlkA, nBlkA * nPairs,
                           1 * nPairs, AdrB, AdrA, &r(iBlockBegKA, 0), c.nStrA,
                           1, 'R', 0, Scale);
        }
      }
      Mem.Free(pInfoA);
    };
    Mem.Free(pInfoB);
  }

  // symmetrize r. the k-block restriction above has broken A/B symmetry.
  if (UseSymmetryAB)
    SymmetrizeCiVector(r);
};

void FHamiltonianTerm2eHub::Contract(FFciVector &r, FFciVector &c,
                                     double Prefactor, FMemoryStack &Mem_) {
  assert(compatible(r, c));
  uint nOrb = r.nOrb;
  while (nOrb != 0 && pHubU[nOrb - 1] == 0)
    --nOrb;
  FOrbStringAdrTable &AdrA = r.AdrA, &AdrB = r.AdrB;
  double pf = BaseFactor * Prefactor;
#pragma omp parallel for schedule(dynamic)
  for (FAdr iStrB = 0; iStrB < AdrB.nStr(); ++iStrB) {
    FOrbPat StrB = AdrB.MakePattern(iStrB);
    for (FAdr iStrA = 0; iStrA < AdrA.nStr(); ++iStrA) {
      FOrbPat StrA = AdrA.MakePattern(iStrA);
      double f = 0;
      for (uint i = 0; i < nOrb; ++i) {
        int na = ((StrA >> i) & 1);
        int nb = ((StrB >> i) & 1);
        f += pHubU[i] * ((na & na) + 2 * (na & nb) + (nb & nb));
      }
      FAdr iStrAB = iStrA + AdrA.nStr() * iStrB;
      r.pData[iStrAB] += f * pf * c.pData[iStrAB];
    };
  }
};

// FFciData *gx = 0;

// void FHamiltonianTerm2e::Contract( FFciVector &r, FFciVector &c, double
// Prefactor, FMemoryStack &Mem )
// {
//     assert(compatible(r,c));
//     // simplest way of getting in the Prefactor at this point...
//     Scale(r.pData, 1.0/Prefactor, r.nData());
//
//
//     // we do:
//     //   r_I += [<I|c^r\alpha c_s\alpha|K> + <I|c^r\beta c_s\beta|K>] *
//     //           (rs|tu) *
//     //          [<K|c^t\alpha c_u\alpha|J> + <K|c^t\beta c_u\beta|J>] c_J *
//     //
//     uint
//         nOrb = r.nOrb,
//         nPairs = nOrb*(nOrb+1)/2;
//     FOrbStringAdrTable
//         &AdrA = r.AdrA,
//         &AdrB = r.AdrB;
//
//     // loop over K determinants
//     for ( FAdr iStrKA = 0; iStrKA < AdrA.nStr(); ++ iStrKA )
//     for ( FAdr iStrKB = 0; iStrKB < AdrB.nStr(); ++ iStrKB )
//     {
//         FOrbPat StrKA = AdrA.MakePattern(iStrKA);
//         FOrbPat StrKB = AdrA.MakePattern(iStrKB);
//         FSubstResult *pSubstA, *pSubstB;
//         uint nSubstA, nSubstB;
//         FormStringSubstsForSpin(pSubstA, nSubstA, 0, AdrA, StrKA, Mem);
//         FormStringSubstsForSpin(pSubstB, nSubstB, 0, AdrB, StrKB, Mem);
//
//         for ( uint isaI = 0; isaI < nSubstA; ++ isaI ) { FSubstResult &saI =
//         pSubstA[isaI]; for ( uint isaJ = 0; isaJ < nSubstA; ++ isaJ ) {
//         FSubstResult &saJ = pSubstA[isaJ];
//             r(saI.iStr, iStrKB) += saI.sign * saJ.sign *
//             gx->GetInt2e(saI.k,saI.l,saJ.k,saJ.l) * c(saJ.iStr, iStrKB);
//         }}
//         for ( uint isbI = 0; isbI < nSubstB; ++ isbI ) { FSubstResult &sbI =
//         pSubstB[isbI]; for ( uint isbJ = 0; isbJ < nSubstB; ++ isbJ ) {
//         FSubstResult &sbJ = pSubstB[isbJ];
//             r(iStrKA, sbI.iStr) += sbI.sign * sbJ.sign *
//             gx->GetInt2e(sbI.k,sbI.l,sbJ.k,sbJ.l) * c(iStrKA, sbJ.iStr);
//         }}
//         for ( uint isaI = 0; isaI < nSubstA; ++ isaI ) { FSubstResult &saI =
//         pSubstA[isaI]; for ( uint isbJ = 0; isbJ < nSubstB; ++ isbJ ) {
//         FSubstResult &sbJ = pSubstB[isbJ];
//             r(saI.iStr, iStrKB) += saI.sign * sbJ.sign *
//             gx->GetInt2e(saI.k,saI.l,sbJ.k,sbJ.l) * c(iStrKA, sbJ.iStr);
//         }}
//         for ( uint isbI = 0; isbI < nSubstB; ++ isbI ) { FSubstResult &sbI =
//         pSubstB[isbI]; for ( uint isaJ = 0; isaJ < nSubstA; ++ isaJ ) {
//         FSubstResult &saJ = pSubstA[isaJ];
//             r(iStrKA, sbI.iStr) += sbI.sign * saJ.sign *
//             gx->GetInt2e(sbI.k,sbI.l,saJ.k,saJ.l) * c(saJ.iStr, iStrKB);
//         }}
//
//         Mem.Free(pSubstA);
//     };
//
//     Scale(r.pData, Prefactor, r.nData());
// };
//
//  ^- this one works... (slow of course), but this should allow for
//     an easier adjustment to sparsity.

// kate: indent-width 4;
