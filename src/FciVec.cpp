/* This program is free software. It comes without any warranty, to the extent
 * permitted by applicable law. You may use it, redistribute it and/or modify
 * it, in whole or in part, provided that you do so at your own risk and do not
 * hold the developers or copyright holders liable for any claim, damages, or
 * other liabilities arising in connection with the software.
 *
 * Developed by Gerald Knizia, 2010--2012.
 */

#include "FciCompat.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <ostream>
#include <stdexcept>
using boost::format;

#include "Fci.h"

void FOrbStringAdrTable::Init(uint nElec_, uint nOrb_) {
  m_nElec = nElec_;
  m_nOrb = nOrb_;

  m_AdrCount = SymDof(m_nOrb, m_nElec);

  MakeStringTable();
};

void FOrbStringAdrTable::MakeStringTable() {
  m_StrTable.reserve(m_AdrCount);
  AddStringsToTableR(0, m_nElec, 0);
  std::sort(m_StrTable.begin(), m_StrTable.end());
};

// recursively add patterns which add up to nElecLeft electons in the nOrbLeft
// last orbitals.
void FOrbStringAdrTable::AddStringsToTableR(FOrbPat OldPat, int nElecLeft,
                                            uint iFirstOrb) {
  if (nElecLeft == 0 && iFirstOrb <= m_nOrb) {
    m_StrTable.push_back(OldPat);
    return;
  };

  for (uint iOrb = iFirstOrb; iOrb < m_nOrb; ++iOrb) {
    FOrbPat NewPat = OldPat | (FOrbPat(1) << iOrb);
    AddStringsToTableR(NewPat, nElecLeft - 1, iOrb + 1);
  };
}

FAdr FOrbStringAdrTable::operator()(FOrbPat BitString) const {
  FAdr adr = std::lower_bound(m_StrTable.begin(), m_StrTable.end(), BitString) -
             m_StrTable.begin();
  assert(m_StrTable[adr] == BitString);
  return adr;
};

// return 0 if number of bits /SET/ in bits# [0..iPos) is even,
// 1 if it is odd
uint StringParityBeforePos(FOrbPat pat, uint iPos) {
  FOrbPat
      // mask out bits at iPos and above.
      tmp = pat & ((FOrbPat(1) << iPos) - 1);

  // count bit parity.
  tmp ^= (tmp >> 32ul); // only last 16 binary digits used.
  tmp ^= (tmp >> 16ul); // only last 8 binary digits used.
  tmp ^= (tmp >> 8ul);  // ...
  tmp ^= (tmp >> 4ul);
  //  note: binary encoded lookup table (i.e.,
  //  return (#magic-number >> (tmp & 4)) & 1;) instead of
  // the following two lines would be faster here.
  tmp ^= (tmp >> 2ul);
  return (tmp ^ (tmp >> 1ul)) & 1;
}

FFciVector::FFciVector(uint nElec_, uint nOrb_, int nSpin_,
                       bool IsSpinProjected_)
    : nElecA((nElec_ + nSpin_) / 2), nElecB((nElec_ - nSpin_) / 2), nOrb(nOrb_),
      AdrA(nElecA, nOrb_), AdrB(nElecB, nOrb_),
      IsSpinProjected(IsSpinProjected_) {
  nStrA = AdrA.nStr();
  nStrB = AdrB.nStr();
  pData = reinterpret_cast<FScalar *>(::malloc(nStrA * nStrB * sizeof(*pData)));
  Clear();
};

FFciVector::~FFciVector() { ::free(pData); };

void FFciVector::operator=(FFciVector const &other) {
  ::free(pData);
  nElecA = other.nElecA;
  nElecB = other.nElecB;
  nOrb = other.nOrb;
  nStrA = other.nStrA;
  nStrB = other.nStrB;
  AdrA = other.AdrA;
  AdrB = other.AdrB;
  IsSpinProjected = other.IsSpinProjected;

  pData = reinterpret_cast<FScalar *>(::malloc(nStrA * nStrB * sizeof(*pData)));
  ::memcpy(pData, other.pData, nStrA * nStrB * sizeof(*pData));
};

void FFciVector::Init0() {
  nElecA = 0;
  nElecB = 0;
  nOrb = 0;
  nStrA = 0;
  nStrB = 0;
  pData = 0;
};

FFciVector::FFciVector(FFciVector const &other) {
  Init0();
  *this = other;
};

void FFciVector::Clear() {
  ::memset(pData, 0, sizeof(*pData) * nStrA * nStrB);
};

double Dot(FFciVector const &a, FFciVector const &b) {
  assert(compatible(a, b));
  return Dot(a.pData, b.pData, a.nData());
};

void Add(FFciVector &r, FFciVector &x, double f) {
  assert(compatible(r, x));
  Add(r.pData, x.pData, f, r.nData());
};

double FFciVector::Norm() const { return std::sqrt(Dot(*this, *this)); };

void FFciVector::Normalize() {
  double f = Norm();
  if (f == 0.0)
    throw std::runtime_error(
        "attempted to normalize a non-normalizable vector.");
  Scale(pData, 1. / f, nStrA * nStrB);
};

// project out components which lie not within the spin-space of TargetS.
void FFciVector::ProjectSpin(FMemoryStack &Mem_, uint TargetS) {
  if (TargetS == 0xffff)
    TargetS = nSpin();
  uint nElec = nElecA + nElecB, Ms2 = nElecA - nElecB;
  assert(Ms2 <= TargetS);
  assert(Ms2 % 2 == TargetS % 2);

  // how this works:
  //    - we loop over all spatial orbital configurations (that is: assign
  //      'doubly occupied', 'singly occupied', 'empty' to each orbital in all
  //      possible ways).
  //    - for each orbital configuration, we get all determinants which have
  //      this orbital configuration. All such determinants share the same
  //      closed-shell orbital part, but differ in their assignment of
  //      alpha/beta spin to the open shells.
  //    - the determinants for a spatial orbital configuration are projected
  //      onto the CSFs with TargetS/M for this orbital configuration. And
  //      then back to determinants
  // Since we actually need to access all determinants, and in rather
  // unfavorable access patterns, this entire thing is rather slow,
  // unfortunatelly. If someone has a clever idea on how to do this better,
  // please drop me a note (cgk at gmx.de).
  for (uint nClosed = 0; nClosed <= (nElec - Ms2) / 2; ++nClosed) {
    void *pBeginOfMemory = Mem_.Alloc(0);
    uint nOpen = nElec - 2 * nClosed;
    if (nOpen > nOrb || nClosed > nOrb)
      continue;

    uint nDets, nCsfs;
    double
        // nDets x nCsfs matrix denoting the expansion coefficients of
        // CSFs in determinants of projection M in the open-shell orbitals.
        *pCsfCoeffs,
        // nDets x nDets matrix denoting the projector onto the CSF space
        // of S/M (i.e., CsfCoeffs * CsfCoeffs^T)
        *pCsfProjector;
    FOrbPat *pAlphaBeg;
    // pAlphaBeg: valid distributions of alpha electrons amongst the open
    // orbitals (all open orbitals which are not alpha are beta)
    MakeCsfExpansion(pCsfCoeffs, pAlphaBeg, nDets, nCsfs, TargetS, Ms2, nOpen,
                     true, Mem_);

    Mem_.Alloc(pCsfProjector, nDets * nDets);
    if (nCsfs != 0)
      Mxm(pCsfProjector, 1, nDets, pCsfCoeffs, 1, nDets, pCsfCoeffs, nDets, 1,
          nDets, nCsfs, nDets);
    // ^- could use Syrk for that.
    else {
      assert((nDets == 0) == (TargetS == 0));
      if (nDets == 0) {
        // there are no open CSFs because there are no open shells.
        // Leave the closed shells alone.
        Mem_.Free(pBeginOfMemory);
        continue;
      } else {
        // combination of spins was not valid. wipe out the configurations.
        memset(pCsfProjector, 0, sizeof(*pCsfProjector) * nDets * nDets);
      }
    }

    // list all patterns of nClosed doubly occupied orbitals
    // in nOrb orbitals and nOpen singly occupied orbitals.
    FOrbPat *pClosedBeg, *pClosedEnd, *pOpenBeg, *pOpenEnd;

    ListBitPatterns(pClosedBeg, pClosedEnd, nClosed, nOrb, Mem_);
    ListBitPatterns(pOpenBeg, pOpenEnd, nOpen, nOrb, Mem_);

    //         FAdr
    //             nCsfTotal = 0,
    //             nDetTotal = 0;

    // loop through all valid orbital configuration combinations.
    {
      std::size_t nOpenPat = (std::size_t)(pOpenEnd - pOpenBeg);
      FMemoryStackArray MemStacks(Mem_);
#pragma omp parallel for schedule(dynamic)
      for (std::size_t iOpenPat = 0; iOpenPat < nOpenPat; ++iOpenPat) {
        FMemoryStack &Mem = MemStacks.GetStackOfThread();
        FOrbPat Open = pOpenBeg[iOpenPat];
        // translate determinant patterns in terms of generic nOpen orbital
        // patterns into orbital patterns for our concrete set of open-shell
        // orbitals.
        FOrbPat *pStrA;
        Mem.Alloc(pStrA, nDets);
        for (std::size_t iDet = 0; iDet < nDets; ++iDet) {
          uint iOpen = 0;
          FOrbPat StrA = 0, AlphaDet = pAlphaBeg[iDet];
          for (uint iOrb = 0; iOrb < nOrb; ++iOrb) {
            FOrbPat OpenBit = (FOrbPat(1) << iOrb);
            if (Open & OpenBit) {
              if (AlphaDet & (FOrbPat(1) << iOpen))
                StrA |= OpenBit;
              iOpen += 1;
            }
          }
          assert(iOpen == nOpen);
          assert(CountBits(StrA) == nElecA - nClosed);

          pStrA[iDet] = StrA;
        };

        for (FOrbPat *pClosed = pClosedBeg; pClosed != pClosedEnd; ++pClosed) {
          FOrbPat Closed = *pClosed;
          if ((Open & Closed) != 0)
            // open/closed patterns have at least one orbital in common.
            continue;
          //                     nCsfTotal += nCsfs;
          //                     nDetTotal += nDets;

          // find all determinants matching the spatial orbital configuration.
          double **pDetAdrs,
              *pCoeffs; // coefficient of this determinant in c.
          Mem.Alloc(pDetAdrs, nDets);
          Mem.Alloc(pCoeffs, 2 * nDets);
          for (uint iDet = 0; iDet < nDets; ++iDet) {
            FOrbPat StrA = Closed | pStrA[iDet],
                    StrB = Closed | (Open ^ pStrA[iDet]);
            assert(CountBits(StrA) == nElecA);
            assert(CountBits(StrB) == nElecB);

            pDetAdrs[iDet] = &pData[AdrA(StrA) + nStrA * AdrB(StrB)];
            pCoeffs[iDet] = *pDetAdrs[iDet];
          };

          // project onto subspace spanned by the CSFs.
          double *pCoeffProj = pCoeffs + nDets;
          Mxv(pCoeffProj, 1, pCsfProjector, nDets, 1, pCoeffs, 1, nDets, nDets);

          // copy back coefficients.
          for (uint iDet = 0; iDet < nDets; ++iDet)
            *pDetAdrs[iDet] = pCoeffProj[iDet];

          Mem.Free(pDetAdrs);
        }
        Mem.Free(pStrA);
      }
    }
    //         std::cout << format("nOpen = %2i ->  nConfCl = %5i   nConfOp =
    //         %5i  nCsf = %3i [%8i tot.]  nDet = %3i [%8i tot.]\n")
    //             % nOpen % (pClosedEnd-pClosedBeg) % (pOpenEnd - pOpenBeg) %
    //             nCsfs % nCsfTotal % nDets %  nDetTotal;
    Mem_.Free(pBeginOfMemory);
  };
};

#include <iostream> // FIXME: remove this.

// project out components not reachable by nSubstOrder single e- substitutions
// from the lowest determinant.
void FFciVector::ProjectToSubstOrder(int nSubstOrder, FMemoryStack &Mem) {
  FOrbPat RefPatA = FOrbPat(-1) >> (8 * sizeof(FOrbPat) - nElecA),
          RefPatB = FOrbPat(-1) >> (8 * sizeof(FOrbPat) - nElecB);
  //     std::cout << format("ref pattern A: %8x  (nElecA = %i)\n") % RefPatA %
  //     nElecA; std::cout << format("ref pattern B: %8x  (nElecB = %i)\n") %
  //     RefPatB % nElecB; std::cout << format("first det: A = %8x  B = %8x\n")
  //     % AdrA.MakePattern(0) % AdrB.MakePattern(0);
  for (FAdr iStrB = 0; iStrB < nStrB; ++iStrB)
    for (FAdr iStrA = 0; iStrA < nStrA; ++iStrA) {
      int
          // actual substitution level of the determinant.
          nDetSubst = 0;
      nDetSubst += CountBits(AdrA.MakePattern(iStrA) ^ RefPatA);
      nDetSubst += CountBits(AdrB.MakePattern(iStrB) ^ RefPatB);
      nDetSubst /= 2; // each substitutions changes two orbitals.
      if ((nSubstOrder > 0 && nDetSubst > nSubstOrder) ||
          (nSubstOrder < 0 && nDetSubst != -nSubstOrder &&
           !(iStrA == 0 && iStrB == 0)))
        (*this)(iStrA, iStrB) = 0;
    };
}

void FFciVector::Print(std::ostream &out, std::string const &VecName,
                       double CutWeight, int Level) const {
  out << format(" Dump of FCI vector %s [print threshold = %f]:") % VecName %
             CutWeight
      << std::endl;
  for (uint iStrB = 0; iStrB < nStrB; ++iStrB)
    for (uint iStrA = 0; iStrA < nStrA; ++iStrA) {
      double f = (*this)(iStrA, iStrB);
      if (std::abs(f) < CutWeight)
        continue;
      std::string DetDesc =
          FmtDet(AdrA.MakePattern(iStrA), AdrB.MakePattern(iStrB), nOrb);
      out << format("   %s[%s] = %f\n") % VecName % DetDesc % f;
    };
  out << std::endl;
};

bool compatible(FFciVector const &a, FFciVector const &b) {
  return a.nStrA == b.nStrA && a.nStrB == b.nStrB;
};

std::string FmtPat(FOrbPat pat, uint nMaxOrb) {
  std::string s(nMaxOrb, ' ');
  assert(pat < (FOrbPat(1) << (nMaxOrb + 1)));
  for (uint i = 0; i < nMaxOrb; ++i)
    s[i] = ((pat & (FOrbPat(1) << i)) != 0) ? '1' : '_';
  return s;
};

std::string FmtDet(FOrbPat patA, FOrbPat patB, uint nMaxOrb) {
  std::string s(nMaxOrb, ' ');
  assert(patA < (FOrbPat(1) << (nMaxOrb + 1)));
  assert(patB < (FOrbPat(1) << (nMaxOrb + 1)));
  for (uint i = 0; i < nMaxOrb; ++i) {
    bool occA = ((patA & (FOrbPat(1) << i)) != 0),
         occB = ((patB & (FOrbPat(1) << i)) != 0);
    char c = '.';
    if (occA && occB)
      c = '2';
    if (!occA && occB)
      c = 'b';
    if (occA && !occB)
      c = 'a';

    s[i] = c;
  }
  return s;
};

uint CountBits(uint64_t bpat) {
  typedef uint64_t T;
  // masks for counting the number of bits set in
  // a consecutive sequence of 1,2,4,8,16,32 bits
  T const mask0 = 0x5555555555555555ULL;
  T const mask1 = 0x3333333333333333ULL;
  T const mask2 = 0x0F0F0F0F0F0F0F0FULL;
  T const mask3 = 0x00FF00FF00FF00FFULL;
  T const mask4 = 0x0000FFFF0000FFFFULL;
  T const mask5 = 0x00000000FFFFFFFFULL;
  // position calculations
  T cnt = bpat;
  cnt = (cnt & mask0) + ((cnt >> 1) & mask0);
  cnt = (cnt & mask1) + ((cnt >> 2) & mask1);
  cnt = (cnt & mask2) + ((cnt >> 4) & mask2);
  cnt = (cnt & mask3) + ((cnt >> 8) & mask3);
  cnt = (cnt & mask4) + ((cnt >> 16) & mask4);
  cnt = (cnt & mask5) + ((cnt >> 32) & mask5);
  return static_cast<uint>(cnt);
}

uint CountBits(uint32_t bpat) {
  typedef uint32_t T;
  // masks for counting the number of bits set in
  // a consecutive sequence of 1,2,4,8,16,32 bits
  T const mask0 = 0x55555555U;
  T const mask1 = 0x33333333U;
  T const mask2 = 0x0F0F0F0FU;
  T const mask3 = 0x00FF00FFU;
  T const mask4 = 0x0000FFFFU;
  // position calculations
  T cnt = bpat;
  cnt = (cnt & mask0) + ((cnt >> 1) & mask0);
  cnt = (cnt & mask1) + ((cnt >> 2) & mask1);
  cnt = (cnt & mask2) + ((cnt >> 4) & mask2);
  cnt = (cnt & mask3) + ((cnt >> 8) & mask3);
  cnt = (cnt & mask4) + ((cnt >> 16) & mask4);
  return static_cast<uint>(cnt);
}

// kate: indent-mode normal; indent-width 4;
