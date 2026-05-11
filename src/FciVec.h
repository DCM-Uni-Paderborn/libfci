/* This program is free software. It comes without any warranty, to the extent
 * permitted by applicable law. You may use it, redistribute it and/or modify
 * it, in whole or in part, provided that you do so at your own risk and do not
 * hold the developers or copyright holders liable for any claim, damages, or
 * other liabilities arising in connection with the software.
 *
 * Developed by Gerald Knizia, 2010--2012.
 */

#ifndef FCIVEC_H
#define FCIVEC_H

#include "FciCompat.h"
#include <string>
#include <vector>

typedef uint64_t
    // encodes FCI spin orbital occupation patterns. Iff bit #i is set,
    // spin-orbital i is occupied with an electron.
    FOrbPat;
typedef std::size_t
    // used for addressing occupation patterns and vector elements.
    FAdr;

std::string FmtPat(FOrbPat pat, uint nMaxOrb);
std::string FmtDet(FOrbPat patA, FOrbPat patB, uint nMaxOrb);

uint64_t binomial_coefficient(uint N, uint k);

// number of symmetric degrees of freedom for dimension nDim.
inline uint64_t SymDof(uint N, uint nDim) {
  uint64_t Res = 1, Den = 1;
  for (uint i = 0; i < nDim; ++i) {
    Res *= (N - i);
    Den *= (i + 1);
  }
  return Res / Den;
};

// auxiliary object for addressing orbital occupation patterns (for one spin).
// Provides addressing of all patterns with a fixed total number of electrons
// in a given number of orbitals.
struct FOrbStringAdrTable {
  FAdr operator()(FOrbPat BitString) const;

  FOrbPat MakePattern(FAdr adr) const {
    assert(adr < m_AdrCount);
    return m_StrTable[adr];
  };

  FOrbStringAdrTable(){};
  FOrbStringAdrTable(uint nElec, uint nOrb) { Init(nElec, nOrb); }

  void Init(uint nElec, uint nOrb);

  FAdr nStr() const { return m_AdrCount; }

  uint nOrb() const { return m_nOrb; }
  uint nElec() const { return m_nElec; }

private:
  //     std::vector<FAdr>
  //         m_AdrTable;
  std::vector<FOrbPat> m_StrTable; // hmpf. should probably not store this.
  FAdr m_AdrCount;
  uint m_nElec, m_nOrb;

  void MakeStringTable();
  void AddStringsToTableR(FOrbPat OldPat, int nElecLeft, uint iFirstOrb);
};

// return 0 if number of bits /SET/ in bits# [0..iPos) is even,
// 1 if it is odd
uint StringParityBeforePos(FOrbPat pat, uint iPos);

typedef double FScalar;

// a FCI vector. The data is stored as a matrix
//    M[iAdrA, iAdrB]
// where iAdrA and iAdrB are indices of orbital occupation strings of alpha/beta
// electrons (see FOrbStringAdrTable). The number of electrons for alpha and
// beta electrons may be different, but the number of orbitals may not. The
// total wave function represented is thus best interpreted as:
//    |I> = AntiSym[|Ialpha> * |Ibeta>]
// where |Ialpha> = c^{i1 alpha} c^{i2 alpha} ... c^{iNA alpha} |0>
// and   |Ibeta> = c^{i1 beta} c^{i2 beta} ... c^{iNA beta} |0>
// are the determinants of alpha/beta electrons (first all alpha, then all
// beta).
struct FFciVector : public FIntrusivePtrDest {
  FFciVector() { Init0(); };
  FFciVector(uint nElec, uint nOrb, int nSpin, bool SpinProjected);
  ~FFciVector();

  uint nElecA, nElecB, nOrb;
  FOrbStringAdrTable AdrA, AdrB;
  uint nStrA, nStrB;
  FScalar *pData;
  bool
      // set if this vector is supposed to be projected onto the
      // spin-space with minimal S for the given Ms2.
      IsSpinProjected;

  std::size_t nData() const { return nStrA * nStrB; };

  uint nElec() const { return nElecA + nElecB; }
  int nSpin() const { return nElecA - nElecB; }

  void Clear();
  double Norm() const;
  void Normalize();

  // project out spin components which lie not within the spin-space of TargetS.
  // If TargetS is not given, it is assumed that TargetS is nElecA - nElecB.
  void ProjectSpin(FMemoryStack &Mem, uint TargetS = 0xffff);
  void ProjectSpinPartial(FFciVector &Out, FMemoryStack &Mem, uint TargetS,
                          uint nPartialOrb);
  // project out determinants not reachable by nSubstOrder single-electron
  // substitutions from the first determinant.
  // nSubstOrder > 0: by ``at most'' nSubstOrder subtitutions.
  // bSubstOrder < 0: by 0 or exactly |nSubstOrder| substitutions.
  void ProjectToSubstOrder(int nSubstOrder, FMemoryStack &Mem);

  void Print(std::ostream &out, std::string const &VecName,
             double CutWeight = 1e-10, int Level = 0) const;

  FScalar &operator()(uint64_t iStrA, uint64_t iStrB) {
    assert(iStrA < nStrA && iStrB < nStrB);
    return pData[iStrA + nStrA * iStrB];
  }
  FScalar operator()(uint64_t iStrA, uint64_t iStrB) const {
    assert(iStrA < nStrA && iStrB < nStrB);
    return pData[iStrA + nStrA * iStrB];
  }
  FScalar &operator[](uint64_t iAdr) {
    assert(iAdr < nStrA * nStrB);
    return pData[iAdr];
  }
  FScalar operator[](uint64_t iAdr) const {
    assert(iAdr < nStrA * nStrB);
    return pData[iAdr];
  }

  void operator=(FFciVector const &other);
  FFciVector(FFciVector const &other);

private:
  void Init0();
};

typedef boost::intrusive_ptr<FFciVector> FFciVectorPtr;

bool compatible(FFciVector const &a, FFciVector const &b);
double Dot(FFciVector const &a, FFciVector const &b);
void Add(FFciVector &r, FFciVector &x, double f);

uint CountBits(uint64_t bpat);
uint CountBits(uint32_t bpat);

#endif // FCIVEC_H

// kate: indent-mode normal; indent-width 4;
