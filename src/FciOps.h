/* This program is free software. It comes without any warranty, to the extent
 * permitted by applicable law. You may use it, redistribute it and/or modify
 * it, in whole or in part, provided that you do so at your own risk and do not
 * hold the developers or copyright holders liable for any claim, damages, or
 * other liabilities arising in connection with the software.
 *
 * Developed by Gerald Knizia, 2010--2012.
 */

#ifndef FCIOPS_H
#define FCIOPS_H

double const ThrNegl = 1e-16;

struct FHamiltonianTerm : public FIntrusivePtrDest {
  double
      // op will be multiplied by this factor (additionally to prefactor on
      // Contract)
      BaseFactor;

  // perform |r> += f * H |c> for some part of the Hamiltonian
  virtual void Contract(FFciVector &r, FFciVector &c, double Prefactor,
                        FMemoryStack &Mem) = 0;

  FHamiltonianTerm() : BaseFactor(1.0){};
  virtual ~FHamiltonianTerm();
};

typedef boost::intrusive_ptr<FHamiltonianTerm> FHamiltonianTermPtr;

// represents h_ij c^i c_j; this h_ij acts as identity on the spin space.
struct FHamiltonianTerm1e : public FHamiltonianTerm {
  void Contract(FFciVector &r, FFciVector &c, double Prefactor,
                FMemoryStack &Mem); // override

  FHamiltonianTerm1e(uint nOrb_, double const *pOp1MatrixA_,
                     double const *pOp1MatrixB_ = 0)
      : nOrb(nOrb_), Spatial(pOp1MatrixB_ == 0), pOp1MatrixA(pOp1MatrixA_),
        pOp1MatrixB(Spatial ? pOp1MatrixA_ : pOp1MatrixB_){};

  uint nOrb;
  bool Spatial;
  double const
      // p[i+nOrb*j] is the h_ij matrix element.
      *pOp1MatrixA,
      *pOp1MatrixB;
};

// represents [1/2] (ij|kl) E^i_j E^k_l.
//    WARNING:
//      o /NOT/ <ij|kl> E^ij_kl (see notes on ModCoreH!)
//      o beware of factor 1/2.
// Integrals for (AA|AA), (BB|BB), and (AA|BB) can be different, but that will
// result in slow execution.
struct FHamiltonianTerm2e : public FHamiltonianTerm {
  void Contract(FFciVector &r, FFciVector &c, double Prefactor,
                FMemoryStack &Mem); // override

  FHamiltonianTerm2e(uint nOrb_, double const *pOp2MatrixAA_,
                     double const *pOp2MatrixBB_ = 0,
                     double const *pOp2MatrixAB_ = 0)
      : nOrb(nOrb_), Spatial(pOp2MatrixBB_ == 0), pOp2MatrixAA(pOp2MatrixAA_),
        pOp2MatrixBB(Spatial ? pOp2MatrixAA_ : pOp2MatrixBB_),
        pOp2MatrixAB(Spatial ? pOp2MatrixAA_ : pOp2MatrixAB_) {
    BaseFactor = 0.5;
  };

  uint nOrb;
  bool Spatial;
  double const
      // [iPairs(i,j) + nPairs*iPair(k,l)] is the (ij|kl) matrix element.
      *pOp2MatrixAA,
      *pOp2MatrixBB, *pOp2MatrixAB;
};

struct FHamiltonianTerm2eHub : public FHamiltonianTerm {
  void Contract(FFciVector &r, FFciVector &c, double Prefactor,
                FMemoryStack &Mem); // override

  FHamiltonianTerm2eHub(uint nOrb_, double const *pHubU_)
      : nOrb(nOrb_), pHubU(pHubU_) {
    BaseFactor = 0.5;
  };

  uint nOrb;
  double const *pHubU;
};

// represents the result of an operation
//    <K[sigma]| c^{sigma k} c_{sigma l} |I[sigma]>
// for given |I[sigma]>, where I[sigma] and |K[sigma]> are
// orbital strings of the given spin projection (i.e., alpha or
// beta strings)
struct FSubstResult {
  char k, l; // k, l references
  signed char sign; // fermionic sign (+/- 1), also when plain char is unsigned
  FOrbPat
      Str; // orbital string K resulting from c^{sigma k} c_{sigma l} |I[sigma]>
  FAdr iStr; // address of string K
};

std::ostream &operator<<(std::ostream &out, FSubstResult const &o);

// form a sparse list of all |K> which can be reached by applying
// c^k_l on string |I>. These are supposed to have the same spin
// projection (e.g., c^{k\alpha}_{l\beta} |I\alpha>.
//
// data is created for all k,l in 0 .. Adr.nOrb-1;
// The nOrb x nOrb matrix p1OpMatrix, if supplied, is used to
// decide whether a specific kl entry is to be processed (if
// the matrix element is 0, it is not)
void FormStringSubstsForSpin(FSubstResult *&pResult, uint &nEntries,
                             FScalar const *p1OpMatrix,
                             FOrbStringAdrTable const &Adr, FOrbPat I,
                             FMemoryStack &Mem);

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
               FMemoryStack &Mem_);

void Add1RdmForSpin(FScalar *pRdm, FScalar *pCoeffL, FScalar *pCoeffR,
                    FOrbStringAdrTable const &Adr1,
                    FOrbStringAdrTable const &Adr2, uint St1, uint St2,
                    FMemoryStack &Mem);

struct FStrInfo {
  FOrbPat Str;
  FSubstResult *pSubst;
  uint nSubst;
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
                      double *pCSum, double Prefactor = 1.0);

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
                            double Prefactor = 1.0);

// +1: positive parity, -1: negative partiy, 0: annihilated.
signed char ApplyCop1(FOrbPat &Out, FOrbPat In, uint iOrb,
                      int iCreateOrDestroy);

void SymmetrizeCiVector(FFciVector &x);

#endif // FCIOPS_H
