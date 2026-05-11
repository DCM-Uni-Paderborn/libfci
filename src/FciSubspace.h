/* This program is free software. It comes without any warranty, to the extent
 * permitted by applicable law. You may use it, redistribute it and/or modify
 * it, in whole or in part, provided that you do so at your own risk and do not
 * hold the developers or copyright holders liable for any claim, damages, or
 * other liabilities arising in connection with the software.
 *
 * Developed by Gerald Knizia, 2010--2012.
 */

#ifndef FCISUBSPACE_H
#define FCISUBSPACE_H

#include "CxPodArray.h"
#include "CxStorageDevice.h"

struct FFciOptions {
  uint DiisDimension, DiisBlockSize;
  std::string DiagBasis, FileName_PSpaceH;
  double ThrVar, ThrPrintC;
  uint nPSpace, iPSpaceRoot, nMaxIt;
  bool
      // if true, project initial guess and trial vectors onto subspaces
      // spanned by the total spin Ms2 (if not given, this just gives the
      // M quantum number; M=0, however, can still be a triplet, quintet,
      // ...)
      ProjectSpin,
      // if true, assume we are dealing with hard-core bosons instead of
      // fermions. Creation/annihilation operators will be assumed to commute
      // instead of anti-commute, and the final wave function is a permanent
      // (but still with determinant-like occupation restrictions: a spin-
      // orbital can be occupied at most once. That's the "hard-core" part).
      BosonicSigns,
      // if true, try absorbing 1e terms into 2e terms. Should only be
      // nonzero if using the Hamiltonian on other than N-electron
      // configuration
      TryAbsorb1e;
  int nEmulatedCiOrder;
  FFciOptions();
};

FAdr const AdrNotFound = static_cast<FAdr>(-1);

struct FHamiltonianData {
  double const *pInt2e_AA, *pInt2e_BB, *pInt2e_AB, // nPairs x nPairs matrices.
      *pInt1e_A, *pInt1e_B; // nOrb x nOrb matrices. May be 0.
  uint nPairs, nOrb;
  bool Absorb1e, Spatial;

  FHamiltonianData(){};
  FHamiltonianData(double *pInt2e_AA_, double *pInt2e_BB_, double *pInt2e_AB_,
                   uint nPairs_, double *pInt1e_A_, double *pInt1e_B_,
                   uint nOrb_, bool Absorb1e_, bool Spatial_);
};

struct FPSpace {
  typedef TArray<FAdr> FConfList;
  FConfList Confs; // determinant indices of p-space configurations.
  typedef TArray<FAdr> FAdrList;
  FAdrList
      // sorted lists of all alpha/beta strings which occur in at
      // least a single configuration.
      HeldAlphaStr,
      HeldBetaStr;
  uint nConfs() const { return Confs.size(); };
  typedef TArray<double> FScalarArray;
  FScalarArray H, // nConf x nConf p-space Hamiltonian
      Ev,         // nConf x nConf eigenvector matrix of Hamiltonian
      Ew;         // nConf eigenvectors
  FOrbStringAdrTable const *pAdrA, *pAdrB;
  FHamiltonianData Data;

  void Init(FConfList const &Confs_, FHamiltonianData const &Data_,
            FOrbStringAdrTable const &AdrA_, FOrbStringAdrTable const &AdrB_,
            FMemoryStack &Mem);
  void MakeH(FMemoryStack &Mem);
  void MakeHEvs();
  bool HaveStr(FAdr iStr, FAdrList const &StrList) const;
  FAdr FindConf(FAdr iConf) const;
  FAdr FindConf(FAdr iStrA, FAdr iStrB) const;
};

// keeps subspace data for davidson diagonalization
struct FSubspaceStates {
  typedef TArray<double> FScalarArray;
  FScalarArray H, S, Ew;
  uint nDim, nMaxDim, iNext, iThis, nDimUsed, nConfP, nDiag, nMaxDiag;
  size_t nAmpLen; // residual and amplitude length must be equal for davidson

  explicit FSubspaceStates(uint nMaxDim_, uint DiisBlockSize = 1 << 20);
  ~FSubspaceStates();

  // 1) add vectors |c> and |s> := H |c> to iterative subspace.
  // 2) diagonalize in iterative subspace and return new amplitude
  //    and sigma vectors.
  // (..so actually, what we refer here to as residual is not actually
  //    the residual, but the \sigma = H * c vector)
  void Apply(double *pThisAmp, double *pThisRes, size_t nLength,
             FPSpace &PSpace, FMemoryStack &Mem);

  void OrthogonalizeUpdate(double *pThisAmp, size_t nLength, FPSpace &PSpace,
                           FMemoryStack &Mem);

private:
  std::size_t BlockSize, nBlocks;
  bool NoIO; // set if everything fits into memory.

  // subspace residuals and amplitudes.
  // arrays of length BlockSize * nMaxDim
  double *pRess, *pAmps;

  ct::FStorageBlock
      // stores the subspace vectors. Layout is as follows:
      //      Amp[0*bz..bz,0..nMaxDim] Res[0*bz..1*bz,0..nMaxDim]
      //      Amp[1*bz..bz,0..nMaxDim] Res[1*bz..2*bz,0..nMaxDim]
      //      ...
      // where bz is the BlockSize member. Last entry may have less
      // data than BlockSize if the sizes don't add up to vector length.
      FileData;

  // reads amplitude/residual block iBlock from storage to pRess/pAmps.
  // returns actual block size (may be be less than BlockSize for the last
  // block).
  void ReadBlock(std::size_t &iOff, std::size_t &iBeg, std::size_t &nSize,
                 uint iBlock);
  // writes amplitude/residual block iBlock to pRess/pAmps.
  void WriteBlock(uint iBlock, uint iFirst = 0, uint iLast = 0xffffffff);

  void GetBlockBoundaries(std::size_t &iOff, std::size_t &iBeg,
                          std::size_t &nSize, uint iBlock);
};

#endif // FCISUBSPACE_H

// kate: indent-width 4
