/* This program is free software. It comes without any warranty, to the extent
 * permitted by applicable law. You may use it, redistribute it and/or modify
 * it, in whole or in part, provided that you do so at your own risk and do not
 * hold the developers or copyright holders liable for any claim, damages, or
 * other liabilities arising in connection with the software.
 *
 * Developed by Gerald Knizia, 2010--2012.
 */

#ifndef FCIDUMP_H
#define FCIDUMP_H

#include "CxPodArray.h"
#include <iosfwd>
using ct::TArray;

enum FIntClass {
  // equal integrals for alpha and beta spin
  INTCLASS_Spatial = 0x01,
  // spatial orbitals are different for alpha and beta spin, but
  // there are no cross terms between e.g., <A|int1e|B>.
  INTCLASS_Unrestricted = 0x02
};

struct FMatrix1 : public TArray<double> {
  typedef TArray<double> FBase;
  uint nRows, nCols;
  enum FInitFlags { INIT_Clear = 0x01 };

  FMatrix1() : nRows(0), nCols(0) {}
  FMatrix1(uint nRows, uint nCols, uint Flags) { Init(nRows, nCols, Flags); }

  void Init(uint nRows_, uint nCols_, uint Flags) {
    nRows = nRows_;
    nCols = nCols_;
    FBase::resize(nRows * nCols);
    if (Flags & INIT_Clear)
      FBase::clear_data();
  }

  double &operator()(uint i, uint k) { return (*this)[i + nRows * k]; }
  double const &operator()(uint i, uint k) const {
    return (*this)[i + nRows * k];
  }
};

inline uint i2e(uint i, uint j, uint k, uint l, uint nPairs) {
  // slow convenience routine for accessing integrals in (ij|kl) format.
  uint ij = std::max(i, j) * (std::max(i, j) + 1) / 2 + std::min(i, j);
  uint kl = std::max(k, l) * (std::max(k, l) + 1) / 2 + std::min(k, l);
  return ij + nPairs * kl;
};

struct FInt2eData : public FMatrix1 {
  void Init(uint nPairs, uint Flags) { FMatrix1::Init(nPairs, nPairs, Flags); }

  // slow convenience routine for accessing integrals in (ij|kl) format.
  double &operator()(uint i, uint j, uint k, uint l) {
    return (*this)[i2e(i, j, k, l, nRows)];
  };
  double const &operator()(uint i, uint j, uint k, uint l) const {
    return (*this)[i2e(i, j, k, l, nRows)];
  };
  using FMatrix1::operator();
};

struct FInt1eData : public FMatrix1 {
  void Init(uint nOrb, uint Flags) { FMatrix1::Init(nOrb, nOrb, Flags); }
};

// Represents integral/problem data from a standard FCIDUMP file.
struct FFciData {
  uint nElec, nOrb, Ms2, iWfSym, nPairs;
  uint nSyOrb[8]; // number of orbtials per symmetry.

  FIntClass IntClass;

  bool
      // if true, 2e integrals are regarded as prefactors of
      // e^r_t e^s_u instead of prefactors of e^rs_tu (default
      // for chemical hamiltonian)
      C1_Integrals;

  TArray<uint> iSyOrb; // symmetry of orbtial #i

  // on integrals: _S means 'spatial'/'restricted' (same integrals for
  // both spin cases). _A means for alpha spin, _B means for beta spin.
  // Which integrals are present is determined by this->IntClass
  FInt1eData
      // core hamiltonian/core fock matrix
      CoreH,
      CoreH_A, CoreH_B;
  TArray<double>
      // if given, used to form the diagonalization basis for --basis==Fock
      FockOp_A, FockOp_B;
  FInt2eData
      // (ij|kl) matrix of two-electron integrals. Format nPairs x nPairs
      Int2e,
      // integrals for (AA|AA), (BB|BB), and (AA|BB) spin-orbitals.
      Int2e_AA, Int2e_BB, Int2e_AB;
  FInt2eData
      // used for calculating p-trace energy if given (hacky).
      // otherwise copy of int2e with elements set to zero or 1/2 according to
      // weight.
      Int2e_ImpPrj,
      Int2e_ImpPrj_AA, Int2e_ImpPrj_BB, Int2e_ImpPrj_AB;
  double fCoreEnergy;

  bool HubbardSys;
  TArray<double> HubU; // (ii|ii) if and only if system is a Hubbard system.
  double
      // these point into CoreH_S or CoreH_A/_B, depending on context.
      *pCoreH_A,
      *pCoreH_B,
      // these point into Int2e_SS or Int2e_**, depending on context.
      *pInt2e_AA, *pInt2e_BB, *pInt2e_AB;

  void ReadDump(std::string const &FileName);
  void CheckSanity();
  void Finalize();
  void ClearAll();
};

void PrintMatrixGen(std::ostream &out, double const *pData, uint nRows,
                    uint nRowStride, uint nCols, uint nColStride,
                    std::string const &Name);

void WriteMatrixToFile(std::string const &FileName, std::string const Desc,
                       double *pData, uint nRows, uint nCols);

#endif // FCIDUMP_H

// kate: indent-width 4
