/* This program is free software. It comes without any warranty, to the extent
 * permitted by applicable law. You may use it, redistribute it and/or modify
 * it, in whole or in part, provided that you do so at your own risk and do not
 * hold the developers or copyright holders liable for any claim, damages, or
 * other liabilities arising in connection with the software.
 *
 * Developed by Gerald Knizia, 2010--2012.
 */

#include "FciCompat.h"
#include <cctype> // for isalnum / isspace
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
using boost::format;
using std::cout;
using std::endl;
#include "Fci.h"

/// prints a general rectangular matrix specified by memory layout.
/// Element m(r,c) is given by pData[ r * nRowStride + c * nColStride ].
void PrintMatrixGen(std::ostream &out, double const *pData, uint nRows,
                    uint nRowStride, uint nCols, uint nColStride,
                    std::string const &Name) {
  uint nFloatWidth = 14, nFloatPrec = 8;
  std::stringstream ColFmt, FloatFmt;
  ColFmt << "%" << nFloatWidth - 3 << "i   ";
  FloatFmt << "%" << nFloatWidth << "." << nFloatPrec << "f";

  if (Name != "") {
    out << "  Matrix " << Name << ", " << nRows << "x" << nCols << "."
        << std::endl;
    if (nRows * nCols == 0)
      return;
  }
  out << "           ";
  for (uint nCol = 0; nCol < nCols; ++nCol)
    out << format(ColFmt.str()) % nCol;
  out << "\n";
  for (uint nRow = 0; nRow < nRows; ++nRow) {
    out << "    " << format("%4i") % nRow << "   ";
    for (uint nCol = 0; nCol < nCols; ++nCol) {
      double const &f = pData[nRow * nRowStride + nCol * nColStride];
      out << format(FloatFmt.str()) % f;
    };
    out << "\n";
  }
  out.flush();
};

void WriteMatrixToFile(std::string const &FileName, std::string const Desc,
                       double *pData, uint nRows, uint nCols) {
  std::ofstream File(FileName.c_str());
  File << format("! %s, %i x %i\n") % Desc % nRows % nCols;
  if (!File.good())
    throw std::runtime_error("failed to open file '" + FileName +
                             "' for writing.");
  for (uint iRow = 0; iRow < nRows; ++iRow) {
    for (uint iCol = 0; iCol < nCols; ++iCol) {
      File << format(" %21.14f") % pData[iRow + nRows * iCol];
    };
    File << "\n";
  }
  //     std::cout << format("*wrote 1-RDM to file %s.") % rdm1 << std::endl;
}

void SkipWhiteSpace(std::istream &str) {
  for (;;) {
    if (isspace(str.peek()))
      str.get();
    else
      break;
  };
};

void Skip(std::istream &str, char const *p) {
  SkipWhiteSpace(str);
  std::stringstream ss;
  for (uint iPos = 0; p[iPos] != 0; ++iPos) {
    int c = str.get();
    ss << (char)c;
    if (c != p[iPos]) {
      throw std::runtime_error("expected '" + std::string(p) +
                               "' in input, but found: '" + ss.str() + "'");
    }
  }
};

void FFciData::ClearAll() {
  nOrb = 0;
  nElec = 0;
  Ms2 = 0;
  iWfSym = 0;

  nPairs = 0;
  for (uint i = 0; i < 8; ++i)
    nSyOrb[i] = 0;

  IntClass = INTCLASS_Spatial;
  C1_Integrals = false;
  HubbardSys = false;
  iSyOrb.clear();

  CoreH.clear();
  CoreH_A.clear();
  CoreH_B.clear();
  Int2e.clear();
  Int2e_AA.clear();
  Int2e_BB.clear();
  Int2e_AB.clear();
  HubU.clear();
  fCoreEnergy = 0.;
};

static void CheckIndexBounds(uint i, uint j, uint k, uint l, uint nOrb) {
  if (i > nOrb || j > nOrb || k > nOrb || l > nOrb) {
    throw std::runtime_error(str(format("Encountered invalid integral index: "
                                        "%i %i %i %i. nOrb is only %i.") %
                                 i % j % k % l % nOrb));
  }
}

static void ReadInt2e(FInt2eData &Int2e, double &fCoreEnergy, bool &HubbardSys,
                      uint nOrb, std::ifstream &File, int iSymAB = 1) {
  uint nPairs = nOrb * (nOrb + 1) / 2;
  Int2e.Init(nPairs, FMatrix1::INIT_Clear);
  // ^- note: in symmetry case not all the data is overwritten, this is
  //          why we must clear & init to zero.
  while (File.good()) {
    double f;
    uint i, j, k, l;
    std::streampos ig = File.tellg();
    File >> f >> i >> j >> k >> l;
    CheckIndexBounds(i, j, k, l, nOrb);
    if (i == 0 && j == 0 && k == 0 && l == 0) {
      // frozen core energy or separator between spin blocks.. depends
      // on context.
      fCoreEnergy = f;
      return;
    } else if (k == 0 && l == 0) {
      // 1e- data... we're not supposed to be here. revert stream
      // position and return.
      File.seekg(ig);
      return;
    } else {
      // 2e- data
      if (k < l)
        std::swap(k, l);
      if (i < j)
        std::swap(i, j);
      if (i != j || i != k || i != l || j != k || j != l || k != l)
        HubbardSys = false;
      uint kl = ((k - 1) * k) / 2 + l - 1;
      uint ij = ((i - 1) * i) / 2 + j - 1;
      if (iSymAB)
        // ^- default case for spatial and AA/BB integrals, but
        // not for unrestricted AB integrals.
        Int2e[kl + nPairs * ij] = f;
      Int2e[ij + nPairs * kl] = f;
    };
  };
};

static void ReadInt1e(FInt1eData &Int1e, double &fCoreEnergy, uint nOrb,
                      std::ifstream &File) {
  Int1e.Init(nOrb, FMatrix1::INIT_Clear);
  while (File.good()) {
    double f;
    uint i, j, k, l;
    File >> f >> i >> j >> k >> l;
    CheckIndexBounds(i, j, k, l, nOrb);
    if (i == 0 && j == 0 && k == 0 && l == 0) {
      // frozen core energy or separator between spin blocks.. depends
      // on context.
      fCoreEnergy = f;
      return;
    } else if (k == 0 && l == 0) {
      if (i < j)
        std::swap(i, j);
      Int1e[(i - 1) + nOrb * (j - 1)] = f;
      Int1e[(j - 1) + nOrb * (i - 1)] = f;
    } else {
      // 2e- data
      throw std::runtime_error(
          "unexpected 2e integral data after 1e integral data.");
    };
  };
};

void FFciData::ReadDump(std::string const &FileName) {
  // format looks like this:
  // &FCI NORB=  9,NELEC=10,MS2= 0,
  //  ORBSYM=1,1,1,1,1,1,1,1,1,
  //  ISYM=1
  // &END

  std::ifstream File(FileName.c_str());
  if (!File.good())
    throw std::runtime_error("Failed to open input file '" + FileName + "'");
  Skip(File, "&FCI");

  ClearAll();
  uint const NotSet = 0xbadc0de;
  nOrb = NotSet;
  nElec = NotSet;

  std::string Field;
  while (1) {
    int iDummy;
    Field.clear();
    // read until next '=', until '&END', or end of input.
    SkipWhiteSpace(File);
    for (;;) {
      int c = File.get();
      if (!File.good())
        throw std::runtime_error("unexpected end of input while reading '" +
                                 FileName + "'");

      if ((char)c == '=' || Field == "&END" || Field == "/")
        break;
      Field.push_back((char)c);
    }
    // cout << "CMD = '" << Field << "'" << std::endl;
    if (Field == "&END" || Field == "/")
      break;
    else if (Field == "NORB")
      File >> nOrb;
    else if (Field == "NELEC")
      File >> nElec;
    else if (Field == "MS2")
      File >> Ms2;
    else if (Field == "ISYM") {
      File >> iWfSym;
      iWfSym -= 1;
    } else if (Field == "ORBSYM") {
      iSyOrb.resize(nOrb, 0);
      for (uint i = 0; i < 8; ++i)
        nSyOrb[i] = 0;
      for (uint i = 0; i < nOrb; ++i) {
        File >> iSyOrb[i];
        --iSyOrb[i];
        nSyOrb[iSyOrb[i]] += 1;
        Skip(File, ",");
      };
    } else if (Field == "IUHF") {
      File >> iDummy;
      if (iDummy == 1)
        IntClass = INTCLASS_Unrestricted;
    } else if (Field == "INTTYPE") {
      File >> iDummy;
      if (iDummy == 1)
        C1_Integrals = true;
    } else {
      File >> iDummy;
      cout << format("WARNING: Wave function/hamitlonian specification not "
                     "regcognized: %s=%i. Ignored!") %
                  Field % iDummy
           << std::endl;
    }

    // skip comma after field, if there should be one. Format is not
    // quite consistent in this regard.
    if ((char)File.peek() == ',')
      File.get();
  };

  SkipWhiteSpace(File);

  if (nOrb == NotSet)
    throw std::runtime_error("input is missing NORB field!");
  if (nElec == NotSet)
    throw std::runtime_error("input is missing NELEC field!");

  cout << format("*File '%s'  nOrb = %i%s  nElec = %i  Ms2 = %i  iWfSym = %i") %
              FileName % nOrb %
              (IntClass == INTCLASS_Spatial ? " (R)" : " (U)") % nElec % Ms2 %
              iWfSym
       << endl;

  // now read the 2e- & 1e- integral information. Symmetry packing
  // is not done here, since we don't actually have the multiplication
  // table.
  //
  // format looks like this:
  //    0.103195741742      2  1  2  1
  //     1.45640645877      2  2  1  1
  //   -0.287855353813E-01  2  2  2  1
  //     1.02254818418      2  2  2  2
  //   -0.342276029299E-17  3  1  1  1
  //    0.517081021881E-18  3  1  2  1
  //   -0.109500851660E-18  3  1  2  2
  //    0.480561321135E-01  3  1  3  1
  //    0.213874381351E-17  3  2  1  1
  //   -0.242735654929E-18  3  2  2  1
  //    0.377671163797E-18  3  2  2  2
  // i.e., i >= j; k >= l.
  //
  // In the unrestricted case, we get data in the following order:
  //   2e AA | 2e BB | 2e AB | 1e A | 1e B | CoreEnergy
  // here each | represents a line of "0.0 0 0 0 0" which is used
  // as separator in the format.

  nPairs = nOrb * (nOrb + 1) / 2;

  fCoreEnergy = 0;
  if (IntClass == INTCLASS_Spatial) {
    HubbardSys = true;
    ReadInt2e(Int2e, fCoreEnergy, HubbardSys, nOrb, File);
    ReadInt1e(CoreH, fCoreEnergy, nOrb, File);
  } else {
    assert(IntClass == INTCLASS_Unrestricted);
    HubbardSys = true;
    ReadInt2e(Int2e_AA, fCoreEnergy, HubbardSys, nOrb, File);
    ReadInt2e(Int2e_BB, fCoreEnergy, HubbardSys, nOrb, File);
    ReadInt2e(Int2e_AB, fCoreEnergy, HubbardSys, nOrb, File, 0);
    if (HubbardSys) {
      // check if hubbard-Us are equal for all the spin cases. If not,
      // disable HubU.
      for (uint i = 0; i < nOrb; ++i)
        if (Int2e_AA(i, i, i, i) != Int2e_AB(i, i, i, i) ||
            Int2e_AA(i, i, i, i) != Int2e_BB(i, i, i, i))
          HubbardSys = false;
    }

    ReadInt1e(CoreH_A, fCoreEnergy, nOrb, File);
    ReadInt1e(CoreH_B, fCoreEnergy, nOrb, File);
    if (File.good()) {
      // ^- depending on the previous cases we may still have or not
      //    have to read the core energy as extra integral.
      double f;
      uint i, j, k, l;
      File >> f >> i >> j >> k >> l;
      if (i != 0 || j != 0 || k != 0 || l != 0)
        throw std::runtime_error("unexpected final integral in UHF case.");
      fCoreEnergy = f;
    }

    //         PrintMatrixGen(cout, &Int2e_AA[0],nPairs,1,nPairs,nPairs,
    //         "Int2e_AA"); PrintMatrixGen(cout,
    //         &Int2e_BB[0],nPairs,1,nPairs,nPairs, "Int2e_BB");
    //         PrintMatrixGen(cout, &Int2e_AB[0],nPairs,1,nPairs,nPairs,
    //         "Int2e_AB");
    //
    //         PrintMatrixGen(cout, &CoreH_A[0],nOrb,1,nOrb,nOrb, "CoreH_A");
    //         PrintMatrixGen(cout, &CoreH_B[0],nOrb,1,nOrb,nOrb, "CoreH_B");
  }

  Finalize();
  cout << endl;
  CheckSanity();
};

void FFciData::Finalize() {
  if (HubbardSys) {
    HubU.resize(nOrb);
    FInt2eData *pInt2e = (IntClass == INTCLASS_Spatial) ? &Int2e : &Int2e_AA;
    for (uint i = 0; i < nOrb; ++i)
      HubU[i] = (*pInt2e)(i, i, i, i);
    //         std::cout << "*HubU = " << fmt::fl(HubU,6,3) << std::endl;
    //         C1_Integrals = true;
  } else
    HubU.clear();

  if (IntClass == INTCLASS_Spatial) {
    pCoreH_A = &CoreH[0];
    pCoreH_B = &CoreH[0];
    pInt2e_AA = &Int2e[0];
    pInt2e_BB = &Int2e[0];
    pInt2e_AB = &Int2e[0];
  } else {
    pCoreH_A = &CoreH_A[0];
    pCoreH_B = &CoreH_B[0];
    pInt2e_AA = &Int2e_AA[0];
    pInt2e_BB = &Int2e_BB[0];
    pInt2e_AB = &Int2e_AB[0];
  }
};

void FFciData::CheckSanity() {
  uint nMaxOrb = 8 * sizeof(FOrbPat), nElecA = (nElec + Ms2) / 2;
  if (nMaxOrb > 64)
    nMaxOrb = 64; // only have binomial_coefficient up to 64.

  if (nElecA > nOrb)
    throw std::runtime_error(
        str(format("Not enough orbitals (nOrb=%i) for %i alpha electrons!") %
            nOrb % nElecA));
  if (nOrb > nMaxOrb)
    throw std::runtime_error(
        str(format("For technical reasons this version of fci cannot use more "
                   "than %i orbitals.") %
            nMaxOrb));
  if ((Ms2 % 2) != (nElec % 2) || Ms2 > nElec)
    throw std::runtime_error(
        str(format("%i electrons cannot have spin S=%i/2.") % nElec % Ms2));
  if (iWfSym != 0)
    cout << "Warning: This version of fci does not support spatial symmetry. "
            "iWfSym is ignored."
         << std::endl;
};

#if 0
int main()
{
    FFciData
        FciData;
    FciData.ReadDump("/tmp/gk343/FCIDUMP");

};
#endif

// kate: indent-mode normal; indent-width 4;
