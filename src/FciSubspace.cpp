/* This program is free software. It comes without any warranty, to the extent
 * permitted by applicable law. You may use it, redistribute it and/or modify
 * it, in whole or in part, provided that you do so at your own risk and do not
 * hold the developers or copyright holders liable for any claim, damages, or
 * other liabilities arising in connection with the software.
 *
 * Developed by Gerald Knizia, 2010--2012.
 */

#include <set>

#include "Fci.h"
FHamiltonianData::FHamiltonianData(double *pInt2e_AA_, double *pInt2e_BB_,
                                   double *pInt2e_AB_, uint nPairs_,
                                   double *pInt1e_A_, double *pInt1e_B_,
                                   uint nOrb_, bool Absorb1e_, bool Spatial_)
    : pInt2e_AA(pInt2e_AA_), pInt2e_BB(pInt2e_BB_), pInt2e_AB(pInt2e_AB_),
      pInt1e_A(pInt1e_A_), pInt1e_B(pInt1e_B_), nPairs(nPairs_), nOrb(nOrb_),
      Absorb1e(Absorb1e_), Spatial(Spatial_) {
  if (Absorb1e) {
    pInt1e_A = 0; // not needed.
    pInt1e_B = 0;
  }
  assert(!Spatial || (pInt1e_A == pInt1e_B && pInt2e_AA == pInt2e_AB &&
                      pInt2e_AA == pInt2e_BB));
}

void FPSpace::Init(FConfList const &Confs_, FHamiltonianData const &Data_,
                   FOrbStringAdrTable const &AdrA_,
                   FOrbStringAdrTable const &AdrB_, FMemoryStack &Mem) {
  Data = Data_;
  Confs = Confs_;
  pAdrA = &AdrA_;
  pAdrB = &AdrB_;

  std::sort(Confs.begin(), Confs.end());

  // make a list of all alpha/beta strings which occur in at least
  // one of the supplied configurations.
  std::set<FAdr> StrsA, StrsB;
  for (FAdr iConf_ = 0; iConf_ < Confs.size(); ++iConf_) {
    FAdr iConf = Confs[iConf_], iStrA = iConf % pAdrA->nStr(),
         iStrB = iConf / pAdrA->nStr();
    StrsA.insert(iStrA);
    StrsB.insert(iStrB);
  }
  HeldAlphaStr.clear();
  HeldAlphaStr.reserve(StrsA.size());
  HeldBetaStr.clear();
  HeldBetaStr.reserve(StrsB.size());
  std::set<FAdr>::const_iterator itAdr;
  for (itAdr = StrsA.begin(); itAdr != StrsA.end(); ++itAdr)
    HeldAlphaStr.push_back(*itAdr);
  for (itAdr = StrsB.begin(); itAdr != StrsB.end(); ++itAdr)
    HeldBetaStr.push_back(*itAdr);

  // form the p-space Hamiltonian.
  MakeH(Mem);
};

FAdr FPSpace::FindConf(FAdr iConf) const {
  FAdrList::const_iterator itp =
      std::lower_bound(Confs.begin(), Confs.end(), iConf);
  if (itp == Confs.end() || *itp != iConf)
    return AdrNotFound;
  return itp - Confs.begin();
};

FAdr FPSpace::FindConf(FAdr iStrA, FAdr iStrB) const {
  return FindConf(iStrA + pAdrA->nStr() * iStrB);
};

bool FPSpace::HaveStr(FAdr iStr, FAdrList const &StrList) const {
  FAdrList::const_iterator itp =
      std::lower_bound(StrList.begin(), StrList.end(), iStr);
  return itp != StrList.end() && *itp == iStr;
};

void FPSpace::MakeH(FMemoryStack &Mem_) {
  // we do:  =
  //  <I|H|J> += [<I|c^r\alpha c_s\alpha|K> + <I|c^r\beta c_s\beta|K>] *
  //              (rs|tu) *
  //             [<K|c^t\alpha c_u\alpha|J> + <K|c^t\beta c_u\beta|J>]
  //
  FAdr nConf = nConfs();
  H.clear();
  H.resize(nConf * nConf, 0.0);
  uint nPairs = Data.nPairs, nOrb = Data.nOrb;

  FMemoryStackArray MemStacks(Mem_);
#pragma omp parallel
  {
    FMemoryStack &Mem = MemStacks.GetStackOfThread();
    //         FMemoryStack2 Mem(200000);
    std::vector<double> H_(nConf * nConf, 0.0);

// Behold the possibly slowest imaginable way of forming matrix
// elements between determinants!

// loop over K determinants
#pragma omp for schedule(dynamic)
    for (FAdr iStrKA = 0; iStrKA < pAdrA->nStr(); ++iStrKA) {
      for (FAdr iStrKB = 0; iStrKB < pAdrB->nStr(); ++iStrKB) {
        bool HaveStrKB = HaveStr(iStrKB, HeldBetaStr),
             HaveStrKA = HaveStr(iStrKA, HeldAlphaStr);
        if (!HaveStrKA && !HaveStrKB)
          continue;

        FOrbPat StrKA = pAdrA->MakePattern(iStrKA);
        FOrbPat StrKB = pAdrB->MakePattern(iStrKB);
        FSubstResult *pSubstA, *pSubstB;
        uint nSubstA, nSubstB;
        FormStringSubstsForSpin(pSubstA, nSubstA, 0, *pAdrA, StrKA, Mem);
        FormStringSubstsForSpin(pSubstB, nSubstB, 0, *pAdrB, StrKB, Mem);

        if (!Data.Absorb1e && HaveStrKA && HaveStrKB) {
          FAdr AdrK = FindConf(iStrKA, iStrKB);
          if (AdrK != AdrNotFound) {
            for (uint isaI = 0; isaI < nSubstA; ++isaI) {
              FSubstResult &saI = pSubstA[isaI];
              FAdr AdrPI = FindConf(saI.iStr, iStrKB);
              if (AdrPI == AdrNotFound)
                continue;
              H_[AdrPI + nConf * AdrK] +=
                  2. * saI.sign * Data.pInt1e_A[saI.k + nOrb * saI.l];
            }
            for (uint isbI = 0; isbI < nSubstB; ++isbI) {
              FSubstResult &sbI = pSubstB[isbI];
              FAdr AdrPI = FindConf(iStrKA, sbI.iStr);
              if (AdrPI == AdrNotFound)
                continue;
              H_[AdrPI + nConf * AdrK] +=
                  2. * sbI.sign * Data.pInt1e_B[sbI.k + nOrb * sbI.l];
            }
          }
          // ^- the 2. * counters the global .5 factor for 2e integrals below.
        }

        if (HaveStrKB)
          for (uint isaI = 0; isaI < nSubstA; ++isaI) {
            FSubstResult &saI = pSubstA[isaI];
            FAdr AdrPI = FindConf(saI.iStr, iStrKB);
            if (AdrPI == AdrNotFound)
              continue;
            for (uint isaJ = 0; isaJ < nSubstA; ++isaJ) {
              FSubstResult &saJ = pSubstA[isaJ];
              FAdr AdrPJ = FindConf(saJ.iStr, iStrKB);
              if (AdrPJ == AdrNotFound)
                continue;
              H_[AdrPI + nConf * AdrPJ] +=
                  saI.sign * saJ.sign *
                  Data.pInt2e_AA[i2e(saI.k, saI.l, saJ.k, saJ.l, nPairs)];
            }
          }
        if (HaveStrKA)
          for (uint isbI = 0; isbI < nSubstB; ++isbI) {
            FSubstResult &sbI = pSubstB[isbI];
            FAdr AdrPI = FindConf(iStrKA, sbI.iStr);
            if (AdrPI == AdrNotFound)
              continue;
            for (uint isbJ = 0; isbJ < nSubstB; ++isbJ) {
              FSubstResult &sbJ = pSubstB[isbJ];
              FAdr AdrPJ = FindConf(iStrKA, sbJ.iStr);
              if (AdrPJ == AdrNotFound)
                continue;
              H_[AdrPI + nConf * AdrPJ] +=
                  sbI.sign * sbJ.sign *
                  Data.pInt2e_BB[i2e(sbI.k, sbI.l, sbJ.k, sbJ.l, nPairs)];
            }
          }
        for (uint isaI = 0; isaI < nSubstA; ++isaI) {
          FSubstResult &saI = pSubstA[isaI];
          FAdr AdrPI = FindConf(saI.iStr, iStrKB);
          if (AdrPI == AdrNotFound)
            continue;
          for (uint isbJ = 0; isbJ < nSubstB; ++isbJ) {
            FSubstResult &sbJ = pSubstB[isbJ];
            FAdr AdrPJ = FindConf(iStrKA, sbJ.iStr);
            if (AdrPJ == AdrNotFound)
              continue;
            H_[AdrPI + nConf * AdrPJ] +=
                saI.sign * sbJ.sign *
                Data.pInt2e_AB[i2e(saI.k, saI.l, sbJ.k, sbJ.l, nPairs)];
          }
        }
        for (uint isbI = 0; isbI < nSubstB; ++isbI) {
          FSubstResult &sbI = pSubstB[isbI];
          FAdr AdrPI = FindConf(iStrKA, sbI.iStr);
          if (AdrPI == AdrNotFound)
            continue;
          for (uint isaJ = 0; isaJ < nSubstA; ++isaJ) {
            FSubstResult &saJ = pSubstA[isaJ];
            FAdr AdrPJ = FindConf(saJ.iStr, iStrKB);
            if (AdrPJ == AdrNotFound)
              continue;
            H_[AdrPI + nConf * AdrPJ] +=
                sbI.sign * saJ.sign *
                Data.pInt2e_AB[i2e(saJ.k, saJ.l, sbI.k, sbI.l, nPairs)];
            //                     H_[AdrPI + nConf*AdrPJ] += sbI.sign *
            //                     saJ.sign *
            //                     Data.GetInt2e(sbI.k,sbI.l,saJ.k,saJ.l);
          }
        }

        Mem.Free(pSubstA);
      };
    };

#pragma omp critical
    Add(&H[0], &H_[0], 1.0, nConf * nConf);
  }

  // 2e integrals are supposed to have a prefactor of .5, since we're summing
  // over all ijkl.
  for (FAdr i = 0; i < H.size(); ++i)
    H[i] *= .5;
  //     PrintMatrixGen(cout, &H[0],
  //     std::min(10ul,nConf),1,std::min(10ul,nConf),nConf, "p-space H: 10 x
  //     10");
};

void FPSpace::MakeHEvs() {
  Ev = H;
  Ew.resize(nConfs(), 0.0);
  Diagonalize(&Ew[0], &Ev[0], nConfs(), nConfs());
  //     cout << format("Lowest p-space root: %16.8f  [%i configurations]") %
  //     Ew[0] % nConfs() << std::endl;
};

FSubspaceStates::FSubspaceStates(uint nMaxDim_, uint DiisBlockSize_)
    : nDim(0), nMaxDim(nMaxDim_), iNext(0), nAmpLen(0), pRess(0), pAmps(0) {
  BlockSize = DiisBlockSize_;
  nDimUsed = 0;
  iNext = 0;
  iThis = 0;
};

FSubspaceStates::~FSubspaceStates() {
  ::free(pRess);
  ::free(pAmps);
};

void FSubspaceStates::GetBlockBoundaries(std::size_t &iOff, std::size_t &iBeg,
                                         std::size_t &nSize, uint iBlock) {
  std::size_t iEnd;
  assert(iBlock < nBlocks);
  iBeg = iBlock * BlockSize;
  iEnd = std::min(iBeg + BlockSize, nAmpLen);
  nSize = iEnd - iBeg;
  if (NoIO)
    iOff = iBeg * nMaxDim; // in-memory offset (into this->pAmps/pRess)
  else
    iOff = 0;
  //     std::cout << format("blk %i: %8i -- %8i [sz = %8i]") % iBlock % iBeg %
  //     iEnd % nSize << std::endl;
}

void FSubspaceStates::ReadBlock(std::size_t &iOff, std::size_t &iBeg,
                                std::size_t &nSize, uint iBlock) {
  GetBlockBoundaries(iOff, iBeg, nSize, iBlock);
  if (NoIO)
    return;
  //     std::cout << format("read %i: %8i -- %8i [sz = %8i] dim = %i") % iBlock
  //     % iBeg % (iBeg+nSize) % nSize % nDim << std::endl;
  std::size_t iBlockFileAdr = iBlock * (2 * nMaxDim * BlockSize);
  FileData.Read(&pAmps[iOff], nDim * nSize, iBlockFileAdr);
  FileData.Read(&pRess[iOff], nDim * nSize, iBlockFileAdr + nMaxDim * nSize);
}

void FSubspaceStates::WriteBlock(uint iBlock, uint iFirst, uint iLast) {
  if (NoIO)
    return;
  iLast = std::min(iLast, nDim);
  std::size_t iBeg, nSize, iOff;
  GetBlockBoundaries(iOff, iBeg, nSize, iBlock);
  std::size_t iBlockFileAdr = iBlock * (2 * nMaxDim * BlockSize);
  //    std::cout << format("write: iblk = %i  ivecs=[%i..%i)") % iBlock %
  //    iFirst % iLast << std::endl;
  FileData.Write(&pAmps[iOff] + iFirst * nSize, (iLast - iFirst) * nSize,
                 iBlockFileAdr + nSize * iFirst);
  FileData.Write(&pRess[iOff] + iFirst * nSize, (iLast - iFirst) * nSize,
                 iBlockFileAdr + nSize * (nMaxDim + iFirst));
};

void PermuteRows(double *pData, std::size_t RowSt, std::size_t ColSt,
                 uint nRows, uint nCols, uint *pPerm, int iDir,
                 FMemoryStack &Mem) {
  if (iDir == 0)
    return;
  assert(iDir == +1 || iDir == -1);
  double *pTmp;
  Mem.Alloc(pTmp, nRows);

  // apply permutation to rows.
  for (uint iCol = 0; iCol < nCols; ++iCol) {
    for (uint iRow = 0; iRow < nRows; ++iRow) {
      if (iDir == +1)
        pTmp[iRow] = pData[pPerm[iRow] * RowSt + iCol * ColSt];
      else
        pTmp[pPerm[iRow]] = pData[iRow * RowSt + iCol * ColSt];
    }
    for (uint iRow = 0; iRow < nRows; ++iRow)
      pData[iRow * RowSt + iCol * ColSt] = pTmp[iRow];
  };

  Mem.Free(pTmp);
};

void PermuteMatrix(double *p, std::size_t ColSt, std::size_t N, uint *pPerm,
                   int iDirRow, int iDirCol, FMemoryStack &Mem) {
  PermuteRows(p, 1, ColSt, N, N, pPerm, iDirRow, Mem);
  PermuteRows(p, ColSt, 1, N, N, pPerm, iDirCol, Mem);
};

// void BasisChange2(double *pOut, uint ldOut,
//                   double const *pIn, uint ldIn,
//                   double const *L, uint ldL, uint nRowsL, uint nColsL,
//                   double const *R, uint ldR, uint nRowsR, uint nColsR,
//                   char Trans, FMemoryStack &Mem)
// {
//     double
//         *pT;
//     Mem.Alloc(pT, nRowsL * nColsR);
//     assert(Trans == 'N' || Trans == 'T');
//     if ( Trans == 'N') {
//         // Make Out = L * C * R^T
//         assert(ldL >= nRowsL && ldR >= nRowsR && ldOut >= nRowsL && ldIn >=
//         nColsL); Mxm( pT, 1, nRowsL,     L , 1, nRowsL,    pIn, 1, ldIn,
//         nRowsL, nColsL, nColsR ); Mxm( pOut , 1, ldOut,   pT, 1, nRowsL, R ,
//         nRowsR, 1,   nRowsL, nColsR, nRowsR );
//     } else {
//         // Make Out = L^T * C * R
//         assert(ldL >= nRowsL && ldR >= nRowsR && ldOut >= nColsL && ldIn >=
//         nRowsL); Mxm(    pT, 1, nRowsL, pIn, 1, ldIn,     R , 1, nRowsR,
//         nRowsL, nRowsR, nColsR ); Mxm(  pOut,   1, ldOut,  L , nRowsL, 1, pT,
//         1, nRowsL,   nColsL, nRowsL, nColsR );
//     }
//     Mem.Free(pT);
// }
//
// #include "fmt.h"

void DiagonalizeGen2(double *pEw, double *pH, uint ldH, double *pS, uint ldS,
                     uint N, FMemoryStack &Mem, uint *pnDimUsed = 0,
                     double ThrS = 1e-10) {
  // construct an orthogonal basis.

  Diagonalize(pEw, pS, ldS, N);
  // throw away overlap eigenvalues below ThrS.
  uint nSkip = 0;
  while (nSkip < N && pEw[nSkip] < ThrS)
    ++nSkip;
  uint N1 = N - nSkip;
  // scale retained eigenvectors with 1/sqrt(Ew) to obtain an orthogonal basis.
  for (uint iEv = nSkip; iEv != N; ++iEv) {
    double f = 1.0 / std::sqrt(pEw[iEv]);
    for (uint iRow = 0; iRow < N; ++iRow)
      pS[iRow + iEv * ldS] *= f;
  }
  if (pnDimUsed != 0)
    *pnDimUsed = N - nSkip;

  // if ( nSkip != 0 ) std::cout << "NSKIP = " << nSkip << std::endl;
  double *pVecs = pS + ldS * nSkip;

  //     std::cout << fmt::fl(pEw,pEw+N,13,8," ") << std::endl;

  // diagonalize H in non-redundant orthogonal basis.
  double *pH1, *pT1;
  Mem.Alloc(pH1, N1 * N1);
  //     PrintMatrixGen(std::cout, pH, N,1, N,ldH, "H");
  //     BasisChange2(pH1,N1,  pH,ldH,  pS+ldS*nSkip,ldS,N,N1,
  //     pS+ldS*nSkip,ldS,N,N1, 'N', Mem);
  Mem.Alloc(pT1, N * N1);
  Mxm(pT1, 1, N, pH, 1, ldH, pVecs, 1, ldS, N, N, N1);
  Mxm(pH1, 1, N1, pVecs, ldS, 1, pT1, 1, N, N1, N, N1);
  Mem.Free(pT1);
  //     PrintMatrixGen(std::cout, pH1, N1,1, N1,N1, "H1/orth");
  Diagonalize(pEw, pH1, N1, N1);

  // transform eigenvectors back to original basis.
  Mxm(pH, 1, ldH, pVecs, 1, ldS, pH1, 1, N1, N, N1, N1);

  // set entries for discarded eigenvectors to zero.
  for (uint iEv = N - nSkip; iEv < N; ++iEv) {
    pEw[iEv] = 0.0;
    for (uint i = 0; i < N; ++i)
      pH[iEv * ldH + i] = 0.0;
    ;
  }
};

FStorageDevicePosixFs Fs;

void FSubspaceStates::Apply(double *pThisAmp, double *pThisRes, size_t nLength,
                            FPSpace &PSpace, FMemoryStack &Mem) {
  if (nDim < nMaxDim) {
    iNext = nDim;
    nDim += 1;
  }

  iThis = iNext;
  if (nLength != nAmpLen) {
    if (nAmpLen != 0)
      throw std::runtime_error(
          "FSubspaceStates: length of subspace vectors inconsistent.");
    nAmpLen = nLength;
    assert(pRess == 0 && pAmps == 0);
    BlockSize = std::min(nAmpLen, 1ul * BlockSize / sizeof(double)); // ~1mb.
    nBlocks = (nAmpLen + (BlockSize - 1)) / BlockSize;
    pRess = (double *)::malloc(nMaxDim * BlockSize * sizeof(double));
    pAmps = (double *)::malloc(nMaxDim * BlockSize * sizeof(double));
    FileData = Fs.AllocNewBlock(2 * nMaxDim * nAmpLen * sizeof(double));
    if (nBlocks != 1)
      std::cout << format("*%i subspace vector pairs stored in %i blocks of %i "
                          "kb [%i mb total].") %
                       nMaxDim % nBlocks % (sizeof(double) * BlockSize >> 10) %
                       ((sizeof(double) * 2 * nMaxDim * nAmpLen) >> 20)
                << std::endl;

    NoIO = nBlocks == 1;
    if (NoIO && omp_get_max_threads() != 1) {
      size_t nSplit = (4 * omp_get_max_threads());
      BlockSize =
          std::max((std::size_t)4096, (BlockSize + nSplit - 1) / nSplit);
      nBlocks = (nAmpLen + (BlockSize - 1)) / BlockSize;
    }

    //         nConfP = PSpace.nConfs();
    nConfP = 0;
    nMaxDiag = nConfP + nMaxDim;
    H.resize(nMaxDiag * nMaxDiag, 0.0);
    S.resize(nMaxDiag * nMaxDiag, 0.0);
    Ew.resize(nMaxDiag, 0.0);
    // copy p-space Hamiltonian.
    for (uint iP = 0; iP < nConfP; ++iP)
      for (uint jP = 0; jP < nConfP; ++jP)
        H[iP + nMaxDiag * jP] = PSpace.H[iP + nConfP * jP];
    for (uint iP = 0; iP < nConfP; ++iP)
      S[iP + nMaxDiag * iP] = 1.0;
  };
  nDiag = nConfP + nDim;

  //     // copy new vectors into their designated storing area.
  //     for ( uint iBlock = 0; iBlock < nBlocks; ++ iBlock ) {
  //         std::size_t iBeg, nSize, iOff;
  //         ReadBlock(iOff, iBeg, nSize, iBlock);
  //         memcpy(&pRess[iOff] + nSize*iThis, pThisRes + iBeg,
  //         sizeof(double)*nSize); memcpy(&pAmps[iOff] + nSize*iThis, pThisAmp
  //         + iBeg, sizeof(double)*nSize);
  //         // calculate new overlap and hamiltonian matrix elements (column
  //         iThis). Mxm(&H[nConfP*(nMaxDiag+1) + iThis*nMaxDiag],1,nMaxDiag,
  //         &pAmps[iOff],nSize,1, pThisRes+iBeg,1,nSize,  nDim,nSize,1,  iBlock
  //         != 0); Mxm(&S[nConfP*(nMaxDiag+1) + iThis*nMaxDiag],1,nMaxDiag,
  //         &pAmps[iOff],nSize,1, pThisAmp+iBeg,1,nSize,  nDim,nSize,1,  iBlock
  //         != 0); WriteBlock(iBlock, iThis, iThis+1);
  //     }
  memset(&H[nConfP * (nMaxDiag + 1) + iThis * nMaxDiag], 0,
         sizeof(H[0]) * nMaxDiag);
  memset(&S[nConfP * (nMaxDiag + 1) + iThis * nMaxDiag], 0,
         sizeof(H[0]) * nMaxDiag);
#pragma omp parallel if (NoIO)
  {
    bool First = true;
    TArray<double> H_(H.size(), 0.0), S_(S.size(), 0.0);

// copy new vectors into their designated storing area.
#pragma omp for
    for (uint iBlock = 0; iBlock < nBlocks; ++iBlock) {
      std::size_t iBeg, nSize, iOff;
      ReadBlock(iOff, iBeg, nSize, iBlock);
      memcpy(&pRess[iOff] + nSize * iThis, pThisRes + iBeg,
             sizeof(double) * nSize);
      memcpy(&pAmps[iOff] + nSize * iThis, pThisAmp + iBeg,
             sizeof(double) * nSize);
      // calculate new overlap and hamiltonian matrix elements (column iThis).
      Mxm(&H_[nConfP * (nMaxDiag + 1) + iThis * nMaxDiag], 1, nMaxDiag,
          &pAmps[iOff], nSize, 1, pThisRes + iBeg, 1, nSize, nDim, nSize, 1,
          !First);
      Mxm(&S_[nConfP * (nMaxDiag + 1) + iThis * nMaxDiag], 1, nMaxDiag,
          &pAmps[iOff], nSize, 1, pThisAmp + iBeg, 1, nSize, nDim, nSize, 1,
          !First);
      First = false;
      WriteBlock(iBlock, iThis, iThis + 1);
    }
#pragma omp critical
    Add(&H[0], &H_[0], 1.0, H.size());
#pragma omp critical
    Add(&S[0], &S_[0], 1.0, H.size());
  }

  // rebuild upper half of H and S.
  for (uint i = 0; i < nDiag; ++i)
    for (uint j = 0; j < i; ++j) {
      H[i + nMaxDiag * j] = H[j + nMaxDiag * i];
      S[i + nMaxDiag * j] = S[j + nMaxDiag * i];
    }

  //     PrintMatrixGen(std::cout, &H[0], nDiag, 1, nDiag, nMaxDiag, "subspace
  //     H"); PrintMatrixGen(std::cout, &S[0], nDiag, 1, nDiag, nMaxDiag,
  //     "subspace S");

  FScalarArray Ev = H, Cs = S;
  DiagonalizeGen(&Ew[0], &Ev[0], nMaxDiag, &Cs[0], nMaxDiag, nDiag);
  nDimUsed = nDiag;

  // store H and S in rotated EV basis.
  H.clear_data();
  S.clear_data();
  for (uint iEv = 0; iEv < nDim; ++iEv) {
    H[(nConfP + iEv) * (1 + nMaxDiag)] = Ew[iEv];
    S[(nConfP + iEv) * (1 + nMaxDiag)] = 1.;
  };

  if (1) {
// transform subspace vectors into (orthogonal) subspace
// hamiltonian eigenvalues. The main reason for this is that
// it allows us to throw away worse vectors when we have to evict
// old subspace vectors to save space.
// WARNING: numerically unstable...
#pragma omp parallel if (NoIO)
    {
      double *pTmp = (double *)::malloc(nMaxDim * BlockSize * sizeof(double));
      if (!pTmp)
        throw std::runtime_error("Davidson: Ran out of memory.");
#pragma omp for
      for (uint iBlock = 0; iBlock < nBlocks; ++iBlock) {
        std::size_t iBeg, nSize, iOff;
        ReadBlock(iOff, iBeg, nSize, iBlock);
        MxmLame(pTmp, 1, nSize, &pAmps[iOff], 1, nSize, &Ev[nConfP], 1,
                nMaxDiag, nSize, nDim, nDim);
        // ^- proper MxM here sometimes crashes OMP'd MKL in newer version.
        //    I don't see why; alignment looks fine to me.
        memcpy(&pAmps[iOff], pTmp, sizeof(double) * nSize * nDim);
        memcpy(pThisAmp + iBeg, pTmp, sizeof(double) * nSize);

        MxmLame(pTmp, 1, nSize, &pRess[iOff], 1, nSize, &Ev[nConfP], 1,
                nMaxDiag, nSize, nDim, nDim);
        memcpy(&pRess[iOff], pTmp, sizeof(double) * nSize * nDim);
        memcpy(pThisRes + iBeg, pTmp, sizeof(double) * nSize);
        WriteBlock(iBlock);
      };
      ::free(pTmp);
    }
    nDim = nDimUsed;
    iNext = std::min(nDim, nMaxDim - 1);
    // ^- worst vector is always at end in this scheme.
  }
}

void FSubspaceStates::OrthogonalizeUpdate(double *pThisAmp, size_t nLength,
                                          FPSpace &PSpace, FMemoryStack &Mem) {
  if (nDim < 2)
    return;
  uint nOrth = nDim;
  if (nDim == nMaxDim)
    nOrth -= 1;

  TArray<double> Sv_(nMaxDiag, 0.0), Hv_(nMaxDiag, 0.0);
  for (uint iBlock = 0; iBlock < nBlocks; ++iBlock) {
    std::size_t iBeg, nSize, iOff;
    ReadBlock(iOff, iBeg, nSize, iBlock);
    //         Mxm(&Hv_[0],1,nMaxDiag, &pAmps[iOff],nSize,1,
    //         pThisRes+iBeg,1,nSize,  nDim,nSize,1,  true);
    Mxm(&Sv_[0], 1, nMaxDiag, &pAmps[iOff], nSize, 1, pThisAmp + iBeg, 1, nSize,
        nOrth, nSize, 1, true);
  }

  for (uint i = 0; i < nOrth; ++i)
    Sv_[i] *= -1;
  PrintMatrixGen(std::cout, &Sv_[0], 1, 1, nOrth, 1, "S[Update] vs Prev");

  for (uint iBlock = 0; iBlock < nBlocks; ++iBlock) {
    std::size_t iBeg, nSize, iOff;
    ReadBlock(iOff, iBeg, nSize, iBlock);
    Mxm(pThisAmp + iBeg, 1, nSize, &pAmps[iOff], 1, nSize, &Sv_[0], 1, nMaxDiag,
        nSize, nOrth, 1, true);
  }
};

// kate: indent-width 4
