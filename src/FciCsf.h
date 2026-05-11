/* This program is free software. It comes without any warranty, to the extent
 * permitted by applicable law. You may use it, redistribute it and/or modify
 * it, in whole or in part, provided that you do so at your own risk and do not
 * hold the developers or copyright holders liable for any claim, damages, or
 * other liabilities arising in connection with the software.
 *
 * Developed by Gerald Knizia, 2010--2012.
 */

#ifndef FCICSF_H
#define FCICSF_H

// make an explicit list of all bit patterns with nElec bits != 0 and all
// bits > nOrb zero.
void ListBitPatterns(FOrbPat *&pBeg, FOrbPat *&pEnd, uint nElec, uint nOrb,
                     FMemoryStack &Mem);

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
                      bool AbsorbSignsForStringOrder, FMemoryStack &Mem);

#endif // FCICSF_H
