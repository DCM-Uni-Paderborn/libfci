/* This program is free software. It comes without any warranty, to the extent
 * permitted by applicable law. You may use it, redistribute it and/or modify
 * it, in whole or in part, provided that you do so at your own risk and do not
 * hold the developers or copyright holders liable for any claim, damages, or
 * other liabilities arising in connection with the software.
 *
 * Developed by Gerald Knizia, 2010--2012.
 */

#ifndef CXFORTRANINT_H
#define CXFORTRANINT_H

#include "FciCompat.h"
// Basic Fortran integer type used by CP2K's BLAS/LAPACK interface.
typedef boost::int32_t FORTINT;

typedef FORTINT const &FINTARG;
typedef double const &FDBLARG;

// BLAS/LAPACK symbols used by the CP2K toolchains are lowercase with one
// trailing underscore.
#define FORT_Extern(lowercasef, UPPERCASEF) lowercasef##_

#endif /* CXFORTRANINT_H */
