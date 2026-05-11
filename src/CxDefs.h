/* This program is free software. It comes without any warranty, to the extent
 * permitted by applicable law. You may use it, redistribute it and/or modify
 * it, in whole or in part, provided that you do so at your own risk and do not
 * hold the developers or copyright holders liable for any claim, damages, or
 * other liabilities arising in connection with the software.
 *
 * Developed by Gerald Knizia, 2010--2012.
 */

#ifndef CXDEFS_H
#define CXDEFS_H

#define AIC_NO_THROW throw()
#define CX_NO_THROW AIC_NO_THROW

// define macros for restricted pointer extensions. By defining a
// pointer as ``restricted'', we promise the compiler that
// the pointer is not aliased in the current scope

#define AIC_RP __restrict__

#define RESTRICT AIC_RP

#undef assert
#undef assert_rt
void AicAssertFail(char const *pExpr, char const *pFile, int iLine);
#define assert_rt(x)                                                           \
  if (x) {                                                                     \
  } else                                                                       \
    AicAssertFail(#x, __FILE__, __LINE__)

#define assert(x) assert_rt(x)

#endif // CXDEFS_H
