/* This program is free software. It comes without any warranty, to the extent
 * permitted by applicable law. You may use it, redistribute it and/or modify
 * it, in whole or in part, provided that you do so at your own risk and do not
 * hold the developers or copyright holders liable for any claim, damages, or
 * other liabilities arising in connection with the software.
 *
 * Developed by Gerald Knizia, 2010--2012.
 */

#include <stdio.h>
#include <stdlib.h> // for exit()

void AicAssertFail(char const *pExpr, char const *pFile, int iLine) {
  printf("\n!! Assertion failed: '%s'\nat %s:%i\n", pExpr, pFile, iLine);
  exit(1);
};
