/* This program is free software. It comes without any warranty, to the extent
 * permitted by applicable law. You may use it, redistribute it and/or modify
 * it, in whole or in part, provided that you do so at your own risk and do not
 * hold the developers or copyright holders liable for any claim, damages, or
 * other liabilities arising in connection with the software.
 */

#ifndef LIBFCI_H
#define LIBFCI_H

#ifdef __cplusplus
extern "C" {
#endif

int libfci_solve(int nspins, int norb, int nelec, int ms2, const double *h1_a,
                 const double *h1_b, const double *eri_aa,
                 const double *eri_ab, const double *eri_bb, int max_iter,
                 double threshold, int pspace, int project_spin, int verbose,
                 double *energy, double *rdm1_a, double *rdm1_b,
                 char *err_msg, int err_msg_len);

int cp2k_fci_solve(int nspins, int norb, int nelec, int ms2,
                   const double *h1_a, const double *h1_b,
                   const double *eri_aa, const double *eri_ab,
                   const double *eri_bb, int max_iter, double threshold,
                   int pspace, int project_spin, int verbose, double *energy,
                   double *rdm1_a, double *rdm1_b, char *err_msg,
                   int err_msg_len);

#ifdef __cplusplus
}
#endif

#endif
