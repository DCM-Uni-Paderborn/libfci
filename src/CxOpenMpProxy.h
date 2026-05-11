/* This program is free software. It comes without any warranty, to the extent
 * permitted by applicable law. You may use it, redistribute it and/or modify
 * it, in whole or in part, provided that you do so at your own risk and do not
 * hold the developers or copyright holders liable for any claim, damages, or
 * other liabilities arising in connection with the software.
 *
 * Developed by Gerald Knizia, 2010--2012.
 */

// include OpenMP header if available or define inline dummy functions
// for OMP primitives if not.
#ifndef CXOPENMPPROXY_H
#define CXOPENMPPROXY_H

#ifdef _OPENMP
#include <omp.h>
#else
inline int omp_get_thread_num() { return 0; } // current thread id
inline void omp_set_num_threads(int) {};
inline int omp_get_max_threads() {
  return 1;
} // total number of threads supposed to be running.

struct omp_lock_t {};
inline void omp_destroy_lock(omp_lock_t *) {};
inline void omp_init_lock(omp_lock_t *) {};
inline void omp_set_lock(omp_lock_t *) {};
#endif

#endif // CXOPENMPPROXY_H
