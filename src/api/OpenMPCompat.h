#ifndef RABBITSKETCH_API_OPENMP_COMPAT_H
#define RABBITSKETCH_API_OPENMP_COMPAT_H

#ifdef _OPENMP
#include <omp.h>
#else
// Pragmas are ignored by non-OpenMP compilers; these small fallbacks keep the
// same code paths correct and single-threaded on portable toolchains.
inline int omp_get_thread_num() noexcept { return 0; }
inline int omp_get_max_threads() noexcept { return 1; }
inline int omp_get_num_threads() noexcept { return 1; }
#endif

#endif // RABBITSKETCH_API_OPENMP_COMPAT_H
