#ifndef OMP_TOOLS_H
#define OMP_TOOLS_H

#ifdef _OPENMP
#include <omp.h> // for tec constraints
#else
int omp_get_max_threads() { return 1; }
int omp_get_thread_num() { return 0; }
#endif // _OPENMP

#endif
