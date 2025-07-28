#ifndef THOMAS_SOLVER_H
#define THOMAS_SOLVER_H

#if __has_include(<hls_x_complex.h>)
#  include <hls_x_complex.h>
#  include <hls_stream.h>
#else
#  include "hls_stub.h"
#endif

#include <cmath>

const int N = 64;
using data_t = float;
using complex_t = hls::x_complex<data_t>;

void thomas_solver(complex_t a[N], complex_t b[N], complex_t c[N],
                   complex_t d[N], complex_t x[N]);

#endif // THOMAS_SOLVER_H
