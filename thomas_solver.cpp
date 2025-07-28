#include "thomas_solver.h"

void thomas_solver(complex_t a[N], complex_t b[N], complex_t c[N],
                   complex_t d[N], complex_t x[N]) {
#pragma HLS INTERFACE s_axilite port=return bundle=CTRL
#pragma HLS INTERFACE m_axi     port=a offset=slave bundle=MA
#pragma HLS INTERFACE m_axi     port=b offset=slave bundle=MB
#pragma HLS INTERFACE m_axi     port=c offset=slave bundle=MC
#pragma HLS INTERFACE m_axi     port=d offset=slave bundle=MD
#pragma HLS INTERFACE m_axi     port=x offset=slave bundle=MX

    complex_t c_prime[N];
    complex_t d_prime[N];
#pragma HLS ARRAY_PARTITION variable=c_prime complete
#pragma HLS ARRAY_PARTITION variable=d_prime complete

    // Forward sweep
    c_prime[0] = c[0] / b[0];
    d_prime[0] = d[0] / b[0];
    for (int i = 1; i < N; ++i) {
#pragma HLS PIPELINE II=1
        complex_t denom = b[i] - a[i] * c_prime[i-1];
        c_prime[i] = c[i] / denom;
        d_prime[i] = (d[i] - a[i] * d_prime[i-1]) / denom;
    }

    // Back substitution
    x[N-1] = d_prime[N-1];
    for (int i = N-2; i >= 0; --i) {
#pragma HLS PIPELINE II=1
        x[i] = d_prime[i] - c_prime[i] * x[i+1];
    }
}
