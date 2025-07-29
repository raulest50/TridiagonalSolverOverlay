#include "thomas_solver.h"

void thomas_solver(complex_t dp,
                   complex_t dp1,
                   complex_t dp2,
                   complex_t off,
                   complex_t b[N],
                   complex_t x[N]) {
#pragma HLS INTERFACE s_axilite port=return bundle=CTRL
#pragma HLS INTERFACE s_axilite port=dp   bundle=CTRL
#pragma HLS INTERFACE s_axilite port=dp1  bundle=CTRL
#pragma HLS INTERFACE s_axilite port=dp2  bundle=CTRL
#pragma HLS INTERFACE s_axilite port=off  bundle=CTRL
#pragma HLS INTERFACE m_axi     port=b    offset=slave bundle=MB
#pragma HLS INTERFACE m_axi     port=x    offset=slave bundle=MX

    complex_t c_prime[N];
    complex_t d_prime[N];
#pragma HLS ARRAY_PARTITION variable=c_prime complete
#pragma HLS ARRAY_PARTITION variable=d_prime complete

    // Forward sweep
    c_prime[0] = off / dp1;
    d_prime[0] = b[0] / dp1;
    for (int i = 1; i < N - 1; ++i) {
#pragma HLS PIPELINE II=1
        complex_t denom = dp - off * c_prime[i-1];
        c_prime[i] = off / denom;
        d_prime[i] = (b[i] - off * d_prime[i-1]) / denom;
    }
    // Last element uses dp2
    d_prime[N-1] =
#pragma HLS PIPELINE II=1
        (b[N-1] - off * d_prime[N-2]) / (dp2 - off * c_prime[N-2]);

    // Back substitution
    x[N-1] = d_prime[N-1];
    for (int i = N-2; i >= 0; --i) {
#pragma HLS PIPELINE II=1
        x[i] = d_prime[i] - c_prime[i] * x[i+1];
    }
}
