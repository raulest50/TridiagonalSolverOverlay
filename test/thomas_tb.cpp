#include "thomas_solver.h"
#include <complex>
#include <vector>
#include <fstream>
#include <iostream>
#include <cassert>

static void read_vec(const char* fn, std::vector<complex_t>& v) {
    std::ifstream in(fn);
    if (!in) {
        std::cerr << "Cannot open " << fn << "\n";
        std::exit(1);
    }
    data_t re, im;
    while (in >> re >> im)
        v.emplace_back(re, im);
}

static void write_vec(const char* fn, const std::vector<complex_t>& v) {
    std::ofstream out(fn);
    for (auto& z : v)
        out << z.real() << ' ' << z.imag() << '\n';
}

// CPU reference using std::complex<float>
static std::vector<std::complex<float>> ref_solver(
        const std::vector<std::complex<float>>& a,
        const std::vector<std::complex<float>>& b,
        const std::vector<std::complex<float>>& c,
        const std::vector<std::complex<float>>& d) {
    int n = d.size();
    std::vector<std::complex<float>> cprime(n);
    std::vector<std::complex<float>> dprime(n);
    std::vector<std::complex<float>> x(n);
    cprime[0] = c[0] / b[0];
    dprime[0] = d[0] / b[0];
    for (int i=1;i<n;i++) {
        std::complex<float> denom = b[i] - a[i]*cprime[i-1];
        cprime[i] = c[i]/denom;
        dprime[i] = (d[i]-a[i]*dprime[i-1])/denom;
    }
    x[n-1]=dprime[n-1];
    for(int i=n-2;i>=0;i--)
        x[i]=dprime[i]-cprime[i]*x[i+1];
    return x;
}

int main() {
    std::vector<complex_t> A, B, C, D;
    read_vec("test/a.dat", A);
    read_vec("test/b.dat", B);
    read_vec("test/c.dat", C);
    read_vec("test/d.dat", D);
    assert(A.size()==N && B.size()==N && C.size()==N && D.size()==N);

    // Run HLS function
    std::vector<complex_t> X(N);
    thomas_solver(A.data(), B.data(), C.data(), D.data(), X.data());

    // Compute reference on CPU for comparison
    std::vector<std::complex<float>> a_f(N), b_f(N), c_f(N), d_f(N);
    for(int i=0;i<N;i++) {
        a_f[i] = std::complex<float>(float(A[i].real()), float(A[i].imag()));
        b_f[i] = std::complex<float>(float(B[i].real()), float(B[i].imag()));
        c_f[i] = std::complex<float>(float(C[i].real()), float(C[i].imag()));
        d_f[i] = std::complex<float>(float(D[i].real()), float(D[i].imag()));
    }
    auto ref = ref_solver(a_f,b_f,c_f,d_f);

    double err=0.0;
    for(int i=0;i<N;i++) {
        double dre = float(X[i].real()) - ref[i].real();
        double dim = float(X[i].imag()) - ref[i].imag();
        err += dre*dre + dim*dim;
    }
    err = std::sqrt(err/N);
    std::cout << "RMS error: " << err << std::endl;

    std::vector<complex_t> ref_out(N);
    for(int i=0;i<N;i++)
        ref_out[i] = complex_t(data_t(ref[i].real()), data_t(ref[i].imag()));
    write_vec("test/golden.dat", ref_out);
    return 0;
}
