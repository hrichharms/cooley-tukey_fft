#include <iostream>
#include <complex>
#include <bits/stdc++.h>
#include <cmath>

using std::complex;
using std::exp;
using std::real;
using std::cout;


typedef complex<double> inum;

const double pi = M_PI;
const inum i(0.0, 1.0);


inum* calculate_fourier_twiddles(int N) {
    inum* fourier_twiddles;
    for (int k=0; k<N; k++) {
        fourier_twiddles[k] = exp(-2.0 * pi * i * double(k) / double(N));
    }
    return fourier_twiddles;
}


inum* fft(const inum* x, int N, int s, const inum* fourier_twiddles) {
    if (N == 1) {
        return new inum(v[0]);
    }
    inum* X = new inum[N];
    int new_N = N/2;
    inum* E = fft(x, new_N, 2*s, fourier_twiddles);
    inum* O = fft(&x[s], new_N, 2*s, fourier_twiddles);
    for (int k=0; k<new_N; k++) {
        X[k] = E[k] + fourier_twiddles[k] * O[k];
        X[k + new_N] = E[k] - fourier_twiddles[k] * O[k];
    }
    delete[] E;
    delete[] O;
    return X;
}


inum* calculate_cosine_twiddles(int N) {
    inum* cosine_twiddles;
    for (int k=0; k<N; k++) {
        cosine_twiddles[k] = 2 * exp(-i * pi * double(k) / (2.0 * N));
    }
    return cosine_twiddles;
}


// double* fct(const double* x, int N, const inum* fourier_twiddles, const inum* cosine_twiddles) {

// }


int main() {
    int N = 16;
    inum* data = new inum[N];

    // test for if all nums are 0.0
    for (int k=0; k<N; k++) {
        data[k] = 0.0;
    }

    return 0;
}