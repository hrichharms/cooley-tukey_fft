#include <iostream>
#include <complex>
#include <bits/stdc++.h>
#include <cmath>

using std::complex;
using std::exp;
using std::real;
using std::cout;


const complex<double> i(0.0, 1.0);


complex<double>* fft(const complex<double>* v, int N, int s) {
    if (N == 1) {
        return new complex<double>(v[0]);
    }
    complex<double>* X = new complex<double>[N];
    complex<double>* ft_1 = fft(v, N/2, 2*s);
    complex<double>* ft_2 = fft(v+s, N/2, 2*s);
    complex<double> p, q;
    ft_1 = fft(v, N / 2, 2 * s);
    ft_2 = fft(v + s, N / 2, 2 * s);
    for (int k=0; k<N/2; k++) {
        p = ft_1[k];
        q = exp(-2 * M_PI * i / (double) N * (double) k) * ft_2[k];
        X[k] = p + q;
        X[k + N/2] = p - q;
    }
    delete[] ft_1;
    delete[] ft_2;
    return X;
}


complex<double>* reorder_sequence(const double* x, int N) {
    complex<double>* v = new complex<double>[N];
    for (int k=0; k<N; k++) {
        if (k < N/2) {
            v[k] = x[2 * k];
        } else {
            v[k] = x[2 * N - 2 * k - 1];
        }
    }
    return v;
}


double* fct(const double* x, int N) {
    complex<double>* V = fft(reorder_sequence(x, N), N, 1);
    double* C = new double[N];
    for (int k=0; k<N/2; k++) {
        V[k] *= 2.0 * exp(-i * M_PI * double(k) / (2.0 * N));
        C[k] = V[k].real();
        C[N - k] = -V[k].imag();
    }
    return C;
}
