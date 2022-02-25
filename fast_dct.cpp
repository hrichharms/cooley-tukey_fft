#include <iostream>
#include <complex>
#include <bits/stdc++.h>
#include <cmath>

using std::complex;
using std::exp;
using std::real;
using std::cout;


const complex<double> i(0.0, 1.0);


complex<double>* fft(const complex<double>* y, int N, int s) {
    if (N == 1) {
        return new complex<double>(y[0]);
    }
    complex<double>* X = new complex<double>[N];
    complex<double>* ft_1 = fft(y, N/2, 2*s);
    complex<double>* ft_2 = fft(y+s, N/2, 2*s);
    complex<double> p, q;
    ft_1 = fft(y, N / 2, 2 * s);
    ft_2 = fft(y + s, N / 2, 2 * s);
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


complex<double>* even_extension(const double* x, int N) {
    complex<double>* y = new complex<double>[2 * N];
    for (int k=0; k<N; k++) {
        y[k] = x[k];
        y[2 * N - k - 1] = x[k];
    }
    return y;
}


double* fct(const double* x, int N) {
    complex<double>* y = even_extension(x, N);
    complex<double>* V = fft(y, N, 1);
    double* C = new double[N];
    for (int k=0; k<N; k++) {
        C[k] = real(V[k]);
    }
    return C;
}


int main() {

    int N = 8;
    double* x = new double[N];
    for (int k=0; k<N; k++) {
        x[k] = k;
    }

    double* C = fct(x, N);
    for (int k=0; k<N; k++) {
        cout << '|' << C[k];
    }
    cout << '|\n';

    return 0;
}
