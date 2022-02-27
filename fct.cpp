#include <iostream>
#include <complex>
#include <bits/stdc++.h>
#include <cmath>
#include <math.h>

using std::complex;
using std::exp;
using std::real;
using std::cos;
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
        return new inum(x[0]);
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
        cosine_twiddles[k] = 2.0 * exp(-i * pi * double(k) / (2.0 * N));
    }
    return cosine_twiddles;
}


// double* fct(const double* x, int N, const inum* fourier_twiddles, const inum* cosine_twiddles) {

// }


int main() {
    int N = 16;
    inum* input_data = new inum[N];
    inum* output_data = new inum[N];

    // I got the tests from here:
    // https://www.stackoverflow.com/questions/25298119/unit-testing-a-discrete-fourier-transformation
    // specifically referencing:
    //      https://www.dsprelated.com/showthread/comp.dsp/71595-1.php

    // calculate twiddle factors for fourier transform
    inum* fourier_twiddles = calculate_fourier_twiddles(N);

    // TEST 1
    cout << "TEST 1:\n\tInput Data:\t";
    for (int k=0; k<N; k++) {
        input_data[k] = 0.0;
        cout << '|' << input_data[k];
    }
    cout << "|\n\tOutput Data:\t";
    output_data = fft(input_data, N, 1, fourier_twiddles);
    for (int k=0; k<N; k++) {
        cout << '|' << output_data[k];
    }
    cout << "|\n\n";

    // TEST 2
    cout << "TEST 2:\n\tInput Data:\t";
    for (int k=0; k<N; k++) {
        input_data[k] = 1.0;
        cout << '|' << input_data[k];
    }
    cout << "|\n\tOutput Data:\t";
    output_data = fft(input_data, N, 1, fourier_twiddles);
    for (int k=0; k<N; k++) {
        cout << '|' << output_data[k];
    }
    cout << "|\n\n";

    // TEST 3
    cout << "TEST 3:\n\tInput Data:\t";
    for (int k=0; k<N; k++) {
        input_data[k] = 1.0 - k % 2 * 2;
        cout << '|' << input_data[k];
    }
    cout << "|\n\tOutput Data:\t";
    output_data = fft(input_data, N, 1, fourier_twiddles);
    for (int k=0; k<N; k++) {
        cout << '|' << output_data[k];
    }
    cout << "|\n\n";

    // TEST 4
    cout << "TEST 4:\n\tInput Data:\t";
    for (int k=0; k<N; k++) {
        input_data[k] = exp(16 * pi * i * double(k) / double(N));
        cout << '|' << input_data[k];
    }
    cout << "|\n\tOutput Data:\t";
    output_data = fft(input_data, N, 1, fourier_twiddles);
    for (int k=0; k<N; k++) {
        cout << '|' << output_data[k];
    }
    cout << "|\n\n";

    // TEST 5
    cout << "TEST 5:\n\tInput Data:\t";
    for (int k=0; k<N; k++) {
        input_data[k] = cos(16 * pi * k / N);
        cout << '|' << input_data[k];
    }
    cout << "|\n\tOutput Data:\t";
    output_data = fft(input_data, N, 1, fourier_twiddles);
    for (int k=0; k<N; k++) {
        cout << '|' << output_data[k];
    }
    cout << "|\n\n";

    // TEST 6
    cout << "TEST 6:\n\tInput Data:\t";
    for (int k=0; k<N; k++) {
        input_data[k] = exp((double(43) / 7) * i * 2.0 * pi * double(k) / double(N));
        cout << '|' << input_data[k];
    }
    cout << "|\n\tOutput Data:\t";
    output_data = fft(input_data, N, 1, fourier_twiddles);
    for (int k=0; k<N; k++) {
        cout << '|' << output_data[k];
    }
    cout << "|\n\n";

    // TEST 7
    cout << "TEST 7:\n\tInput Data:\t";
    for (int k=0; k<N; k++) {
        input_data[k] = cos((double(43) / 7) * 2 * pi * k / N);
        cout << '|' << input_data[k];
    }
    cout << "|\n\tOutput Data:\t";
    output_data = fft(input_data, N, 1, fourier_twiddles);
    for (int k=0; k<N; k++) {
        cout << '|' << output_data[k];
    }
    cout << "|\n\n";

    // return successful exit code
    return 0;
}
