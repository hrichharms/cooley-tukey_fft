/*
Radix-2 Out-of-Place DIT FFT Algorithm for 1D Real Input

TODO:
    - implement IDFT
    - replace complex number operations with macros
    - add SIMD build support
*/

#include <complex>
#include <math.h>
#include <iostream>

using std::cos;
using std::sin;


typedef double real;
typedef std::complex<real> complex;


const double pi=3.141592653589793238462643383279502884197169399375105820974944;


/*
Calculates twiddle factors (complex roots of unity) for an N/2-point DFT
and the associated post-processing for a real-valued 1D input
*/
complex* calculate_fft_twiddles(int N) {
    complex* twiddles = new complex[N/2];
    double phase;
    for (int k=0; k<N; k++) {
        phase = -2 * pi * k / N;
        twiddles[k] = complex(cos(phase), sin(phase));
    }
    return twiddles;
}


/*
Combine the outputs of two DFTs
*/
void r2_butterfly(
    complex* twiddles,
    complex* output,
    int stride,
    int m
) {
    complex* output2 = output + m;
    complex t;
    do {
        t = *output2 * *twiddles;
        *output2 = *output - t;
        *output += t;

        twiddles += stride * 2;
        output++;
        output2++;
    } while (--m);
}


/*
Radix-2 Cooley-Tukey FFT
*/
void fft_recursive(
    complex* twiddles,
    const complex* input,
    int n,
    complex* output,
    int stride
) {
    int m = n/2;

    if (m == 1) {
        output[0] = input[0];
        output[1] = input[stride];
    } else {
        fft_recursive(twiddles, input, m, output, 2*stride);
        fft_recursive(twiddles, input+stride, m, output+m, 2*stride);
    }

    r2_butterfly(twiddles, output, stride, m);
}


/*
Collapses real input of size N into complex sequence of size N/2,
calculates the N/2-point DFT, then extracts the N-point DFT of the
original N-point real input sequence
*/
void fft(
    complex* twiddles,
    const real* input,
    int N,
    complex* output
) {

    // collapse real input into N/2-point complex sequence
    complex* input_complex = new complex[N];
    for (int k=0; k<N/2; k++) {
        input_complex[k] = complex(input[2 * k], input[2 * k + 1]);
    }

    // perform N/2-point FFT on complex sequence
    fft_recursive(twiddles, input_complex, N/2, output, 1);

    // derive N-point DFT of input data from N/2-point DFT
    complex T, Tc, c3, c4, c5;
    for (int k=0; k<N/4; k++) {
        T = output[k];
        Tc = complex(output[N/2 - k].real(), -output[N/2 - k].imag());

        c3 = T + Tc;
        c4 = T - Tc;
        c5 = c4 * twiddles[k];

        output[k] = complex(0.5 * (c3.real() + c5.real()), 0.5 * (c3.imag() + c5.imag()));
        output[N/2 - k] = complex(0.5 * (c3.real() - c5.real()), 0.5 * (c5.imag() - c3.imag()));

    }
    for (int k=0; k<N; k++) {
        std::cout << '|' << output[k];
    }
    std::cout << "|\n";
}


int main() {

    int N = 4;
    real* input = new real[N];
    complex* output = new complex[N];

    complex* twiddles = calculate_fft_twiddles(N);

    input[0] = 0;
    input[1] = 1;
    input[2] = 1;
    input[3] = 0;

    fft(twiddles, input, N, output);

    return 0;
}
