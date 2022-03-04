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

using std::cout;
using std::endl;
using std::conj;
using std::cos;
using std::sin;


typedef double real;
typedef std::complex<real> complex;


const complex i = complex(0.0, 1.0);
const double pi=3.141592653589793238462643383279502884197169399375105820974944;


/*
Calculates twiddle factors (complex roots of unity) for an N/2-point DFT
and the associated post-processing for a real-valued 1D input
*/
complex* calculate_fft_twiddles(int N, bool inverse=0) {
    complex* twiddles = new complex[N];
    double phase;
    for (int k=0; k<N; k++) {
        phase = -2.0 * pi * k / N;
        if (inverse) {phase *= -1;}
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

        twiddles += 2 * stride;
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
    complex* output,
    bool inverse = 0
) {

    // collapse real input into N/2-point complex sequence
    complex* input_complex = new complex[N];
    for (int k=0; k<N/2; k++) {
        input_complex[k] = complex(input[2 * k], input[2 * k + 1]);
    }

    // perform N/2-point FFT on complex sequence
    fft_recursive(twiddles, input_complex, N/2, output, 1);

    // derive N-point DFT of input data from N/2-point DFT
    output[N/2] = output[0].real() - output[0].imag();
    output[0] = output[0].real() + output[0].imag();
    output[3*N/4] = output[N/4];
    output[N/4] = conj(output[N/4]);
    complex T, Tc, c3, c4, c5;
    for (int k=1; k<N/4; k++) {
        T = output[k];
        Tc = conj(output[N/2 - k]);

        c3 = T + Tc;
        c4 = T - Tc;
        c5 = i * twiddles[k] * c4;

        output[k] = 0.5 * (c3 - c5);
        output[N/2 - k] = conj(0.5 * (c3 + c5));
        if (inverse) {
            output[k] *= 1.0 / N;
            output[N/2 - k] *= 1.0 / N;
        }
        output[N - k] = conj(output[k]);
        output[N/2 + k] = conj(output[N/2 - k]);

    }

}


/*
Calculates twiddle factors (complex roots of unity) for an N-point DCT
*/
complex* calculate_fct_twiddles(int N) {
    complex* twiddles = new complex[N];
    double phase;
    for (int k=0; k<N; k++) {
        phase = -pi * k / (2.0 * N);
        twiddles[k] = complex(cos(phase), sin(phase));
    }
    return twiddles;
}


/*

*/
void fct(
    const complex* twiddles,
    real* input,
    int N,
    real* output,
    complex* fft_twiddles,
    complex* fft_buffer
) {

    // re-order input sequence (VERY UGLY RIGHT NOW)
    for (int k=0; k<N; k++) { // copy input sequence to fft_buffer
        fft_buffer[k] = input[k];
    }
    for (int k=0; k<N; k++) { // calculate re-ordered sequence
        input[k % 2 * N / 2 + k / 2] = fft_buffer[k].real();
    }
    for (int k=0; k<N; k++) {
        cout << '|' << input[k];
    }
    cout << "|\n";

    // compute FFT of re-ordered input
    fft(fft_twiddles, input, N, fft_buffer, 1);
    for (int k=0; k<N; k++) {
        cout << '|' << fft_buffer[k];
    }
    cout << "|\n";

    // multiply DFT output sequence by DCT twiddle factors and extract
    // real-valued DCT output sequence
    for (int k=0; k<N/2; k++) {
        output[k] = (twiddles[k] * fft_buffer[k]).real();
    }

}


int main() {

    int N = 16;
    real* input = new real[N];
    real* output = new real[N];
    complex* fft_buffer = new complex[N];

    complex* fft_twiddles = calculate_fft_twiddles(N);
    complex* fct_twiddles = calculate_fct_twiddles(N);

    input[0] = 0;
    input[1] = 1;
    input[2] = 1;
    input[3] = 0;
    input[4] = 1;
    input[5] = 1;
    input[6] = 1;
    input[7] = 0;
    input[8] = 0;
    input[9] = 1;
    input[10] = 0;
    input[11] = 1;
    input[12] = 1;
    input[13] = 0;
    input[14] = 0;
    input[15] = 0;
    // for (int k=0; k<N; k++) {
    //     input[k] = k;
    // }

    for (int k=0; k<N; k++) {
        cout << '|' << input[k];
    }
    cout << "|\n";
    fct(fct_twiddles, input, N, output, fft_twiddles, fft_buffer);

    for (int k=0; k<N; k++) {
        cout << '|' << output[k];
    }
    cout << endl;

    return 0;
}
