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
void fct_recursive(
    const complex* twiddles,
    real* input,
    real* output,
    complex* buffer
) {

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

}
