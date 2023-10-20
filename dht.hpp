#include <complex>
#include <cmath>
/*
based on code from:
http://wwwa.pikara.ne.jp/okojisan/otfft-en/optimization1.html
*/

#include <complex>
#include <cmath>

//typedef std::complex<float> complex_t;

template <typename complex_t>
void fft0(int n, int s, bool eo, complex_t* x, complex_t* y)
// n  : sequence length
// s  : stride
// eo : x is output if eo == 0, y is output if eo == 1
// x  : input sequence(or output sequence if eo == 0)
// y  : work area(or output sequence if eo == 1)
{
    using T = decltype(x->real());
    const int m = n/2;

    const T theta0 = T(2)*M_PI/n;

    if (n == 2) {
        complex_t* z = eo ? y : x;
        for (int q = 0; q < s; q++) {
            const complex_t a = x[q + 0];
            const complex_t b = x[q + s];
            z[q + 0] = a + b;
            z[q + s] = a - b;
        }
    }
    else if (n >= 4) {
        for (int p = 0; p < m; p++) {
            const complex_t wp = complex_t(std::cos(p*theta0), -std::sin(p*theta0));
            for (int q = 0; q < s; q++) {
                const complex_t a = x[q + s*(p + 0)];
                const complex_t b = x[q + s*(p + m)];
                y[q + s*(2*p + 0)] =  a + b;
                y[q + s*(2*p + 1)] = (a - b) * wp ;
            }
        }
        fft0(n/2, 2*s, !eo, y, x);
    }
}

template <typename T>
void dht(T * begin, int n)
{
    std::vector<std::complex<T>> y(n);
    std::vector<std::complex<T>> x(begin, begin + n);

    fft0(n, 1, 0, x.data(), y.data());
    for (int i=0; i<n; ++i)
    {
        begin[i] = x[i].real() - x[i].imag();
    }
}

template <typename T>
void idht(T * begin, int n)
{
    dht(begin, n);
    T scale = T(1) / n;
    for (int i=0; i < n; ++i)
    {
        begin[i] *= scale;
    }
}
