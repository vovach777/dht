#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <climits>
#include <cfloat>
#include <cassert>
#include <complex>
#include <string>


template <typename complex_t>
void fft0(int n, int s, bool eo, complex_t* x, complex_t* y)
// n  : sequence length
// s  : stride
// eo : x is output if eo == 0, y is output if eo == 1
// x  : input sequence(or output sequence if eo == 0)
// y  : work area(or output sequence if eo == 1)
{
    const int m = n/2;
    const float theta0 = 2*M_PI/n;

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
            const complex_t wp = complex_t(cos(p*theta0), -sin(p*theta0));
            for (int q = 0; q < s; q++) {
                const complex_t a = x[q + s*(p + 0)];
                const complex_t b = x[q + s*(p + m)];
                y[q + s*(2*p + 0)] =  a + b;
                y[q + s*(2*p + 1)] = (a - b) * wp;
            }
        }
        fft0(n/2, 2*s, !eo, y, x);
    }
}



namespace Hartley {

template <typename T>    
inline void dht0(int n, T * real, float scale ) // Fourier transform
{
    assert( n>0 && (n & (n - 1))==0 );
    static_assert(std::is_floating_point<T>::value); 
    std::vector< std::complex<T> > x(real,real+n);
    std::vector< std::complex<T> > y(n);
    
    fft0(n, 1, 0, x.data(), y.data());    

    for (int i = 0; i < n; i++) {
        const auto v = x[i] * scale; // (x[i] + std::conj(x[n - i])) * scale;
        real[i] = (v.real() - v.imag());
    }
}
        

    template <typename T>    
    void dht(int n, T * real ) 
    {
        dht0(n, real, T(1) / n);
    }

    template <typename T>    
    void  idht(int n, T * real )
    {        
       dht0(n, real, T(1));
    }  

  
    template <typename T>    
    void dht(std::vector<T> &real) 
    {
        dht0(real.size(), real.data(), T(1) / real.size());
    }

    template <typename T>    
    void  idht(std::vector<T> &real)
    {        
        dht0(real.size(), real.data(), T(1));
    }  


}

    using Hartley::dht;
    using Hartley::idht;

template <typename container>
void print_signal(container && c)
{
    for (auto value : c)
    {
        std::cout   << std::fixed
                    << std::setprecision(4)
                    << std::setw(8)
                    << value
                    << std::endl;
    }
    std::cout << std::endl;
}

int main() {
    
 
    std::vector<float> signal{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
  
    print_signal(signal);

    dht(signal);
    print_signal(signal);

    idht(signal);
 
    print_signal(signal);


}

