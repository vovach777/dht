#include <cmath>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include  <complex>
#include <string>
#include <sstream>
#include <iostream>
#include <vector>


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


std::string show(std::complex<float> v) 
{
    float max = std::max( std::abs( v.real() ),  std::abs( v.imag() ) );
    std::stringstream ss;
    ss << std::fixed << std::setprecision(1) << v;
    return max < 0.001 ? "" :  ss.str();
        
}


inline float sinus(int freq, int sampleRate, int i) {
    return std::sin( (2.0 * M_PI * (i % sampleRate) )*float(freq)/sampleRate );
}

template <typename F>
inline float gen(F && f, int freq, int sampleRate, int i) {    
    return f( (2.0 * M_PI * (i % sampleRate) )*float(freq)/sampleRate );
}


inline float genarg(int freq, int sampleRate, int i) {    
    return ( (2.0 * M_PI * (i % sampleRate) )*float(freq)/sampleRate ); 
}




std::vector<std::complex<float>> create_fir_for_overlap_save(int size_n)
{
    int fir_size = size_n;
    std::vector<std::complex<float>> fir(size_n);
    // for (int i = 0; i < fir_size; ++i)
    // {
    //     float x = M_PI*4 * (float(i-(size_n/2) + 1)/size_n); 
    //     fir[i] =  (std::abs(x) > 0.0001 ? sin(x) / x : 1.0f);
    // }

    for (int i = 0; i < fir_size; i++) {
        fir[i] = { float(i)/fir_size, float(-i)/fir_size + 0.5f };
    }    
    
    return fir;

}

//0,0,0,0.00003,0,-0.00008,0,0.00016,0,-0.00028,0,0.00043,0,-0.00061,0,0.00084,0,-0.00111,0,0.00143,0,-0.00180,0,0.00222,0,-0.00271,0,0.00326,0,-0.00388,0,0.00459,0,-0.00539,0,0.00629,0,-0.00733,0,0.00851,0,-0.00988,0,0.01148,0,-0.01337,0,0.01565,0,-0.01847,0,0.02208,0,-0.02688,0,0.03367,0,-0.04414,0,0.06271,0,-0.10553,0,0.31812,0.50000,0.31812,0,-0.10553,0,0.06271,0,-0.04414,0,0.03367,0,-0.02688,0,0.02208,0,-0.01847,0,0.01565,0,-0.01337,0,0.01148,0,-0.00988,0,0.00851,0,-0.00733,0,0.00629,0,-0.00539,0,0.00459,0,-0.00388,0,0.00326,0,-0.00271,0,0.00222,0,-0.00180,0,0.00143,0,-0.00111,0,0.00084,0,-0.00061,0,0.00043,0,-0.00028,0,0.00016,0,-0.00008,0,0.00003,0,0,0,


template <typename iterator>
class Slice {
    private:
    iterator begin_;
    iterator end_;
    public:
    using value_type = decltype(*(iterator{}));    
    using difference_type = decltype(std::distance(iterator{},iterator{}));
    explicit Slice(iterator begin_, iterator end_) : begin_(begin_), end_(end_) {};
    template <typename container>
    explicit Slice(container && c) : Slice(c.begin(), c.end()) {};
    constexpr auto begin() const noexcept { return begin_; }
    constexpr auto end() const noexcept { return end_; }
    constexpr auto slice( difference_type from=0, difference_type to=0) const noexcept {
        auto it_from = (from < 0) ? end_ + from : begin_ + from;
        auto it_to = (to <= 0) ?  end_ + to : begin_ + to;
        if (std::distance(it_from, it_to) < 0 )
            std::swap(it_from, it_to);
        return Slice(it_from,it_to);
    }
  constexpr size_t size() const noexcept {
        return std::distance(begin_, end_);
    }

    constexpr bool empty() const noexcept {
        return begin_ == end_;
    }    
};

template<typename T>
Slice(T) -> Slice<decltype(T{}.begin())>;


// template <typename Iterator>
// void processRange(Iterator begin, Iterator end) {

//     //using T = typename std::iterator_traits<Iterator>::value_type; //работает
//     using T = typename std::remove_reference<decltype(*begin)>::type;
//     if constexpr (std::is_same<T,float>::value )
//     {
//         std::cout << "float\n";
//     } else
//     if constexpr (std::is_same<T,double>::value )
//     {
//         std::cout << "double\n";        
//     } else
//     if constexpr (std::is_same<T,std::complex<float> >::value ) 
//     {
//         std::cout << "complex<float>\n";
//     } else
//     if constexpr (std::is_same<T,std::complex<double> >::value ) 
//     {
//         std::cout << "complex<double>\n";
//     } else
//     {
//         std::cout << "unknown\n";
//     }

// }

// template <typename Iterator, 
//         typename = typename std::enable_if<std::is_same<typename std::iterator_traits<Iterator>::value_type,
//             double
//         >::value>::type>
// void processRange(Iterator begin, Iterator end) {
//         std::cout << "double" << std::endl;
// }


namespace Hartley {

template <bool scale_f,typename T>    
auto __dht(int n, const T * real ) // Fourier transform
{
    //assert( n>0 && (n & (n - 1))==0 );
    static_assert(std::is_floating_point<T>::value); 

    std::vector< std::complex<T> > x(real,real+n);
    std::vector<std::complex<T>> y(n);
    std::vector<T> output(n);
    fft0(n, 1, 0, x.data(), y.data());    

   // Зеркальное отражение
    const auto scale = T(1) / n;
    for (int i = 0; i < n; i++) {
        const auto v = x[i] * (scale_f ? scale : T(1)); // (x[i] + std::conj(x[n - i])) * scale;
        output[i] = (v.real() - v.imag());
    }
    return output;
}
        

    template <typename T>    
    auto dht(int n, const T * real ) 
    {
        return __dht<true>(n, real);
    }


    template <typename T>    
    auto idht(const std::vector<T> & y)
    {    
        
        return __dht<false>(y.size(), y.data());
    }    
}



int main() {
    constexpr int size_n = 8;
    std::vector<float> v(size_n);
    for (int i = 0; i < size_n; ++i) {
        v[i] = i;// sinus(3,size_n,i);// (rand() ) / double(RAND_MAX) + 1;  //sinus(3,size_n,i)*0.5 + 0.123456; //gen(cosf,1,size_n,i) + gen(sinf,2,size_n,i) + gen(cosf,3,size_n,i) + gen(sinf,4,size_n,i);
        
    }

   for (auto val : v) {
        
        std::cout << std::fixed << std::setprecision(6) 
                   << std::setw(40) <<  val
                   
                   << std::endl;
    }

    std::cout << "--------\n";

    auto t = Hartley::dht(v.size(), v.data());
    
    for (auto val : t) {
        
        std::cout << std::fixed << std::setprecision(6) 
                   << std::setw(40) <<  val
                   
                   << std::endl;
    }

    std::cout << "--------\n";

    auto t_inv = Hartley::idht(t);
    for (auto val : t_inv) {
        
        std::cout << std::fixed << std::setprecision(6) 
                   << std::setw(40) <<  val
                   
                   << std::endl;
    }

    
}
