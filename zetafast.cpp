#include "zeta.hpp"
#include <functional>
#include <boost/math/tools/roots.hpp>
#include <iostream>
#include <cmath>
using namespace boost::math::constants;
int resolveV(mpc_complex s, mpfr_float accuracy)
{
    if (s.real() >= 1)
    {
        return ceil(log(8 / accuracy)).convert_to<int>();
    }
    auto findV = [&s, &accuracy](double x) 
    {
        double evaluation = x 
        - std::max((1 - s.real().convert_to<double>()) / 2, 0.) 
        * std::log(1 / 2 + x + s.imag().convert_to<double>()) 
        - log(8 / accuracy.convert_to<double>());
        double derivative = 1
        - (std::max((1 - s.real().convert_to<double>()) / 2, 0.) 
        / (1 / 2 + x + s.imag().convert_to<double>()));
        auto result = std::make_tuple(evaluation, derivative);
        return result;
    };
    double guess = 10;
    double min = 0;
    double max = 50000;
    return ceil(boost::math::tools::newton_raphson_iterate(findV, guess, min, max, 10));
};

number<backends::mpfr_float_backend<0U, allocate_dynamic>, et_on> resolveCapitalN(
    int v, 
    number<backends::mpfr_float_backend<0U, allocate_dynamic>, et_on> tau)
{
    return sqrt(1.11 * (1 + ((0.5 + tau) / v)));
}

mpc_complex d_series(
    mpc_complex s, 
    int v, 
    number<backends::mpfr_float_backend<0U, allocate_dynamic>, et_on> N,
    mpfr_float accuracy)
{
    auto findGammaSum = [&s, &v, &N] (int n)
    {
        return pow(n, -s) * boost::math::gamma_q(v, n/N);
    };
    return utils::compute_sum(findGammaSum, accuracy);
}

mpc_complex e_miu_1_series(
    int m,
    mpc_complex s,
    number<backends::mpfr_float_backend<0U, allocate_dynamic>, et_on> N, 
    int sign, 
    int v)
{   
    mpc_complex sum = 0;
    int w = 0;
    while (w <= (v-1))
    {
        sum += utils::computeBinomialCoefficient(s-1, w) 
        * pow((m + (mpc_complex(0,sign)/(2 * pi<mpf_float>() * N))),(s-1-w)) 
        * ((mpc_complex(0,-sign)/pow((2 * pi<mpf_float>() * N),w)));
        w++;
    }
    return pow(m,(s - 1)) - sum;
}

mpc_complex e_1_series(
    mpc_complex s,
    number<backends::mpfr_float_backend<0U, allocate_dynamic>, et_on> N, 
    int sign, int v, 
    mpfr_float accuracy)
{
    auto find_e_1_series = [&s, &N, &sign, &v, &accuracy](double m)
    {
        return e_miu_1_series(m, s, N, sign, v);
    };
    return pow((2*pi<mpf_float>()), s-1) 
    * boost::math::tgamma(1 - s) 
    * exp(mpc_complex(0,sign * pi<mpf_float>()/2)*(1-s)) 
    * utils::compute_sum(find_e_1_series, accuracy);
}

mpc_complex gamma_function(
    mpc_complex s,
    int v,
    number<backends::mpfr_float_backend<0U, allocate_dynamic>, et_on> N)
{
    return (boost::math::tgamma(1-s + v) / ((1 - s) * boost::math::tgamma(v))) * pow(N, 1-s);
}

mpc_complex zetafast(mpc_complex s, mpfr_float accuracy)
{
    int v = resolveV(s, accuracy);
    number<backends::mpfr_float_backend<0U, allocate_dynamic>, et_on> N = resolveCapitalN(v, s.imag());
    std::cout<<v<<std::endl;
    std::cout<<N<<std::endl;
    if(0<= s.real() && s.real() <= 2 && s.imag() > 0 && accuracy <= 0.05)
    {
        double c = 0.5;
        //return d_series_lambda(s,v,N) + e_1_series_m_s(math.ceil(N),s, N, v) - gamma_function(s, v, N) + c*accuracy
    }
    if(s.imag() == 0)
    {
        double c = 0.5;
        //return d_series_lambda(s,v,N) - gamma_function(s,v,N) + c*accuracy
    }
    return d_series(s,v,N, accuracy) 
    + e_1_series(s, N, 1, v, accuracy) 
    + e_1_series(s, N, -1, v, accuracy) 
    - gamma_function(s, v, N);
}

int main() {
    mpc_complex s;
    mpfr_float epsilon = 10e-6;
    s.assign(50, 100);
    std::cout<<zetafast(s, epsilon)<<std::endl;
}