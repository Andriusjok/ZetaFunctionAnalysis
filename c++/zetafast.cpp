#include "zeta.hpp"
#include <functional>
#include <boost/math/tools/roots.hpp>
#include <iostream>

static const double lambda_var = 3.151;

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

mpfr_float q_gamma(int v, mpfr_float m)
{
    mpfr_float sum = 0;
    for (unsigned w = 0; w <= v-1; w++)
    {
        sum += pow(m,w) / boost::math::factorial<mpf_float>(w) * pow(e<mpf_float>(), -m);
    }
    return sum;
}

mpfr_float resolveCapitalN(
    int v, 
    mpfr_float tau)
{
    return sqrt(1.11 * (1 + ((0.5 + tau) / v)));
}

mpc_complex d_series(
    mpc_complex s, 
    int v, 
    mpfr_float N,
    mpfr_float accuracy)
{
    auto findGammaSum = [&s, &v, &N] (int n)
    {
        return pow(n, -s) * q_gamma(v, n/N);
    };
    return utils::compute_sum(findGammaSum, accuracy);
}

mpc_complex e_miu_1_series(
    int m,
    mpc_complex s,
    mpfr_float N, 
    int sign, 
    int v)
{   
    mpc_complex sum = 0;
    for (int w = 0; w <= (v-1); w++)
    {
        sum += utils::computeBinomialCoefficient(s-1, w) 
        * pow(m + (mpc_complex(0,sign)/(2 * pi<mpf_float>() * N)), (s-1-w)) 
        * pow(mpc_complex(0,-sign)/(2*pi<mpf_float>()*N),w);
    }
    return pow(m,(s - 1)) - sum;
}

mpc_complex e_1_series(
    mpc_complex s,
    mpfr_float N, 
    int sign, int v, 
    mpfr_float accuracy)
{
    auto find_e_1_series = [&s, &N, &sign, &v, &accuracy](double m)
    {
        return e_miu_1_series(m, s, N, sign, v);
    };
    return pow(2*pi<mpf_float>(),s-1)
        * utils::gamma(1 - s)
        * pow(e<mpf_float>(),mpc_complex(0,sign * pi<mpf_float>()/2)*(1-s)) 
        * utils::compute_sum(find_e_1_series, accuracy);
}

mpc_complex d_series_lambda(
    mpc_complex s, 
    int v,
    mpfr_float N)
{
    auto f = [&s, &N, &v](int n)
    {
        return pow(n, -s) * q_gamma(v, n/N);
    };
    mpc_complex sum {0, 0};
    for (auto n = 1; n <= ceil(lambda_var*v*N).convert_to<int>(); n++)
    {
        sum+=f(n);
    }
    return sum;
}

mpc_complex e_1_series_m_s(int capital_M, mpc_complex s, mpfr_float N, int v)
{
    mpc_complex sum {0,0};
    for (auto m = 1; m <= capital_M; m++)
    {
        sum+=e_miu_1_series(m, s, N, 1, v);
    }
    return pow(2*pi<mpf_float>(),(s-1))
    * utils::gamma(1-s)
    * pow(e<mpf_float>(),mpc_complex(0,pi<mpf_float>()/2)*(1-s))
    * sum;
}
mpc_complex gamma_function(
    mpc_complex s,
    int v,
    mpfr_float N)
{
    return (utils::gamma(1-s + v) / ((1 - s) * utils::gamma(v))) * pow(N, 1-s);
}
namespace zeta
{
    mpc_complex zetafast(mpc_complex s, mpfr_float accuracy)
    {
        int v = resolveV(s, accuracy);
        mpfr_float N = resolveCapitalN(v, s.imag());
        if(0<= s.real() && s.real() <= 2 && s.imag() > 0 && accuracy <= 0.05)
        {
            return d_series_lambda(s,v,N) + e_1_series_m_s(ceil(N).convert_to<int>(),s, N, v) - gamma_function(s, v, N);
        }
        if(s.imag() == 0)
        {
            return d_series_lambda(s,v,N) - gamma_function(s,v,N);
        }
        return d_series(s,v,N, accuracy) 
        + e_1_series(s, N, 1, v, accuracy) 
        + e_1_series(s, N, -1, v, accuracy) 
        - gamma_function(s, v, N);
    }
}