#include "zeta.hpp"
#include <functional>
#include <boost/math/tools/roots.hpp>
#include <iostream>

const static boost::math::normal distribution;

mpfr_float c(int k, mpfr_float k1, mpfr_float miu, mpfr_float sigma_n)
{
    if (k < k1)
    {
        return 1;
    }
    else
    {
        return 1 - cdf(distribution, (k - miu)/sigma_n);
    }
}

namespace zeta
{
    mpc_complex zetaNA(mpc_complex s, int d, int m, int j)
    {
        mpfr_float n, miu, sigma_n, z, k0, k1;
        int p, k;
        z = boost::math::quantile(distribution,1 - pow(10, -d));
        n = pi<mpf_float>() / 2 * s.imag() + (d + m) * log(10);
        if (j == 1)
        {
            n = (n + log(2) - log(log(2))) / (log(3 + sqrt(8)));
            miu = n / sqrt(2);
            sigma_n = sqrt(n) / pow(32, 0.25);
        }
        else
        {
            n = (n + log(2) - log(log(2))) / log(2);
            miu = n/2;
            sigma_n = sqrt(n)/2 ;
        }
        k0 = ceil(miu + z*sigma_n);
        k1 = miu - z*sigma_n;
        mpc_complex capital_s{0, 0};
        p = -1;
        for (int k = 0 ; k < k0; k++)
        {
            auto var = p*c(k, k1, miu, sigma_n) * exp(-s.real()*log(k+1)) * mpc_complex(cos(s.imag()*log(k+1)), -sin(s.imag()*log(k+1)));
            p = -p;
            capital_s += p*c(k, k1, miu, sigma_n) * exp(-s.real()*log(k+1)) * mpc_complex(cos(s.imag()*log(k+1)), -sin(s.imag()*log(k+1)));
        }
        return capital_s / (1 - 2 * exp(-s.real() * log(2)) * mpc_complex(cos(s.imag() * log(2)), -sin(s.imag() * log(2))));
    }
}