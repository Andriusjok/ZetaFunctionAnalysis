#include "zeta.hpp"
#include <functional>
#include <boost/math/tools/roots.hpp>
#include <iostream>
const static boost::math::normal distribution;

mpfr_float C(mpfr_float k0, mpfr_float k1, mpfr_float miu_n, mpfr_float sigma_n, int k)
{
    if (k0 < k1)
    {
        return 1;
    }
    return 1 - cdf(distribution, (k - miu_n)/sigma_n);
}

mpfr_float omega(mpfr_float sigma, mpfr_float t)
{
    return 1-(pow(2,2-sigma))*cos(t*log(2))+pow(4,1-sigma);
}
namespace zeta
{
    mpc_complex zetaAMTCT(mpc_complex s, int d)
    {
        mpfr_float sigma = s.real();
        mpfr_float t = s.imag();
        mpfr_float n = pi<mpfr_float>()* abs(t) 
            + (1 + 2*sigma)*log(abs(t)) 
            + 2*log(boost::math::tgamma(sigma)) 
            - 2 * log(abs(1 - pow(2,(1-sigma)))) 
            + 2 *d*log(10) 
            + log(2) 
            - log(pi<mpfr_float>());
        n = n/(2*log(3+sqrt(8)));
        mpfr_float epsilon = pow(10, -d);
        mpfr_float miu_n = n/sqrt(2);
        mpfr_float sigma_n = sqrt(n)/pow(32,0.25);
        mpfr_float z = boost::math::quantile(distribution,1-epsilon);
        mpfr_float k0 = ceil(miu_n + (z*sigma_n));
        mpfr_float k1 = floor(miu_n - (z*sigma_n));
        mpfr_float ReZ = 0;
        mpfr_float ImZ = 0;
        for (auto k = 0; k<= k0; k++)
        {
            mpfr_float a = pow(-1,k)*C(k0,k1,miu_n,sigma_n,k) / pow(k+1,sigma);
            mpfr_float b1 = t*log(k+1.);
            mpfr_float b2 = t*log((k+1)/2.);
            ReZ+= a*(cos(b1)-pow(2,1-sigma)*cos(b2));
            ImZ+= a*(sin(b1)-pow(2,1-sigma)*sin(b2));
        }
        return mpc_complex(ReZ, -ImZ)/omega(sigma, t);
    }
}