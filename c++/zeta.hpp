#include <numbers>
#include <cmath>
#include <complex>
#include <boost/multiprecision/mpc.hpp>
#include <boost/math/distributions/normal.hpp>
#include <algorithm>
using namespace boost::multiprecision;
using namespace boost::math::constants;
namespace zeta
{
    mpc_complex zetafast(mpc_complex s, mpfr_float accuracy);
    mpc_complex zetaAMTCT(mpc_complex s, int d);
    mpc_complex zetaNA(mpc_complex s, int d, int m, int j);
}
namespace utils
{
    inline mpc_complex compute_sum(std::function<mpc_complex(int n)> function, mpfr_float tolerance)
    {
      int n = 1;
      mpc_complex computed_sum;
      mpc_complex previous_sum;
      while ((abs(computed_sum-previous_sum).real() > tolerance and abs(computed_sum-previous_sum).imag() > tolerance) or n == 1)
      {
        previous_sum = computed_sum;
        computed_sum += function(n);
        n=n+1;
      }
      return computed_sum;
    }
    inline mpc_complex computeBinomialCoefficient(mpc_complex n, int k) {
        if (n == k || k == 0) {
            return 1;
        }
        mpc_complex out = n - k + 1;
        for (int i = 1; i < k; ++i) {
            out = out * (n - k + 1 + i) / (i + 1);
        }
        return out;
    }
    static const int g=7;
    static const double p[g+2] = {0.99999999999980993, 676.5203681218851,
    -1259.1392167224028, 771.32342877765313, -176.61502916214059,
    12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6,
    1.5056327351493116e-7};

    inline mpc_complex gamma(mpc_complex s)
    {
        if ( s.real()<0.5 ) {
            return pi<mpf_float>() / (sin(pi<mpf_float>()* s)*gamma(1.0-s));
        }
        s--;
        mpc_complex x = p[0];
        for (int i=1; i<g+2; i++) {
            x += p[i]/(s+mpc_complex(i,0));
        }
        mpc_complex t = s + (g + 0.5);
            return sqrt(2*pi<mpf_float>()) * pow(t,s + 0.5) * exp(-t) * x;
    }
}