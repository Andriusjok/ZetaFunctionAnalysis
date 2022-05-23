#include <numbers>
#include <cmath>
#include <complex>
#include <boost/multiprecision/mpc.hpp>
#include <algorithm>
using namespace boost::multiprecision;
namespace zeta
{
    mpc_complex zetafast(mpc_complex s, float accuracy);
    mpc_complex zetaAMTCT(mpc_complex s, float d);
    mpc_complex zetaAMT(mpc_complex s, float d);
}
namespace utils
{
    mpc_complex compute_sum(std::function<mpc_complex(int n)> function, mpfr_float tolerance)
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
}