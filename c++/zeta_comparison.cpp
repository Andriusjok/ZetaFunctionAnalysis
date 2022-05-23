#include "zeta.hpp"
#include <iostream>
#include <chrono>

template <typename Resolution = std::chrono::duration<double,std::micro>>
class Stopwatch {
  typedef std::chrono::steady_clock Clock;
  private:
    std::chrono::time_point<Clock> last;
  public:
    void reset() noexcept {
      last = Clock::now();
    }   
    Stopwatch() noexcept {
      reset();
    }   
    auto operator()() const noexcept {// returns time in Resolution
      return Resolution(Clock::now() - last).count() * 0.001;
    }
    ~Stopwatch() {
      std::cout << (*this)() << "\n";
    }
};

int main() {
    mpc_complex s = mpc_complex(0.5, 18457.135778615);
    int digits = 6;
    float accuracy = pow(10,-digits);
    {
    std::cout << "Zetafast:" << std::endl;
    Stopwatch sw;
    mpc_complex zetafastValue = zeta::zetafast(s,accuracy);
    std::cout << zetafastValue << std::endl;
    }
    {
    std::cout << "ZetaAMT:" << std::endl;
    Stopwatch sw;
    mpc_complex zetaamtct = zeta::zetaAMTCT(s,digits);
    std::cout << zetaamtct << std::endl;
    }
    {
    std::cout << "ZetaNAMB:" << std::endl;
    Stopwatch sw;
    mpc_complex zetanamb = zeta::zetaNA(s,digits, 1, 1);
    std::cout << zetanamb << std::endl;
    }
    {
    std::cout << "ZetaNABLC:" << std::endl;
    Stopwatch sw;
    mpc_complex zetanablc = zeta::zetaNA(s,digits, 1, 2);
    std::cout << zetanablc << std::endl;
    }
}