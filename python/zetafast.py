import math
from sum import compute_sum, compute_complex_sum
from mpmath import mp 
from mpmath.calculus.optimization import Newton

lambda_var = 3.151

def q_gamma(v, m):
    w = 0
    sum = 0
    while w <= v-1:
        sum += m**w / mp.factorial(w) * (math.e ** -m)
        w+=1
    return sum

def resolve_capital_n(v, tau):
    return 1.11 * (1 + ((0.5 + tau) / v)) ** 0.5

def resolve_v(sigma, tau, accuracy):
    def f(x):
        return x - max((1-sigma)/2, 0) * math.log(1/2 + x + tau) - math.log(8/accuracy)
    if sigma >= 1:
        return mp.ceil(math.log(8/accuracy))
    newton = next(Newton(mp, f= f, x0 = [5]).__iter__())
    return mp.ceil(newton[0])

def d_series(s, v, N, accuracy):
    def f(n):
        return n**-s * q_gamma(v, n / N)
    return compute_sum(f, accuracy)

def d_series_lambda(s, v, N):
    def f(n):
        return n**-s * q_gamma(v, n / N)
    n = 1
    sum = 0
    while n <= math.ceil(lambda_var*v*N):
        sum+=f(n)
        n+=1
    return sum

def complex_exponent(exponent):
    return mp.mpf(math.e)**mp.mpf(exponent.real)*mp.mpc(mp.mpc(math.cos(exponent.imag),math.sin(exponent.imag)))

def e_1_series(s, N, sign, v, accuracy):
    def f(m):
        return  e_miu_1_series(m, s, N, sign, v)
    return ((2*math.pi) ** s-1) * mp.gamma(1 - s) * complex_exponent(mp.mpc(0,sign * math.pi/2)*(1-s)) * compute_complex_sum(f, accuracy)

def e_miu_1_series(m, s, N, sign, v):
    sum = 0
    w = 0
    while w <= (v-1):
        sum += mp.binomial(s-1, w) * ((m + (mp.mpc(0,sign)/(2 * math.pi * N))) ** (s-1-w)) * ((mp.mpc(0,-sign)/(2*math.pi*N))**w)
        w+=1
    return (m**(s - 1)) - sum

def e_1_series_m_s(capital_M, s, N, v):
    sum = 0
    m = 1
    while m <= capital_M:
        sum+=e_miu_1_series(m, s, N, 1, v)
        m = m+1
    return ((2*math.pi)**(s-1))*mp.gamma(1-s)*complex_exponent(mp.mpc(0,math.pi/2)*(1-s))*sum

def gamma(s, v, N):
    return (mp.gamma(1-s + v) / ((1 - s) * mp.gamma(v))) * N**(1-s)

def zetafast(s, accuracy):
    v = resolve_v(s.real, s.imag, accuracy)
    N = resolve_capital_n(v, s.imag)
    if 0 <= s.real and s.real <= 2 and s.imag > 0 and accuracy <= 0.05:
        return d_series_lambda(s,v,N) + e_1_series_m_s(mp.ceil(N),s, N, v) - gamma(s, v, N)
    else:
        if s.imag == 0:
            return d_series_lambda(s,v,N) - gamma(s,v,N)
        else:
            return d_series(s,v,N, accuracy) + e_1_series(s, N, 1, v, accuracy) + e_1_series(s, N, -1, v, accuracy) - gamma(s, v, N)