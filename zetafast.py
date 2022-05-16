import math
from sum import compute_sum
from scipy.optimize import fsolve
from mpmath import mp 

accuracy = 10e-3   
lambda_var = 3.151

def q_gamma(v, m):
    w = 0
    w_limit = v-1
    sum = 0
    while w < w_limit:
        sum += m**w / mp.factorial(w) * (math.e ** -m)
        w+=1
    return sum

def resolve_capital_n(v, tau):
    return 1.11 * ((1 + ((0.5 + tau) / v)) ** 0.5)

def resolve_v(sigma, tau):
    def f(x):
        return x - (max((1-sigma)/2, 0) * math.log(1/2 + x + tau)) - math.log(8/accuracy)
    return math.ceil(fsolve(f, 1))

def d_series(s, v, N):
    def f(n):
        return n**-s * q_gamma(v, n / N)
    return compute_sum(f)

def d_series_lambda(s, v, N):
    def f(n):
        return n**-s * q_gamma(v, n / N)
    n = 1
    sum = 0
    while n < math.ceil(lambda_var*v*N):
        sum+=f(n)
        n+=1
    return sum

def complex_exponent(exponent):
    return mp.mpf(math.e)**mp.mpf(exponent.real)*mp.mpc(complex(math.cos(exponent.imag),math.sin(exponent.imag)))

def e_1_series(s, N, sign):
    def f(m):
        e_miu_1_series(m, s, N, sign)
    return ((2*math.pi) ** s-1) * mp.gamma(1 - s) * complex_exponent(complex(0,sign)*(math.pi/2)*(1-s)) * compute_sum(f)

def e_miu_1_series(m, s, N, sign):
    def f(w):
        return mp.binomial(s-1, w) * ((m + (sign/(2* math.pi * N)))**s-1-w) * ((complex(0,-sign)/(2*math.pi*N))**w)
    return (m**(s - 1)) - compute_sum(f)

def e_1_series_m_s(capital_M, s):
    sum = 0
    m = 1
    while capital_M <= m:
        s+=e_1_series_m_s(m, s)
        m+=1
    return ((2*math.pi)**(s-1))*mp.gamma(1-s)*complex_exponent(complex(0,math.pi/2)*(1-s))*sum

def gamma(s, v, N):
    return (mp.gamma(1-s + v) / ((1 - s) * mp.gamma(v))) * N**(1-s)

def zetafast(s):
    v = resolve_v(s.real, s.imag)
    N = resolve_capital_n(v, s.imag)
    #if 0 <= s.real and s.real <= 2 and s.imag > 0 and accuracy <= 0.05:
    #    c = 0.5
    #    return d_series_lambda(s,v,N) + e_1_series_m_s(math.ceil(N),s) - gamma(s, v, N) + c*accuracy
    #else:
    if s.imag == 0:
        c = 0.5
        return d_series_lambda(s,v,N) - gamma(s,v,N) + c*accuracy
    else:
        return d_series(s,v,N) + e_1_series(s, N, 1) + e_1_series(s, N, -1) - gamma(s, v, N)