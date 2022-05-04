import math
from newton import newton_method
import scipy.special
from sum import compute_sum

accuracy = 0.05

def q_gamma(v, m):
    w = 0
    w_limit = v-1
    sum = 0
    while w < w_limit:
        sum += m**w / math.factorial(w) * (math.e ** -m)
        w+=1
    return sum

def resolve_capital_n(v, tau):
    return 1.11 * ((1 + ((0.5 + tau) / v)) ** 0.5)

def resolve_v(sigma, tau):
    def f(x):
        return (x - max((1-sigma)/2, 0) * math.log(1/2 + x + tau))/math.log(8/accuracy) - 1
    return newton_method(f, 1, accuracy)

def d_series(s, v, N):
    def f(n):
        return 1/(n**s) * q_gamma(v, n / N)
    return compute_sum(f)

def e_1_series(s, N, sign):
    def f(m):
        returne_m_sum += e_miu_1_series(m, s, N, sign)
    compute_sum(f)
    return ((2*math.pi) ** s-1) * scipy.special.gamma(1 - s) * (math.e ** (sign*(math.pi/2)*(1-s)))

def e_miu_1_series(m, s, N, sign):
    def f(w):
        return (scipy.special.gamma(s) / (math.factorial(w)*scipy.special.gamma((s-1)-w + 1))) * ((m + (sign/(2* math.pi * N)))**s-1-w) * ((-sign/(2*math.pi*N))**w)
    return (m**(s - 1)) - compute_sum(f)

def gamma(s, v, N):
    return (scipy.special.gamma(1-s + v) / ((1 - s) * math.gamma(v))) * N**(1-s)

def zetafast():
    s = 0.5 + 14.134725141734693790457251983562j
    v = resolve_v(s.real, s.imag)
    #print(v)
    N = resolve_capital_n(v, s.imag)
    #print (N)
    d_series_value = d_series(s, v, N)
    #print (d_series_value)
    e_series_positive = e_1_series(s, N, 1)
    #print(e_series_positive)
    e_series_negative =  e_1_series(s, N, -1)
    #print(e_series_negative)
    gamma_value = gamma(s, v, N)
    #print(gamma_value)
    zeta = d_series_value + e_series_positive + e_series_negative - gamma_value
    s= s+0.0005j
    print(s)
    print(zeta)

zetafast()