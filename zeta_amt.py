from statistics import NormalDist
import mpmath as mp

def omega(sigma, t):
    return 1 - 2**(2-sigma)*mp.cos(t*mp.log(2))+4**(1-sigma)

def C(n, k):
    return 1 - NormalDist().cdf(mp.floor((8**(1/4))*(k*mp.sqrt(2)-n)/mp.sqrt(n)))

def zeta_amt(sigma, t, d):
    n = 1.3*d + 0.9*abs(t)
    ReZ = 0
    ImZ = 0
    k = 0
    while k <= mp.floor(n-1):
        a = (-1)**k * C(n, k) / ((k + 1)**sigma)
        b1 = t*mp.log(k+1)
        b2 = t*mp.log((k+1)/2)
        ReZ += a*(mp.cos(b1)-(2**(1-sigma))*mp.cos(b2))
        ImZ += a*(mp.sin(b1)-(2**(1-sigma))*mp.sin(b2))
        k+=1
    return mp.mpc(ReZ, -ImZ)/omega(sigma, t)