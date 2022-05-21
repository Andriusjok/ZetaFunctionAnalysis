import mpmath as mp
from statistics import NormalDist

def omega(sigma, t):
    return 1 - (2**(2-sigma))*mp.cos(t*mp.log(2))+4**(1-sigma)

def C(k0, k1, miu_n, sigma_n, k):
    if (k0 < k1):
        return 1
    return 1 - NormalDist().cdf((k - miu_n)/sigma_n)

def zeta_amtct(sigma, t, d):
    sigma = mp.mpf(sigma)
    t = mp.mpf(t)
    n = mp.pi * abs(t) + (1 + 2 * sigma)*mp.log(abs(t)) + 2*mp.log(mp.gamma(sigma)) - 2*mp.log(abs( 1 - 2 ** (1-sigma))) + 2*d*mp.log(10) + mp.log(2) - mp.log(mp.pi)
    n = n/(2*mp.log(3+mp.sqrt(8)))
    epsilon = 10 ** -d
    miu_n = n/mp.sqrt(2)
    sigma_n = mp.sqrt(n)/32**0.25
    z = NormalDist().inv_cdf(1-epsilon)
    k0 = mp.ceil(miu_n + (z*sigma_n))
    k1 = mp.floor(miu_n - z*sigma_n)
    ReZ = 0
    ImZ = 0
    k = 0
    while k <= k0:
        a = (-1)**k * C(k0, k1, miu_n, sigma_n, k) / ((k + 1)**sigma)
        b1 = t*mp.log(k+1)
        b2 = t*mp.log((k+1)/2)
        ReZ += a*(mp.cos(b1)-(2**(1-sigma))*mp.cos(b2))
        ImZ += a*(mp.sin(b1)-(2**(1-sigma))*mp.sin(b2))
        k+=1
    return mp.mpc(ReZ, -ImZ)/omega(sigma, t)