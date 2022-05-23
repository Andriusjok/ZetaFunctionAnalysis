import mpmath as mp
from statistics import NormalDist

def c(k, k1, miu, sigma_n):
    if k < k1:
        return 1
    else:
        return 1 - NormalDist().cdf((k - miu)/sigma_n)

def zeta_na(sigma, t, d, m, j):
    n = (mp.pi/2)*t + (d + m)*mp.log(10)
    z = NormalDist().inv_cdf(1-10**-d)
    if (j == 1):
        n = (n + mp.log(2) - mp.log(mp.log(2))) / (mp.log(3 + mp.sqrt(8)))
        miu = n/mp.sqrt(2)
        sigma_n = mp.sqrt(n) / (32 ** 0.25)
    else:
        n = (n + mp.log(2) - mp.log(mp.log(2))) / mp.log(2)
        miu = n/2
        sigma_n = mp.sqrt(n)/2
    k0 = mp.ceil(miu + z*sigma_n)
    k1 = miu - z*sigma_n
    capital_s = 0
    p = -1
    k = 0
    while k <= k0:
        p = -p
        capital_s += p*c(k, k1, miu, sigma_n) * mp.exp(-sigma*mp.log(k+1)) * mp.mpc(mp.cos(t*mp.log(k+1)), -mp.sin(t*mp.log(k+1)))
        k+=1
    return capital_s / (1 - 2 * mp.exp(-sigma * mp.log(2)) * mp.mpc(mp.cos(t * mp.log(2)), -mp.sin(t * mp.log(2))))