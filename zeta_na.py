import math
from scipy.stats import norm

def c(k, k1, miu, sigma):
    if k < k1:
        return 1
    else:
        return 1 - norm.cdf((k - miu)/sigma)

def zeta_na(sigma, t, d, m, j):
    n = (math.pi/2)*t + (d + m)*math.log(10)
    z = norm.cdf(1-10**(-d))**-1
    if (j == 1):
        n = (n + math.log(2) - math.log(math.log(2)))/(math.log(3 + math.sqrt(8)))
        miu = n/math.sqrt(2)
        sigma = math.sqrt(n) * (32 ** 0.25)
    else:
        n = (n + math.log(2) - math.log(math.log(2)))/math.log(2)
        miu = n/2
        sigma = math.sqrt(2)/2
    k0 = math.ceil(miu + z*sigma)
    k1 = miu - z*sigma

    capital_s = 0 + 0j
    p= -1
    for k in range (0, math.floor(k0)):
        p = -p
        capital_s = capital_s + p*c(k, k1, miu, sigma)*math.exp(-sigma*math.log(k+1))*complex(math.cos(t*math.log(k+1)), -math.sin(t*math.log(k+1)))
    return capital_s / (1 - 2*math.exp(-sigma * math.log(2))) * complex(math.cos(t*math.log(2)), -math.sin(t*math.log(t* math.log(2))))