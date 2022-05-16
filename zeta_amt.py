from sympy.ntheory.factor_ import totient
import math

def totatives(n):
    phi = int(n > 1 and n)
    for p in range(2, int(n ** .5) + 1):
        if not n % p:
            phi -= phi // p
            while not n % p:
                n //= p
    #if n is > 1 it means it is prime
    if n > 1: phi -= phi // n 
    return phi

def omega(sigma, t):
    return 1 - 2**(2-sigma)*math.cos(t*math.log(2))+4**(1-sigma)

def C(n, k):
    argument = (8**(1/4))*(k*math.sqrt(2)-n)/math.sqrt(n)
    if (argument < 1):
        return 1
    return 1 - totatives(math.floor(argument))

def zeta_amt(sigma, t, d):
    n = 1.3*d + 0.9*abs(t)
    ReZ = 0
    ImZ = 0
    for k in range(0,math.floor(n-1)):
        a = (-1)**k * C(n, k) / ((k + 1)**sigma)
        b1 = t*math.log(k+1)
        b2 = t*math.log((k+1)/2)
        ReZ += a*(math.cos(b1)-(2**(1-sigma))*math.cos(b2))
        ImZ += a*(math.sin(b1)-(2**(1-sigma))*math.sin(b2))
    return complex(ReZ, -ImZ)/omega(sigma, t)