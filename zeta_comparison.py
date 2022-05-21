import zetafast
import zeta_amtct
import zeta_amt
import zeta_na
import time
import mpmath as mp
 
#s = 0.5 + 9.134725141734693790457251983562j
#s = 0.5 + 21.022039638771554992628479593896j
#s = 0.5 + 25.010857580145688763213790992562j
#s = 2 + 0j
#s =  0.5 + 18457.135778615j
#s = 5 + 2j

#s = mp.mpc(3, 1000)
s = mp.mpc(1.19, 1001)
#d
d = 10e-3
#m
m = 2

start = time.perf_counter()
zetafast_value = zetafast.zetafast(s, d)
end = time.perf_counter()
zetafast_time = end - start

start = time.perf_counter()
zeta_amtct_value = zeta_amtct.zeta_amtct(s.real, s.imag, d)
end = time.perf_counter()
zeta_amtct_time = end - start

start = time.perf_counter()
zeta_amt_value = zeta_amt.zeta_amt(s.real, s.imag, d)
end = time.perf_counter()
zeta_amt_time = end - start

start = time.perf_counter()
mpmath_value = mp.zeta(s)
end = time.perf_counter()
mpath_time = end - start

print("Algorithm:{} executed in {} with value: {}".format("zetafast", zetafast_time, zetafast_value))
print("Algorithm:{} executed in {} with value: {}".format("zeta_amtct", zeta_amtct_time, zeta_amtct_value))
print("Algorithm:{} executed in {} with value: {}".format("zeta_amt", zeta_amt_time, zeta_amt_value))
print("Algorithm:{} executed in {} with value: {}".format("mpmath", mpath_time, mpmath_value))