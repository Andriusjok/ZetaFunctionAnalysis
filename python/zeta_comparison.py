import zetafast
import zeta_amtct
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
#s = mp.mpc(50000, 20000)
s = mp.mpc(1.5, 500)
#d
d = 6
accuracy = 10**-d
#m
m = 1

start = time.perf_counter()
zetafast_value = zetafast.zetafast(s, accuracy)
end = time.perf_counter()
zetafast_time = end - start

start = time.perf_counter()
zeta_amtct_value = zeta_amtct.zeta_amtct(s.real, s.imag, d)
end = time.perf_counter()
zeta_amtct_time = end - start

start = time.perf_counter()
zeta_namb_value = zeta_na.zeta_na(s.real, s.imag, d, m, 1)
end = time.perf_counter()
zeta_namb_time = end - start

start = time.perf_counter()
zeta_nablc_value = zeta_na.zeta_na(s.real, s.imag, d, m, 2)
end = time.perf_counter()
zeta_nablc_time = end - start

start = time.perf_counter()
mpmath_value = mp.zeta(s)
end = time.perf_counter()
mpath_time = end - start

print("Algorithm:{} executed in {} with value: {}".format("zetafast", zetafast_time, zetafast_value))
print("Algorithm:{} executed in {} with value: {}".format("zeta_amtct", zeta_amtct_time, zeta_amtct_value))
print("Algorithm:{} executed in {} with value: {}".format("zeta_namb", zeta_namb_time, zeta_namb_value))
print("Algorithm:{} executed in {} with value: {}".format("zeta_nablc", zeta_nablc_time, zeta_nablc_value))
print("Algorithm:{} executed in {} with value: {}".format("mpmath", mpath_time, mpmath_value))