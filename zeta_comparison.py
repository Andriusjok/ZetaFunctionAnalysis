import zetafast
import zeta_amt
import zeta_na
import time
 
#s = 0.5 + 9.134725141734693790457251983562j
#s = 0.5 + 21.022039638771554992628479593896j
#s = 0.5 + 25.010857580145688763213790992562j
#s = 2 + 0j
s =  0.5 + 18457.135778615j
start = time.perf_counter()
zetafast_value = zetafast.zetafast(s)
end = time.perf_counter()
zetafast_time = end - start

#d
d = 2
#m
m = 2

start = time.perf_counter()
zeta_amt_value = zeta_amt.zeta_amt(s.real, s.imag, d)
end = time.perf_counter()
zeta_amt_time = end - start

start = time.perf_counter()
zeta_namb_value = zeta_na.zeta_na(s.real, s.imag, d, m, 1)
end = time.perf_counter()
zeta_namb_time = end - start

start = time.perf_counter()
zeta_nablc_value = zeta_na.zeta_na(s.real, s.imag, d, m, 2)
end = time.perf_counter()
zeta_nablc_time = end - start

print("Algorithm:{} executed in {} with value: {}".format("zetafast", zetafast_time, zetafast_value))
print("Algorithm:{} executed in {} with value: {}".format("zeta_amt", zeta_amt_time, zeta_amt_value))
print("Algorithm:{} executed in {} with value: {}".format("zeta_namb", zeta_namb_time, zeta_namb_value))
print("Algorithm:{} executed in {} with value: {}".format("zeta_nablc", zeta_nablc_time, zeta_nablc_value))