import mpmath as mp

def compute_sum(f, tolerance):
    n=1
    computed_sum = 0.
    previous_sum = 0.
    while abs(computed_sum-previous_sum) > tolerance or n == 1:
        previous_sum = computed_sum
        computed_sum += f(n)
        n=n+1
    return computed_sum

def compute_complex_sum(f, tolerance):
    n=1
    computed_sum = mp.mpc(0, 0)
    previous_sum = mp.mpc(0, 0)
    while (abs(computed_sum-previous_sum).real > tolerance and abs(computed_sum-previous_sum).imag > tolerance) or n == 1:
        previous_sum = computed_sum
        computed_sum += f(n)
        n=n+1
    return computed_sum