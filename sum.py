def compute_sum(f, tolerance=1e-5):
    n=1.
    computed_sum = 0.
    previous_sum = 0.

    while ((computed_sum - previous_sum).real > tolerance or (computed_sum - previous_sum).imag > tolerance):
        previous_sum = computed_sum
        computed_sum += f(n)
        n=n+1.
    return computed_sum