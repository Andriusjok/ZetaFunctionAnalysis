from scipy import misc

def dx(f, x):
    return abs(0-f(x))

def newton_method(f, x0, e):
    delta = dx(f, x0)
    while delta > e:
        x0 = x0 - f(x0)/misc.derivative(f, x0)
        delta = dx(f, x0)
    return x0