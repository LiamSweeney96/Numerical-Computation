import numpy as np

def newton(f, df, x0, tol):
    x = x0
    y = f(x)
    it = 0
    while abs(y) > tol:   # iterate until less than or eq tol
        x = x - y / df(x)  # apply one Newton iteration
        y = f(x)           # reevaluate f at new estimate
        it = it + 1

    return x, it

# CHANGE ME
def f(t):
    return t*t - 100000000 * t + 1

# CHANGE ME
def df(t):
    return 2*t - 100000000

# CHANGE ME
x, it = newton(f, df, 0, 1.e-8)
print(f"Newton method: {x} after {it} iterations")
np.testing.assert_allclose(abs(f(x)), 0.0, atol=1.0e-8)


"""
    newton (f, df, START POINT, TOL (EG: 10^-4))
            
"""