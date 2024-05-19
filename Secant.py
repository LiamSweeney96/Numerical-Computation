import numpy as np

def secant(f, x0, x1, tol):
    x = x1
    it = 0
    while abs(f(x)) > tol:   # iterate until less than or eq tol
        x = x - f(x1) *(x1-x0) / (f(x1) - f(x0))  # apply one Newton iteration
        x0 = x1
        x1 = x
        it = it + 1

    return x, it

def f(x):
    return x*x*x - 6.0*x*x +9.0*x


x, it = secant(f, 4.0, 5.0, 1.e-6)
print(f"The secant method: {x} after {it} iterations")
np.testing.assert_allclose(abs(f(x)), 0.0, atol=1.0e-6)
