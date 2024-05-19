import numpy as np

eps = np.finfo(float).eps
dx = np.sqrt(eps)

x0 = 1.0
df_approx = ((x0+dx)**3 - x0**3) / dx
abs_error = abs(df_approx - 3)
rel_error = abs_error / 3

print("dx =", dx)
print("df_approx =", df_approx)
print("abs_error =", abs_error)
print("rel_error =", rel_error)
