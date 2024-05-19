from pprint import pprint
import numpy as np

def jacobi(A, b, n=3, x=None):
    """Solves the equation Ax=b via the Jacobi iterative method."""
    # Create an initial guess if needed                                                                                                                                                            
    if x is None:
        x = np.zeros(len(A[0]))

    # Create a vector of the diagonal elements of A                                                                                                                                                
    # and subtract them from A                                                                                                                                                                     
    D = np.diag(A)
    R = A - np.diagflat(D)

    # Iterate for N times                                                                                                                                                                          
    for i in range(n):
        x = (b - np.dot(R,x)) / D
    return x

A = np.array([[2, 1],[1, 2]])
b = np.array([3, 3])
guess = np.array([0, 0])

sol = jacobi(A, b, n=3, x=guess)

print ("A:")
pprint(A)

print ("b:")
pprint(b)

print ("x:")
pprint(sol)


print("Residual: ", np.matmul(A,guess)-b)