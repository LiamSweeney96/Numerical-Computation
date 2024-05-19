import numpy as np

def lower_triangular_solve(A, b):
    """
    Solve the system  A x = b  where A is assumed to be lower triangular,
    i.e. A(i,j) = 0 for j > i, and the diagonal is assumed to be nonzero,
    i.e. A(i,i) != 0.
    
    The code checks that A is lower triangular and converts A and b to
    double precision before computing.

    ARGUMENTS:  A   lower triangular n x n array
                b   right hand side column n-vector

    RETURNS:    x   column n-vector solution
    """

    # we should take care to ensure that arrays are stored with the correct type - float!
    A = A.astype(np.float64)
    b = b.astype(np.float64)
     
    # check sizes of A and b match appropriately
    nb=len(b)
    n, m = A.shape
    if n != m:
        raise ValueError(f'A is not a square matrix!')
    if n != nb:
        raise ValueError(f'shapes of A and b do not match!')
    
    # checks whether A is lower triangular
    for i in range(n):
        for j in range(i+1,n):
            if not np.isclose(A[i, j], 0.0):
                raise ValueError(f'A is not lower triangular!')

    # checks whether A has zero diagonal element
    for i in range(n):
        if np.isclose(A[i, i], 0.0):
            raise ValueError(f'A[{i}, {i}] is zero')
    
    # create a new array to store the results
    x = np.empty_like(b)
    
    # perform forward substitution
    x[0] = b[0] / A[0, 0]
    for i in range(1,n):
        x[i] = b[i] / A[i, i]
        for j in range(i):
            x[i] = x[i] - A[i,j]*x[j]/A[i, i]
        
    return x

def Jacobi_iteration(A, b, max_iteration, x0 = None):
    # we should take care to ensure that arrays are stored with the correct type - float!
    A = A.astype(np.float64)
    b = b.astype(np.float64)
     
    # check sizes of A and b match appropriately
    nb=len(b)
    n, m = A.shape
    if n != m:
        raise ValueError(f'A is not a square matrix!')
    if n != nb:
        raise ValueError(f'shapes of A and b do not match!')

    # check diagonal is non zero
    for i in range(n):
        if np.isclose(A[i, i], 0):
            raise ValueError(f'A[{i}, {i}] is zero')

    # construct iteration matrices
    P=np.zeros([n,n])    # matrix P = D^{-1}(L+U)
    p=np.zeros(n)        # vector p = D^{-1} b
    for i in range(n):
        p[i]=b[i]/A[i,i] 
        for j in range(n):
             P[i,j] = A[i,j]/A[i,i]
        P[i,i] = 0
        
    #create a new array to store the results, initialised as zero
    if x0 is None:
        x = np.zeros_like(b)
    else:
        x = x0.copy()
    
    # perform iteration x <- p - P * x
    for it in range(max_iteration):
        xnew = np.empty_like(x)
        for i in range(n):
            xnew[i] = p[i]
            for j in range(n):
                xnew[i] -= P[i, j] * x[j]
        x = xnew.copy()
                
    return x

def Gauss_Seidel_iteration(A, b, max_iteration, x0 = None):
    # we should take care to ensure that arrays are stored with the correct type - float!
    A = A.astype(np.float64)
    b = b.astype(np.float64)
     
    # check sizes of A and b match appropriately
    nb=len(b)
    n, m = A.shape
    if n != m:
        raise ValueError(f'A is not a square matrix!')
    if n != nb:
        raise ValueError(f'shapes of A and b do not match! ')

    for i in range(n):
        if np.isclose(A[i, i], 0):
            raise ValueError(f'A[{i}, {i}] is zero')

    # do not construct iteration matrices explicitly
    LD = np.zeros_like(A)
    U = np.zeros_like(A)
    for i in range(n):
        for j in range(n):
            if i < j:
                U[i, j] = A[i, j]
            else:
                LD[i, j] = A[i, j]
    
    # p = (L + D)^{-1} b --> found by solving triangular system
    # (L + D) p = b
    p = lower_triangular_solve(LD, b)
      
    #create a new array to store the results, initialised as zero
    if x0 is None:
        x = np.zeros_like(b)
    else:
        x = x0.copy()
        
    # perform iteration x <- p - P * x
    # (L+D)(xnew - p) = U*x
    Ux = np.empty_like(x)
    for it in range(max_iteration):
        for i in range(n):
            Ux[i] = 0.0
            for j in range(i+1, n):
                Ux[i] += U[i, j] * x[j]
        Px = lower_triangular_solve(LD, Ux)
        x = p - Px
                
    return x

A = np.array([[2, -1, 0], [-1, 1, 3], [1, 0, 1]])
b = np.array([1, 3, 2])
x_exact = np.array([0, 0, 0])


# numpy linear solver
x0 = np.linalg.solve(A,b)
print("Solution by numpy solver:", x0)

print('\n')

x = Jacobi_iteration(A, b, 3)
print("Solution by Jacobi iteration: ",x)
print("Absolute Error:", x_exact - x)

print("Residual: ", np.matmul(A,x)-b)

print('\n')

x = Gauss_Seidel_iteration(A, b, 1)
print("Solution by Gauss Seidel iteration: ",x)
print("Absolute Error:", x_exact - x)
print("Residual: ", np.matmul(A,x)-b)


"""

    ERROR:
    
    SQRT((EXACT SOLUTION - APPROX SOLUTION)**2 + (EXACT SOLUTION - APPROX SOLUTION)**2)

    RESIDUAL:
        
    1. CONVERT CODE RESULT RESIDUAL TO FRACTION
    2. SQUARE THE FRACTIONS
    3. ADD THE FRACTIONS
    4. SQUARE ROOT THE RESULT
    
    
    (1/2, 7/2, 3/2) = (-3.5, 4.5, 0)
    
    SQRT( (7/2)^2 + (9/2)^2 + (0)^2 ) = SQRT ( 130/4 )
    
    RESIDUAL = 5.701 (3 D.P.)

"""

