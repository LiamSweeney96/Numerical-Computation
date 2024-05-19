
def f(x):
    return x**2 - 100000000 * x + 1

def bisection_method(a, b, tol):
    iteration = 0
    
    if f(a)*f(b) > 0:
        #end function, no root.
        print("No root found.")
    else:
        while (b - a)/2.0 > tol:
            iteration += 1
            midpoint = (a + b)/2.0
            if f(midpoint) == 0:
                return(midpoint) #The midpoint is the x-intercept/root.
            elif f(a)*f(midpoint) < 0: # Increasing but below 0 case
                b = midpoint
            else:
                a = midpoint
                
        print("Iteration:", iteration)
        return(midpoint)

answer = bisection_method(0, 0.1, 1.e-8)

print("Answer:", answer)