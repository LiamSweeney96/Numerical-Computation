def number_systems(b, t, l, u):
    
    mantissa = [9] * t
    total = 0
    count = 1
    
    for i in range(t):
        if i == 0:
            total += mantissa[i]
        else:
            total *= 10
    
    for i in range(max(u-l,l-u)):
        count += 1
        
    total = total * count * t
    
    return total + 2

def largest(b, t, l, u):
    
    array = [9] * t
    join = ''.join(str(x) for x in array)
    return float(join)

def smallest():
    num = 0.1
    return num * (10**l)
    

b = 10
t = 3
l = -3
u = 3

x = number_systems(b, t, l, u)
y = largest(b, t, l, u)
z = smallest()
             
print ('Total numbers:', x)
print ('Largest number:', y)
print ('Smallest number:', z)
print ('Second smallest number', z+(z/10))
print ('Smallest Difference:', z/10)