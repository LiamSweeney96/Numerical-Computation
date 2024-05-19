def best_fit (next_value):
    
    """
        
    UPDATE ARRAYS
    
    """
    
    x = [-1, 0, 1]
    y = [0.5, 1.2, 1.5]
    
    meanX = 0
    meanY = 0
    countX = 0
    countY = 0
    
    for i in range(len(x)):
        meanX += x[i]
        countX += 1
    
    for i in range(len(y)):
        meanY += y[i]
        countY += 1
    
    meanX = meanX / countX
    meanY = meanY / countY
    
    sumSquares = 0
    sumProducts = 0
    
    # (X - M_x)^2
    for i in range(len(x)):
        sumSquares += (x[i] - meanX)**2
    
    # (X - M_x)(Y - M_y)
    for i in range(len(x)):
        sumProducts += (x[i] - meanX) * (y[i] - meanY)
        
    b = sumProducts / sumSquares
    
    a = meanY - (b * meanX)
    
    newY = (b * next_value) + a
    
    print (newY)

best_fit(4)    
    
    