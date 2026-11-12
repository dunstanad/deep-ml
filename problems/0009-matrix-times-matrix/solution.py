def matrixmul(a:list[list[int|float]],
              b:list[list[int|float]])-> list[list[int|float]]:
    
    ar = len(a)
    ac = len(a[0])

    br = len(b)
    bc = len(b[0])

    if ac != br:
        return -1
    
    c = [ [0 for _ in range(bc)] for _ in range(ar)]

    for i in range(ar):
        for j in range(bc):
            for k in range(ac):
                c[i][j] += a[i][k] * b[k][j]
                        
    
    return c