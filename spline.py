import numpy as np
from collections import Counter

def knot_insert(knots, coeffs, d, new_knot):
    """
    Return: new_knots, new_coeff
    """

    j = next( (idx for idx,knot in enumerate(knots) if knot > new_knot), len(knots))
    #print(j)

    new_coeff = []

    for i in range(j-d + 1):
        #print(i, knots[i:i+d])
        new_coeff.append(coeffs[i])

    for i in range(j - d, j):
        t = (new_knot - knots[i])/(knots[i+d] - knots[i])
        new_coeff.append(coeffs[i]*(1-t) + coeffs[i+1]*t)

    for i in range(j, len(knots) - d + 1):
        #print(i, knots[i:i+d])
        new_coeff.append(coeffs[i])
        

    new_knots = list(knots)
    new_knots.insert(j, new_knot)
    return new_knots, new_coeff

def nurbs2bezier(knots, coeffs, d):
    """
    Return: bezier coeffs, no duplicated coeffs
    """

    # check knots multiplicity
    multi = Counter()
    for i in knots:
        multi[i] += 1
    # print(multi)

    assert(multi.most_common(1)[0][1] <= d)

    new_knots = list(knots)
    new_coeff = list(coeffs)
    
    for knot, m in multi.items():
        for i in range(d - m):
            new_knots, new_coeff = knot_insert(new_knots, new_coeff, d , knot)

    return new_knots, new_coeff
    
    
    

if __name__ == '__main__':
    #knots = [-4, -4, -4, -3,-2,-1,0,1,2,3, 4, 4, 4]
    #coeff = [0,0, 0, 0,36,0,0,0,0]
    
    knots = [-1, -1, -1, 0, 1,2,3,4 , 5,5,5]
    coeff = [0,0,0, 0, 6,0,0,0 , 0]

    print(nurbs2bezier(knots, coeff, 3))
    
    # for i in range(3):
    #     knots, coeff = knot_insert(knots, coeff, 3 , 2)
    #     print(knots, coeff)
