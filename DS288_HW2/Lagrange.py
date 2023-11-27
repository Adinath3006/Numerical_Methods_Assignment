# The Lagrange Interpolation polynomial method.

# Refer Numerical Analysis by Richard L. Burden (9th Edition) Chapter 3 Section 1

#The argument A gets the input data matrix of the form: A = [[x1, x2, x3, ..., xn],[f(x1), f(x2), f(x3), ..., f(xn)]]

import sympy as sym

def Lagrange(A,val):
    n = len(A[0]) # The total number of nodes 
    i = 0 # Iterating variable
    P = [0 for x in range(n)] # Creating an array to store the terms of the polynomial
    Poly = 0 # Intitailizing the interpolation polynomial
    x = sym.Symbol('x')
    while i<n:
       # Computing the Lagrange coefficient 
        j = 0 # Iterating variable 
        Lag_coeff = 1 # Intitailizing the Lagrange coefficient
        while j<n:
            if i == j:
                Lag_coeff = Lag_coeff
            else:
                Lag_coeff = ((x-A[0][j])/(A[0][i]-A[0][j]))*Lag_coeff
            j = j+1
        # Computing the ith term of the interpolating polynomial and storing it in an array
        P[i] = (Lag_coeff)*(A[1][i])
        i = i+1     
    m = 0 # Iterating variable
    while m<n:
        # Sum the ith terms to obtain the desired interpolation polynomial 
        Poly = Poly+P[m] 
        m = m+1
    print("-------------------------------------------------")
    print("Lagrange Interpolation Polynomial")
    print(sym.N(sym.simplify(Poly),4))
    print(sym.N(Poly.subs(x,val),8))

