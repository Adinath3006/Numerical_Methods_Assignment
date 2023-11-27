# The Hermite Interpolation polynomial method.

# Refer Numerical Analysis by Richard L. Burden (9th Edition) Chapter 3 Section 4

#The argument A gets the input data matrix of the form: A = [[x1, x2, x3, ..., xn],
#                                                            [f(x1), f(x2), f(x3), ..., f(xn)],
#                                                            [f'(x1), f'(x2), f'(x3), ..., f'(xn)]]

import sympy as sym
import numpy as np

# Computing the Hermite Polynomial using the divide differnce method

# Intitializing z_i's and storing it in an array
def Modif_Dividediff(A,k,i):
    n = len(A[0])
    z = np.zeros((2*n,1))
    p = 0
    while p < n:
        z[2*p][0] = A[0][p]
        z[2*p+1][0] = A[0][p]
        p =p+1
    # The base case of iteration
    if k==0: 
        if i%2 == 0:
            f_diff = A[1][int(i/2)]
        else:
            f_diff = A[1][int((i-1)/2)]
    elif k==1:
        if i%2 == 0:
            f_diff = A[2][int(i/2)]
        else:
            f_diff = (Modif_Dividediff(A,k-1,i+1)-Modif_Dividediff(A,k-1,i))/(z[i+k][0]-z[i][0])
    else:
        f_diff = (Modif_Dividediff(A,k-1,i+1)-Modif_Dividediff(A,k-1,i))/(z[i+k][0]-z[i][0])
    return f_diff

def Hermite_Divided_Diff(A,val):
    n = len(A[0])
    z = np.zeros((2*n,1))
    i = 0
    while i < n:
        z[2*i][0] = A[0][i]
        z[2*i+1][0] = A[0][i]
        i =i+1
    x = sym.Symbol('x')
    i = 0
    m = 2*len(A)
    while i <m:
        if i == 0:
            Poly = Modif_Dividediff(A,i,0)
        else:
            j = 0
            while j<=i:
                if j == 0:
                    Ith_term = Modif_Dividediff(A,i,0)
                else:
                    Ith_term = Ith_term*(x-z[j-1][0])
                j = j+1
            Poly = Poly+Ith_term
        i = i+1
    print("-------------------------------------------------")
    print("Hermite polynomial computed using Divided difference method")
    print(sym.simplify(Poly))
    print(Poly.subs(x,val))
# 

def Lagrange_Coeff(A,i):
    x = sym.Symbol('x')
    n = len(A[0])
    m = 0
    Lag_coeff = 1
    while m<n:
        if i == m:
            Lag_coeff = Lag_coeff
        else:
            Lag_coeff = ((x-A[0][m])/(A[0][i]-A[0][m]))*Lag_coeff
        m = m+1
    return Lag_coeff

def H1_nj(A,j):
    x = sym.Symbol('x')
    return (1-2*(x-A[0][j])*((sym.diff(Lagrange_Coeff(A,j))).subs(x,A[0][j])))*((Lagrange_Coeff(A,j))**2)

def H2_nj(A,j):
    x = sym.Symbol('x')
    return (x-A[0][j])*((Lagrange_Coeff(A,j))**2)

def Hermite_Lagrange(A,val):
    Poly = 0
    i = 0
    x = sym.Symbol('x')
    n = len(A[0])-1
    while i <= n:
        Poly =  Poly + (A[1][i]*H1_nj(A,i)) + (A[2][i]*H2_nj(A,i))
        i = i+1
    print("-------------------------------------------------")
    print("Hermite polynomial computed using Lagrange Interpolation method")
    print(sym.simplify(Poly))
    print(Poly.subs(x,val))