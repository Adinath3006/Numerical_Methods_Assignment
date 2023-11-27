import numpy as np
import sympy as sym

x = sym.symbols('x')

def MatrixA(coeff,a,b,n):
    h = (b-a)/(n+1)
    i = 0
    A = np.zeros((n, n))
    while i<n:
        j = 0
        while j<n:
            if abs(i-j) > 1:
                A[i][j] = 0
            else:
                if i > j:
                    A[i][j] = -1 - (h/2)*(coeff[0].subs(x,a+(i+1)*h))
                elif j > i:
                    A[i][j] = -1 + (h/2)*(coeff[0].subs(x,a+(i+1)*h))
                else:
                    A[i][j] = 2 + (h**2)*(coeff[1].subs(x,a+(i+1)*h))
            j = j+1
        i =i+1
    return A

def Matrixb(coeff,a,b,n,init):
    h = (b-a)/(n+1)
    i = 0
    B = np.zeros((n, 1))
    while i<n:
        if i == 0:
            B[i][0] = -(h**2)*(coeff[2].subs(x,a+(i+1)*h)) + (1 + (h/2)*(coeff[0]).subs(x,a+(i+1)*h))*(init[0])
        elif i == n-1:
            B[i][0] = -(h**2)*(coeff[2].subs(x,a+(i+1)*h)) + (1 - (h/2)*(coeff[0]).subs(x,a+(i+1)*h))*(init[1])
        else:
            B[i][0] = -(h**2)*(coeff[2].subs(x,a+(i+1)*h))
        i = i+1
    return B

def MatrixW(coeff,a,b,n,init):
    soln = np.linalg.inv(MatrixA(coeff,a,b,n)).dot(Matrixb(coeff,a,b,n,init))
    print("")
    print("-------------------------------------------------")
    print("             Finite-Difference Method            ")
    print("-------------------------------------------------")
    print("")
    print(soln)
    return soln
