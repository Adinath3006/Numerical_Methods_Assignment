# The set of function in this code computes the interpolation polynomial P_n(x) is the Lagrange polynomial that agrees with the 
# function 'f' at the points x1, x2, x3, ..., xn (nodes), when xi's are equally spaced with each other we can use FORWARD DIFFERENCE METHOD
# to represent the polynomial in the form
#                                               x_(i+1) - x(i) = h   where i = 0,...,n-1
#                                                       Let      x = x_(0) + s*h
#
#                       P_n(x) = f[x0] + (sC1)*(1!)*(h^1)*f[x_(0),x_(1)] + ... + (sCn)*(n!)*(h^n)*f[x_(0),....,x_(n)]
#                                                   
#                                         where sCi is the combination formula --> (s)*(s-1)*...*(s-i+1)/(i!)
#                       
#                                      P_n(x) = f[x0] + (sC1)*Delta_f^1(x0) + ... + (sCn)*Delta_f^n(x0)
#                                                                                                          
#                                                                                                         Delta_f^i(x0)
# Where the ai's are the constants found bu computing the ith divided difference: f[x_(0),....,x_(0+i)] = -------------
#                                                                                                          i!*(h^i)  


#The argument A gets the input data matrix of the form: A = [[x1, x2, x3, ..., xn],[f(x1), f(x2), f(x3), ..., f(xn)]]

import sympy as sym
import math

#                                                                                       ( f[x_(i+1),....,x_(i+k)] - f[x_(i),....,x_(i+k-1)] )
# Defining a function to compute the kth divided difference --> f[x_(i),....,x_(i+k)] = -----------------------------------------------------
#                                                                                                         x_(i+k) - x_(i)

def Dividediff(A,k,i):
    # The base case of iteration
    if k==0: 
        f_diff = A[1][i]
    else:
        f_diff = (Dividediff(A,k-1,i+1)-Dividediff(A,k-1,i))/(A[0][i+k]-A[0][i])
    return f_diff


# Defining a function to compute the forward difference Delta_f^k(x0)

def ForwardDiff_Poly(A,val):
    n = len(A[0])
    s = sym.symbols('s')
    i = 1
    Poly = A[1][0]
    K_term = 1
    h = A[0][1]-A[0][0]
    while i < n:
        K_term = K_term*(s+1-i)/i
        K_coeff = Dividediff(A,i,0)*(sym.factorial(i))*(h**i)
        Poly = Poly + K_coeff*K_term
        i = i+1
    print("-------------------------------------------------")
    print("Newton's Forward Difference Method")
    print(Poly)
    print(Poly.subs(s,val))