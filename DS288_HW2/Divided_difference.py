# The set of function in this code computes the interpolation polynomial P_n(x) is the Lagrange polynomial that agrees with the 
# function 'f' at the points x1, x2, x3, ..., xn (nodes), using DIVIDED DIFFERENCE METHOD we represent the polynomial in the form

#                       P_n(x) = a0 + a1*(x-x_0) + a2*(x-x_0)(x-x_1) + ... + an*(x-x_0)*...*(x-x_(n-1)) 

# Where the ai's are the constants found bu computing the ith divided difference: f[x_(0),....,x_(0+i)]


#The argument A gets the input data matrix of the form: A = [[x1, x2, x3, ..., xn],[f(x1), f(x2), f(x3), ..., f(xn)]]


import sympy as sym
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

# Defining a function to compute the ith term of the polynomial --> a_i(x-x_0)....(x-x_(i-1)) 
# 
# where a_i is the ith divided difference f[x_(0),....,x_(0+i)]

def I_Term(A,i):
    # Writing the expression symbolically
    x = sym.Symbol('x')
    # The 1st term is calculated
    if i == 0:
        i_term = Dividediff(A,0,0)
    # The ith term is calculated, i != 1
    else:
        # We perform a loop to compute the product a_i(x-x_0)....(x-x_(i-1))
        j = 0 # Iterative variable 
        while j<=i:
            # The initial term of the product
            if j == 0:
                i_term = Dividediff(A,i,0)
            # Iteratively multipying successive terms to the product
            else:
                i_term = i_term*(x-A[0][j-1])
            # Incrementing the iterative variable to continue the loop
            j = j+1
    return i_term

# Defining a function to compute the polynomial P_n(x) (as shown above) symbolically and compute the value at the given point

# The argument 'val' tells the value at which the polynomial must be computed

def Dividepoly(A,val):
    x = sym.Symbol('x') # Used to compute the polynomial symbolically
    i = 0 # Iterative variable
    n = len(A[0]) # The total number of nodes
    Poly = 0 # Intializing the Polynomial with value 0
    # Compute the sum using a loop, which gives the desired polynomial
    while i<n:
        Poly = Poly + I_Term(A,i)
        i = i+1 # Incrementing the iterative variable to continue the loop
    print("-------------------------------------------------")
    print("Newton Divided Difference Polynomial")    
    print(Poly.subs(x,val)) 
    print(Poly)