# Refer NUMERICAL ANALYSIS by Richard L. Burden; Chapter 4: Numerical Differentiation and Integration; 
# Section 4.7: Gaussian Quadrature
# Notations used as in the above book

import sympy as sym
from Legendre_Poly import P

# The function P evaluates the Legendre Polynomial 
# n is the order of the Gaussian Quadrature
# lowlim = lower limit in the integration, a
# upplim = upper limit in the integration, b
# f is the function to be integrated

# Function to calculate the required coefficients needed in evaluating the integrals

def Coefficient(n,i,Roots):

    # Declaring the symbolic variable 'x'
    x = sym.symbols('x')

    # Calculating the coefficient, c_n,i for the corresponding roots r_n,i
    
    # Computing the Product term
    j = 1
    Product = 1
    while j <= n:
        if i != j:
            Product = Product*((x-Roots[j-1])/(Roots[i-1]-Roots[j-1]))
        else:
            Product = Product
        j = j+1
    
    # Integrate the Product to get the value of the coefficient
    from sympy import integrate
    Coeff = integrate(Product,(x,-1,1))
    
    return Coeff
    
# Function to evaluate the Integral using the Gaussian Quadrature Method    

def Gaussian_Quadrature(n,f,lowlim,upplim):

    # Initializing the variable which stores the evaluated integral 
    I = 0
    
    # Roots of the Legendre polynomial of nth order stored as an array
    Roots = sym.solve(P(n))
    
    # The Integral after transforming into an integral over [-1,1]
    i = 1
    while i <= n:
        I = I + ((upplim - lowlim)/2)*(Coefficient(n,i,Roots))*(f(((upplim - lowlim)*Roots[i-1] + (upplim + lowlim))/2))
        i = i+1
    
    print("")
    print("-------------------------------------------------")
    print("           GAUSSIAN QUADRATURE ( n =",n,")")
    print("-------------------------------------------------")
    print("")
    print(" The integral evaluated in the interval [",lowlim,",",upplim,"] in this method is:", sym.N(I,10))
    return sym.N(I,10)
    