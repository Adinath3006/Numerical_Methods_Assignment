import numpy as np
from Bisection import Bisection
from Newton import ModifNewton, Newton
from RegulaFalsi import Regula

from Secant import Secant

# Defining the function
def f(x):
    return x**3-9*x+1

def f1(x):
    return np.cos(x)-x*np.exp(x)

def f2(x):
    return x +np.exp(-(x**2))*np.cos(x)

def f3(x):
    return (x +np.exp(-(x**2))*np.cos(x))**2

N = 100 # Maximum iterating variable
#PROBLEM 1

# Declaring initial values for the variables 
a = 2
b = 4 
TOL = 10**(-3) # Tolerance value
print("")
print("PROBLEM 1")
Bisection(a,b,TOL,N,f)
#PROBLEM 2

print("")
print("PROBLEM 2")
# Defining the initial values 
x0 = 0
x1 = 1
rTOL = 10**(-6)
Regula(x0,x1,rTOL,N,f1)
Secant(x0,x1,rTOL,N,f1)

#PROBLEM 3

print("")
print("PROBLEM 3")
print("")
print("Part (a)")
Newton(0,10**(-6),100,f2)
Secant(x0,x1,10**(-6),N,f2)
ModifNewton(0,10**(-6),100,f2)

print("")
print("Part (b)")
Newton(0,10**(-6),100,f3)
Secant(x0,x1,10**(-6),N,f3)
ModifNewton(0,10**(-6),100,f3)