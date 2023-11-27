from cmath import exp
import numpy as np
import math as mt
from Muller import Muller

# Problem 1
print("")
print("Problem 1")
print("")
# Defining the function  
def f(x):
    return np.cos(x)-(x*exp(x))
# Computing the value
Muller(f,-1.0,0.0,1.0,7)

#Problem 2
print("")
print("Problem 2")
print("")
from Fixed_Point import Fixed_Point
# Defining the function
def g(x):
    return x+(exp(-x))-np.sin(x)
# Computing the value
Fixed_Point(g,10**(-6),0.5,50)

# Problem 3
print("")
print("Problem 3")
print("")
# Importing the function
from Lagrange import Lagrange
# Input data
A = [[300, 304, 305, 307],
    [2.4771, 2.4829, 2.4843, 2.4871]]
val3 = 301
# Computing the value
Lagrange(A,val3)    

# Problem 4
print("")
print("Problem 4")
print("")
# Part (a)
# Importing the function
from Divided_difference import Dividepoly
# Input data 
Ba = [[0.5,  1.5, 3.0, 5.0, 6.5, 8.0],
    [1.625, 5.875, 31, 131.0, 282.125, 521.0]]
val4a = 7
# Computing the value
Dividepoly(Ba,val4a)

# Part (b)
# Importing the function
from Forward_difference import ForwardDiff_Poly
# Input data
B_b = [[0.1, 0.2, 0.3, 0.4, 0.5],
    [1.40, 1.56, 1.76, 2.00, 2.28]]
val4b1 = 0.25
val4b2 = 0.35
# Computing the value and the polynomial
ForwardDiff_Poly(B_b,(val4b1-0.1)/0.1)
ForwardDiff_Poly(B_b,(val4b2-0.1)/0.1)

C = [[1.0, 1.3, 1.6, 1.9, 2.2],[0.7651977, 0.6200860, 0.4554022, 0.2818186, 0.1103623]]
#Polynomial(C,1.1)

#Problem 5
print("")
print("Problem 5")
print("")
from Cubic_Spline import PolySpline
# Input data
D = [[3.0, 4.5, 7.0, 9.0], [2.5, 1.0, 2.5, 0.5]]
# Computing the function
PolySpline(D,7.5)
PolySpline(D,4)

#Problem 6
print("")
print("Problem 6")
print("")
from Hermite import Hermite_Divided_Diff, Hermite_Lagrange
# Input data
E = [[1.3, 1.6, 1.9], [0.6200860, 0.4554022, 0.2818186], [-0.5220232, -0.5698959, -0.5811571]]
# Computing the function
Hermite_Lagrange(E,1.5)
Hermite_Divided_Diff(E,1.5)

