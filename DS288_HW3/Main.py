import numpy as np

import sympy as sym

# Problem 1
print("-------------------------------------------------")
print("                    PROBLEM 1                    ")
print("-------------------------------------------------")
# Input Data
Data1 = [[0, 0.4, 0.8, 1.2, 1.6], [2.90, 3.10, 3.56, 4.60, 6.70]]

# Output Data
from Least_Square_Poly_Fit import Least_Square_Fit_Poly

exp1 = Least_Square_Fit_Poly(Data1,2)
exp2 = Least_Square_Fit_Poly(Data1,3)

from Curve_Fit import Fit_Symbolic_Exp

import matplotlib.pyplot as plt

plt.figure(1)
Fit_Symbolic_Exp(Data1,exp1,0,0,1000,1)
Fit_Symbolic_Exp(Data1,exp2,0,0,1000,1)
plt.title('Least Square Fit')
plt.xlabel('x')
plt.ylabel('y')
plt.show()

print("-----------------------------------------------------------------------------")

# Problem 2
print("-------------------------------------------------")
print("                    PROBLEM 2                    ")
print("-------------------------------------------------")

# Input Data
Data2 = [[0.1, 0.2, 0.3, 0.4], [0.76, 0.58, 0.44, 0.35]]

# Output Data
from Least_Square_Fit import Least_Square_Fit

a = sym.symbols('a')
b = sym.symbols('b')
x =sym.symbols('x')

param = [a,b]
func = [sym.exp(-3*x), sym.exp(-2*x)]

Least_Square_Fit(Data2,param,func)

print("-----------------------------------------------------------------------------")

# Problem 3
print("-------------------------------------------------")
print("                    PROBLEM 3                    ")
print("-------------------------------------------------")
# Input Data
Data = [[1.8, 1.9, 2.0, 2.1, 2.2], [10.889365, 12.703199, 14.778112, 17.148957, 19.855030]]

# Output 
from Differentiation import Five_Point_Method, Three_Point_Method, Second_Derivative_Midpoint_Formula
Three_Point_Method(Data,0,2.0,0)
Three_Point_Method(Data,1,2.0,0)
Three_Point_Method(Data,1,2.0,1)
Five_Point_Method(Data,0,2.0,0)
Five_Point_Method(Data,1,2.0,0)
Five_Point_Method(Data,1,2.0,1)
val1 = Second_Derivative_Midpoint_Formula(Data,2.0)
val2 = Second_Derivative_Midpoint_Formula(Data,1.9)
val3 = Second_Derivative_Midpoint_Formula(Data,2.1)
Data_deriv = [[1.9,2.0,2.1],[val2[0],val1[0],val3[0]]]
Three_Point_Method(Data_deriv,0,2.0,0)
print("-----------------------------------------------------------------------------")


# Problem 4
print("-------------------------------------------------")
print("                    PROBLEM 4                    ")
print("-------------------------------------------------")
# Input Data
# Defining the function 
lowlim = 0
upplim = np.pi
def f(x):
    return (np.sin(x))**2

# Output
from Composite_Integration import Composite_Integration

Composite_Integration(f,0,lowlim,upplim,2)
Composite_Integration(f,0,lowlim,upplim,4)
Composite_Integration(f,0,lowlim,upplim,8)
Composite_Integration(f,0,lowlim,upplim,16)

Composite_Integration(f,1,lowlim,upplim,2)
Composite_Integration(f,1,lowlim,upplim,4)
Composite_Integration(f,1,lowlim,upplim,8)
Composite_Integration(f,1,lowlim,upplim,16)

Composite_Integration(f,2,lowlim,upplim,5)

print("-----------------------------------------------------------------------------")

# Problem 5
print("-------------------------------------------------")
print("                    PROBLEM 5                    ")
print("-------------------------------------------------")
# Input Data

# Defining the fucntion which needs to be integrated
import numpy as np
def f(x):
    return (x**2)*((np.e)**x)

from Gaussian_Quad import Gaussian_Quadrature
# Part (a)
print("------------------- Part (a) --------------------")

# Output 
val1 = Gaussian_Quadrature(8,f,-1,1)
print("")

# Part (b)
print("------------------- Part (b) --------------------")

# Output 
val21 = Gaussian_Quadrature(4,f,-1,0)
val22 = Gaussian_Quadrature(4,f,0,1)
print("")
print(" The integral evaluated in the interval [",-1,",",1,"] is:", val21+val22)
print("")

# Part (c)
print("------------------- Part (c) --------------------")

# Output 
val31 = Gaussian_Quadrature(2,f,-1,-0.5)
val32 = Gaussian_Quadrature(2,f,-0.5,0)
val33 = Gaussian_Quadrature(2,f,0,0.5)
val34 = Gaussian_Quadrature(2,f,0.5,1)
print("")
print(" The integral evaluated in the interval [",-1,",",1,"] is:", val31+val32+val33+val34)