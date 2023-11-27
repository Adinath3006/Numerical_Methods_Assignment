import math
import numpy as np

# Derivative computing function at x numerically
def Diff(f,x): 
    h = 1e-5
    return (f(x+h)-f(x-h))/(2*h)

#Double derivative computing function at x numerically
def dDiff(f,x): 
    h = 1e-5
    return (f(x+2*h)+f(x-2*h)-2*f(x))/(4*h**2)

# Defining the function which performs the Newton's method
def Newton(x0,rTOL,N,f):
    print("")
    print("                Newton's Method")
    i = 1 # Iterating variable
    while i <= N: # The iterating variable is less than the limit, N, if it is true the loop procceds
        fprime = Diff(f,x0) # Derivative of the function at x0
        x1 = x0-(f(x0)/fprime) # The x intercept of the line with slope "fprime" at 'x0'
        if f(x1)==0 or math.isclose(x1,x0,rel_tol = rTOL): # x can either be the proper root or satisfy the relative error condition
            print("The approximate solution obtained in",i,"iteration is",x1)
            print("The root of the equation with relative tolerance of",rTOL,"is",x1)
            print("The number of iterations taken to achieve this value is",i)
            break
        else: # If the "if" condition is not satisfied then it should proceed to try the x-intercept formed by the line at f(x1) and slope f'(x1)
            print("The approximate solution obtained in",i,"iteration is",x1)
            x0 = x1 # Change the value of x0 to x1
        i = i+1 # Incrementing the iterating variable
        if i > N: # When the maximum number of iterations is reached and we still don't get the solution, we display a message
            print("")
            print("Maximum number of iterations reached!")
            print("No root is found upto",N,"iterations is found!")
            break
        else:
            continue

# Defining the function which performs the Modified Newton's method
def ModifNewton(x0,rTOL,N,f):
    print("")
    print("            Modified Newton's Method")
    i = 1 # Iterating variable
    while i <= N: # The iterating variable is less than the limit, N, if it is true the loop procceds
        fprime = Diff(f,x0) # Derivative of the function at x0
        fdprime = dDiff(f,x0) # Double derivative of the function at x0
        x1 = x0-(f(x0)*fprime/(fprime**2-(f(x0)*fdprime))) # The x intercept of the line with slope "fprime/(fprime**2-(f(x0)*fdprime))" at 'x0'
        if f(x1)==0 or math.isclose(x1,x0,rel_tol = rTOL): # x can either be the proper root or satisfy the relative error condition
            print("The approximate solution obtained in",i,"iteration is",x1)
            print("The root of the equation with relative tolerance of",rTOL,"is",x1)
            print("The number of iterations taken to achieve this value is",i)
            break
        else: # If the "if" condition is not satisfied then it should proceed to try the x-intercept formed by the line at f(x1) and the iterative slope
            print("The approximate solution obtained in",i,"iteration is",x1)
            x0 = x1 # Change the value of x0 to x1
        i = i+1 # Incrementing the iterating variable
        if i > N: # When the maximum number of iterations is reached and we still don't get the solution, we display a message
            print("")
            print("Maximum number of iterations reached!")
            print("No root is found upto",N,"iterations is found!")
            break
        else:
            continue