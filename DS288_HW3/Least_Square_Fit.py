# param ---> Parameters of the fit function, this should be an array of the form [a1, a2, ... , an] where a1,a2,.. are all symbolic variables
# func ---> Functions of the fit function, this should be an array of the form [f1(x), f2(x), ... , fn(x)]

def Least_Square_Fit(A,param,func):
    
    import numpy as np
    import sympy as sym
    from sympy.solvers import solve
    
    x = sym.symbols('x')
    
    # Creating the fit function from the arguments, exp = a1*f1(x) + a2*f2(x) + ... + an*fn(x)
    i = 0
    exp = 0
    while i<len(param):
        exp = exp + param[i]*func[i]
        i = i+1
    
    # Creating the error function
    j = 0
    E = 0
    while j<len(A[0]):
        E = E + ( exp.subs(x,A[0][j]) - A[1][j] )**2
        j = j+1
    
    # Creating an array of system of equations which on solving gives the desired parameters
    i = 0
    system = [0] * len(param)
    while i<len(param):
        system[i] = sym.diff(E,param[i])
        i = i+1
    
    Parameter = sym.solve(system,dict = True)
    
    # The final fit function
    Fin_exp = exp.subs(Parameter[0])

    print("")
    print("-------------------------------------------------")
    print("               LEAST SQUARE FIT ")
    print("-------------------------------------------------")
    print("")

    print("")
    print("The final fit function is:", Fin_exp)
    return Fin_exp