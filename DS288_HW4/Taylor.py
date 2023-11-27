import numpy as np
import sympy as smp
import matplotlib.pyplot as plt

def diff_f(f,p):
    x = smp.symbols('x')
    y = smp.symbols('y')
    f_p = f
    for i in np.arange(0,p):
        f_p = smp.diff(f_p,x) + smp.diff(f_p,y)*f

    return f_p

def Taylor_Method(a,b,h,f,s,n): 
    
    N = int((b-a)/h)
    X = np.zeros(N+1)
    Y = np.zeros(N+1)

    x = smp.symbols('x')
    y = smp.symbols('y')
    T_n = f


    for k in np.arange(1,n+1):
        T_n = T_n + (h**(k)/smp.factorial(k+1))*(diff_f(f,k))
    
    for i in np.arange(0,N+1):
        X[i] = a + i*h
    
    Y[0] = s

    for j in np.arange(0,N):
        Y[j+1] = Y[j] + h*T_n.subs({x:X[j],y:Y[j]})

    print("")
    print("-------------------------------------------------")
    print("             Taylor order 4  (h =",h,")          ")
    print("-------------------------------------------------")
    print("")
    return X,Y
