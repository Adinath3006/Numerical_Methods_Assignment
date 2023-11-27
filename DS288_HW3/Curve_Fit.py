def Fit_Symbolic_Exp(A,f,x_lowlim,x_upplim,n_x,method):
    
    # Method = 0 ==> Plot the function in the given range [x_lowlim, x_upplim]
    # Method = 1 ==> Plot the function and A[1](i) in the range A[0] 
    
    import numpy as np
    import sympy as sym
    import matplotlib.pyplot as plt
    
    x = sym.symbols('x')
    
    # [x_lowlim, x_upplim] ---> The x range which needs to be plotted
    # n_x ---> No. of points to be evaluated
    
    if method == 0:
        X = np.linspace(x_lowlim, x_upplim, n_x)
        f_X = np.zeros([n_x,1])
    
        # Creating an array which stores all the values f(x_i)
        i = 0
        while i < n_x:
            f_X[i] = f.subs(x,X[i])
            i = i+1
        
        plt.plot(X,f_X)
        plt.show()
    
    elif method == 1:
        X = np.linspace(A[0][0], A[0][len(A[0])-1], n_x)
        f_X = np.zeros([n_x,1])
    
        # Creating an array which stores all the values f(x_i)
        i = 0
        while i < n_x:
            f_X[i] = f.subs(x,X[i])
            i = i+1
        
        plt.plot(X,f_X)
        plt.plot(A[0],A[1],'o')
        
    else:
        print("Invalid input! Please choose method 0 or 1")