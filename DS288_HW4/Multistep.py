def MultistepMethod(f,a,b,y_0,h,method,disp):

    # disp = 0 ---> No title displayed, else title is displayed

    # Method = 0 --->  fourth-order Adams-Bashforth technique
    # Method = 1 --->  fourth-order Adams-Moulton technique

    from RK4 import Runge_Kutta

    if method == 0:
        
        import numpy as np

        # Initializing the parameters
        N = int((b-a)/h)
        t_i = a
        
        # Creating a array which stores the values w_i
        array_w = np.zeros(N+1)

        i = 0
        while i < 4:
            array_w[i] = Runge_Kutta(f,a,b,y_0,h,0)[i]
            i = i+1
        
        j = 3
        while j<N:
            array_w[j+1] = array_w[j] + (h/24)*( 55*f(a+j*h,array_w[j]) - 59*f(a+(j-1)*h,array_w[j-1]) + 37*f(a+(j-2)*h,array_w[j-2]) - 9*f(a+(j-3)*h,array_w[j-3]) )
            j = j+1
        
        if disp == 0:
            return array_w
        else:
            print("")
            print("-------------------------------------------------")
            print("         4-step Adam-Bashforth  (h =",h,")       ")
            print("-------------------------------------------------")
            print("")
            print(array_w)
            return array_w
        
    elif method == 1:

        import numpy as np
        import sympy as sym

        # Initializing the parameters
        N = int((b-a)/h)
        t_i = a
        
        # Creating a array which stores the values w_i
        array_w = np.zeros(N+1)

        i = 0
        while i < 3:
            array_w[i] = Runge_Kutta(f,a,b,y_0,h,0)[i]
            i = i+1
        
        j = 2
        while j<N:  

            w = sym.symbols('w')

            sum =  19*f(a+j*h,array_w[j]) - 5*f(a+(j-1)*h,array_w[j-1]) + f(a+(j-2)*h,array_w[j-2])
            eqn = w - array_w[j] - (h/24)*( 9*f(a+(j+1)*h,w) + sum )
            soln = sym.solve(eqn)
            array_w[j+1] = soln[0]
            j = j+1

        if disp == 0:
            return array_w
        else:
            print("")
            print("-------------------------------------------------")
            print("         3-step Adam-Moulton  (h =",h,")       ")
            print("-------------------------------------------------")
            print("")
            print(array_w)
            return array_w
    
    else:
        print("Invalid input! Please select the correct method - 0/1")

