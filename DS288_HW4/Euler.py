# Refer NUMERICAL ANALYSIS by Richard L. Burden; Chapter 5: Initial-Value Problems for Ordinary Differential Equations; 
# Section 5.1: Euler's Method 
# Notations used as in the above book

def Euler_Method(f,a,b,y_0,h,disp):

    # disp = 0 ---> No title displayed, else title is displayed

    import numpy as np
    from scipy.integrate import odeint

    N = int((b-a)/h)

    # Creating a array which stores the values w_i
    array_w = np.zeros(N+1)

    t = a

    # Computing w_i's recursively 

    i = 0
    while i <= N:
        if i == 0:
            array_w[0] = y_0
        else:
            array_w[i] = array_w[i-1] + h*f(t,array_w[i-1])
            t = a + i*h
        i = i+1
    
    # Output to be displayed

    if disp == 0:
        print(array_w)
        return array_w
    else:
        print("")
        print("-------------------------------------------------")
        print("             Euler's Method (h =",h,")           ")
        print("-------------------------------------------------")
        print("")
        print(array_w)
        return array_w




# ---------------------------------------------------------------------------------------------------------------------------------------

# Refer NUMERICAL ANALYSIS by Richard L. Burden; Chapter 5: Initial-Value Problems for Ordinary Differential Equations; 
# Section 5.4: Runge-Kutta Methods ---> Modified Euler Method 
# Notations used as in the above book

def Modif_Euler_Method(f,a,b,y_0,h,disp):

    # disp = 0 ---> No title displayed, else title is displayed

    import numpy as np
    from scipy.integrate import odeint

    N = int((b-a)/h)

    # Creating a array which stores the values w_i
    array_w = np.zeros(N+1)

    t = a

    # Computing w_i's recursively 

    i = 0
    while i <= N:
        if i == 0:
            array_w[0] = y_0
        else:
            array_w[i] = array_w[i-1] + (h/2)*( f(t,array_w[i-1]) +f(t+h,array_w[i-1]+h*f(t,array_w[i-1])) )
            t = a + i*h
        i = i+1
    
    # Output to be displayed

    if disp == 0:
        return array_w
    else:
        print("")
        print("-------------------------------------------------")
        print("       Modified Euler's Method  (h =",h,")       ")
        print("-------------------------------------------------")
        print("")
        print(array_w)
        return array_w
