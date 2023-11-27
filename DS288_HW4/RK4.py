# Refer NUMERICAL ANALYSIS by Richard L. Burden; Chapter 5: Initial-Value Problems for Ordinary Differential Equations; 
# Section 5.1: Runge-Kutta Methods 
# Notations used as in the above book

def Runge_Kutta(f,a,b,y_0,h,disp):
    import numpy as np

    # disp = 0 ---> No title displayed, else title is displayed

    N = int((b-a)/h)

    # Creating a array which stores the values w_i
    array_w = np.zeros(N+1)

    t = a

    i = 0
    while i <= N:
        if i == 0:
            array_w[0] = y_0
        else:
            k_1 = h*f(t,array_w[i-1])
            k_2 = h*f(t+h/2,array_w[i-1]+k_1/2)
            k_3 = h*f(t+h/2,array_w[i-1]+k_2/2)
            k_4 = h*f(t+h,array_w[i-1]+k_3)
            array_w[i] = array_w[i-1] + (1/6)*(k_1 + 2*k_2 + 2*k_3 + k_4)
            t = a + i*h
        i = i+1
    
    # Output to be displayed

    if disp == 0:
        return array_w
    else:
        print("")
        print("-------------------------------------------------")
        print("          Runge-Kutta Order 4  (h =",h,")        ")
        print("-------------------------------------------------")
        print("")
        print(array_w)
        return array_w
    