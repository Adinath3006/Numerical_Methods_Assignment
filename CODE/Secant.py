import math

# Defining the function which performs the secant method
def Secant(a,b,rTOL,N,f):
    print("")
    print("                 Secant Method")
    x_curr = b # The variable to store the value of x_(n)
    x_prev = a # The variable to store the value of x_(n-1)
    i = 1 # iterating variable
    while i <= N: # The iterating variable is less than the limit, N, if it is true the loop procceds
        m = (f(x_curr)-f(x_prev))/(x_curr-x_prev) # Slope of the line formed by f(x_(n-1)) and f(x_(n))
        x_next = x_prev - (f(x_prev)/m) # The variable to store the value of x_(n+1) which is the intercept of the line formed by f(x_(n-1)) and f(x_(n))
        if f(x_next)==0 or abs(x_next-x_curr)/abs(x_next)<=rTOL: # x can either be the proper root or satisfy the relative error condition
            print("The approximate solution obtained in",i,"iteration is",x_curr)
            print("The root of the equation with relative tolerance of",rTOL,"is",x_next)
            print("The number of iterations taken to achieve this value is",i)
            break
        else: # If the "if" condition is not satisfied then it should proceed to try the x-intercept formed by f(x_(n)) and f(x_(n+1))
            print("The approximate solution obtained in",i,"iteration is",x_curr)
            x_prev = x_curr # Change the value of x_(n-1) to x_n
            x_curr = x_next # Change the value of x_(n) to x_(n+1)
        i = i+1 # Incrementing the iterating variable
        if i > N: # When the maximum number of iterations is reached and we still don't get the solution, we display a message
            print("")
            print("Maximum number of iterations reached!")
            print("No root is found upto",N,"iterations is found!")
            break
        else:
            continue