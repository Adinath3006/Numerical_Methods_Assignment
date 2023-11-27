import math

def Bisection(a,b,TOL,N,f):
    print("")
    print("            Bisection Method")
    # Check if the initial conditions satisfy the condition 
    if f(a)*f(b)>0:
        print("The initial conditions do not contain the root between!")
    else:
        i = 1 # Iterating variable
        while i <= N: # The iterating variable is less than the limit, N, if it is true the loop procceds
            x_curr = (a+b)/2 # The mid point is chosen to be checked if it is a root
            if f(x_curr) == 0 or abs(b-a)/abs(b)<TOL:  # x can either be the proper root or satisfy the error condition
                print("The approximate solution obtained in",i,"iteration is",x_curr)
                print("The root of the equation with relative tolerance of",TOL,"is",x_curr)
                print("The number of iterations taken to achieve this value is",i)
                break
            else: # If not then update the current values accordingly for the next run
                print("The approximate solution obtained in",i,"iteration is",x_curr)
                if f(a)*f(x_curr)<0:
                    b = x_curr
                else:
                    a = x_curr               
            i = i+1
            if i > N: # When the maximum number of iterations is reached and we still don't get the solution, we display a message
                    print("")
                    print("Maximum number of iterations reached!")
                    print("No value satisfying the condition is found!")
                    break
            else:
                continue
