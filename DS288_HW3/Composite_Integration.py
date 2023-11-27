# Refer NUMERICAL ANALYSIS by Richard L. Burden; Chapter 4: Numerical Differentiation and Integration; 
# Section 4.4: Composite Numerical Integration and Section 4.5: Romberg Integration
# Notations used as in the above book

# 'f' is the function that is to be integrated in the limits 'lowlim' to 'upplim'

def Composite_Integration(f,method,lowlim,upplim,n):
    
    # Method = 0 ==> Composite Trapezoidal rule
    # Method = 1 ==> Composite Simpson's rule
    # Method = 2 ==> Romberg Integration
    
    if method == 0:
        h = (upplim - lowlim)/n
        
        # Computing the sum f(x_1) + ... + f(x_n-1) and x_i = lowlim + i*h
        Sum_f_x_i = 0
        i = 1
        while i<=n-1:
            x_i = lowlim + i*h
            Sum_f_x_i = Sum_f_x_i + f(x_i)
            i = i+1
            
        # Computing the integral
        Integral = (h/2)*(f(lowlim) + 2*Sum_f_x_i + f(upplim))
        
        print("")
        print("-------------------------------------------------")
        print("       COMPOSITE TRAPEZOIDAL RULE ( n =", n,")")
        print("-------------------------------------------------")
        print("")
        
        print("")
        print("The integral is evaluated to be:", Integral)
        return Integral
    
    elif method == 1:
        
        # n must be even for this method 
        if n%2 == 0:
            h = (upplim - lowlim)/n
        
            # Computing the sum f(x_2) + ... + f(x_n-2) and x_i = lowlim + i*h
            Sum_f_x_i = 0
            i = 1
            while i<= n/2 -1:
                x_i = lowlim + 2*i*h
                Sum_f_x_i = Sum_f_x_i + f(x_i)
                i = i+1
            
            # Computing the sum f(x_1) + ... + f(x_n-1) and x_i = lowlim + i*h
            Sum_f_x_j = 0
            j = 1
            while j<= n/2:
                x_j = lowlim + (2*j-1)*h
                Sum_f_x_j = Sum_f_x_j + f(x_j)
                j = j+1
                
            # Computing the integral
            Integral = (h/3)*(f(lowlim) + 2*Sum_f_x_i + 4*Sum_f_x_j + f(upplim))
            
            print("")
            print("-------------------------------------------------")
            print("        COMPOSITE SIMPSON'S RULE ( n =", n,")")
            print("-------------------------------------------------")
            print("")
        
            print("")
            print("The integral is evaluated to be:", Integral)
            return Integral
        
        else:
            print("")
            print("-------------------------------------------------")
            print("        COMPOSITE SIMPSON'S RULE ( n =", n,")")
            print("-------------------------------------------------")
            print("")
            
            print("")
            print("Invalid input! n must be an even integer")

    # The algorithm used for Romberg Integration is followed from the book

    elif method == 2:

        import numpy as np

        print("")
        print("-------------------------------------------------")
        print("              ROMBERG INTEGRATION")
        print("-------------------------------------------------")
        print("")
        
        h = upplim - lowlim
        i = 2
        R2j = np.zeros([n])
        R1j = np.zeros([n])
        R1j[0] = (h/2)*( f(upplim) + f(lowlim) )
        print(f'R_1,{1}:',R1j[0])
        while i<=n:
            Sum = 0
            for k in np.arange(1,2**(i-2)+1):
                Sum = Sum + f(lowlim + (k-0.5)*h)
            R2j[0] = (1/2)*(R1j[0] + h*Sum)
            j = 2
            print(f'R_2,{1}:',R2j[0])
            while j<=i:
                R2j[j-1] = R2j[j-2] + (R2j[j-2] - R1j[j-2])/(4**(j-1) - 1)
                print(f'R_2,{j}:',R2j[j-1])
                j = j+1
            h = h/2
            R1j = R2j
            i = i+1
    
    else:
        print("Invalid input! Please choose method 0 or 1 or 2")