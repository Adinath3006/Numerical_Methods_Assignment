# Refer NUMERICAL ANALYSIS by Richard L. Burden; Chapter 4: Numerical Differentiation and Integration; 
# Section 4.1: Numerical Differentiation 
# Notations used as in the above book

# The argument A gets the input data matrix of the form: A = [[x1, x2, x3, ..., xn],[f(x1), f(x2), f(x3), ..., f(xn)]]

# ----------------------------------------------------------------------------------------------------------------------------------------

import sympy as sym

# The Five Point Formula 
def Five_Point_Method(A,method,val,end):

    # end = 0 ==> Left endpoint 
    # end = 1 ==> Right endpoint

    # Method = 0 ==> Midpoint Formula
    # Method = 1 ==> Endpoint Formula

    # Midpoint formula
    if method == 0:

        print("")
        print("-------------------------------------------------")
        print("           FIVE POINT MIDPOINT FORMULA           ")
        print("-------------------------------------------------")
        print("")

        if len(A[0]) < 5:
            print("Invalid input!")
        else:
            
            # Check if the input value is valid and this method can be used to evaluate the Derivative at that point

            # Calculates the index of the input variable
            j = 0
            while j < len(A[0]):
                if val == A[0][j]:
                    break
                else:
                    j = j+1

            # Checks if the index lies in the acceptable range
            
            # If the value is not in the input array then the following message is displayed
            if j > len(A[0])-1:
                print("Invalid input!")
                print("Input value out of the range of the input data")

            # If the value lies in the region near the endpoints such that the midpoint formula cannot be used then the following message is displayed
            elif j < 2 or len(A[0])-j <3:
                print("Invalid input!")
                print("Evaluation of the derivative using this method for the given value is not possible!")
                
            # If yes then it proceeds to evaluate the derivative
            else:
                
                # Iteratively calculate derivative using this method for increasing values of h (if possible) 
                i = 1
                while i > 0:
                    if j+2*i < len(A[0]) and j-2*i>=0:
                        h = val - A[0][j-i]
                        Deriv_f = (1/(12*h))*(A[1][j-2*i] - 8*A[1][j-i] + 8*A[1][j+i] - A[1][j+2*i])
                        print("The derivative evaluated at",val,"using this method with h =",(sym.N(h,7)),"is:",sym.N(Deriv_f,8))
                        i = i+1
                    else:
                        break
    
    # Endpoint formula
    elif method == 1:

        print("")
        print("-------------------------------------------------")
        print("           FIVE POINT ENDPOINT FORMULA           ")
        print("-------------------------------------------------")
        print("")

        if len(A[0]) < 5:
            print("Invalid input!")
        else:
            
            # Check if the input value is valid and this method can be used to evaluate the Derivative at that point

            # Calculates the index of the input variable
            j = 0
            while j < len(A[0]):
                if val == A[0][j]:
                    break
                else:
                    j = j+1

            # Checks if the index lies in the acceptable range
            
            # Left Endpoint is selected
            if end == 0:
                print("----------------- LEFT ENDPOINT -----------------")
                print("")
                # If the value is not in the input array then the following message is displayed
                if j > len(A[0])-1:
                    print("Invalid input!")
                    print("Input value out of the range of the input data")
            
                # If the value, with index j, and the point with index j+4 is not in the input data set then the following message is displayed
                elif j + 4 > len(A[0])-1:
                    print("Invalid input!")
                    print("Evaluation of the derivative using this method for the given value is not possible!")

                # If yes then it proceeds to evaluate the derivative
                else:
                    
                    # Iteratively calculate derivative using this method for increasing values of h (if possible) 
                    i = 1
                    while i > 0:
                        if j+4*i < len(A[0]):
                            h = A[0][j+i] - val
                            Deriv_f = (1/(12*h))*(-25*A[1][j] + 48*A[1][j+1] - 36*A[1][j+2] + 16*A[1][j+3] - 3*A[1][j+4])
                            print("The derivative evaluated at",val,"using this method with h =",(sym.N(h,7)),"is:",sym.N(Deriv_f,8))
                        else:
                            break
            
            # Right Endpoint is selected
            elif end == 1:
                print("---------------- RIGHT ENDPOINT -----------------")
                print("")
                # If the value is not in the input array then the following message is displayed
                if j > len(A[0])-1:
                    print("Invalid input!")
                    print("Input value out of the range of the input data")
            
                # If the value, with index j, and the point with index j-4 is not in the input data set then the following message is displayed
                elif j - 4 < 0:
                    print("Invalid input!")
                    print("Evaluation of the derivative using this method for the given value is not possible!")

                # If yes then it proceeds to evaluate the derivative
                else:
                    
                    # Iteratively calculate derivative using this method for increasing values of h (if possible) 
                    i = 1
                    while i > 0:
                        if j-4*i >= 0:
                            h = A[0][j-i] - val
                            Deriv_f = (1/(12*h))*(-25*A[1][j] + 48*A[1][j-1] - 36*A[1][j-2] + 16*A[1][j-3] - 3*A[1][j-4])
                            print("The derivative evaluated at",val,"using this method with h =",(sym.N(h,7)),"is:",sym.N(Deriv_f,8))
                        else:
                            break
            
            else:
                print("Invalid endpoint selected!")
    
    else:
        print("Invalid method selected!")

# ----------------------------------------------------------------------------------------------------------------------------------------

def Three_Point_Method(A,method,val,end):
    # end = 0 ==> Left endpoint 
    # end = 1 ==> Right endpoint

    # Method = 0 ==> Midpoint Formula
    # Method = 1 ==> Endpoint Formula

    # Midpoint formula
    if method == 0:

        print("")
        print("-------------------------------------------------")
        print("          THREE POINT MIDPOINT FORMULA           ")
        print("-------------------------------------------------")
        print("")

        if len(A[0]) < 3:
            print("Invalid input!")
        else:
            
            # Check if the input value is valid and this method can be used to evaluate the Derivative at that point

            # Calculates the index of the input variable
            j = 0
            while j < len(A[0]):
                if val == A[0][j]:
                    break
                else:
                    j = j+1

            # Checks if the index lies in the acceptable range
            
            # If the value is not in the input array then the following message is displayed
            if j > len(A[0])-1:
                print("Invalid input!")
                print("Input value out of the range of the input data")

            # If the value lies in the region near the endpoints such that the midpoint formula cannot be used then the following message is displayed
            elif j < 1 or len(A[0])-j <2:
                print("Invalid input!")
                print("Evaluation of the derivative using this method for the given value is not possible!")
                
            # If yes then it proceeds to evaluate the derivative
            else:

                # Iteratively calculate derivative using this method for increasing values of h (if possible) 
                i = 1
                while i > 0:
                    if j+i < len(A[0]) and j-i>=0:
                        h = val - A[0][j-i]
                        Deriv_f = (1/(2*h))*(A[1][j+i] - A[1][j-i])
                        print("The derivative evaluated at",val,"using this method with h =",(sym.N(h,7)),"is:",sym.N(Deriv_f,8))
                        i = i+1
                    else:
                        break
    
    # Endpoint formula
    elif method == 1:

        print("")
        print("-------------------------------------------------")
        print("          THREE POINT ENDPOINT FORMULA           ")
        print("-------------------------------------------------")
        print("")

        if len(A[0]) < 3:
            print("Invalid input!")
        else:
            
            # Check if the input value is valid and this method can be used to evaluate the Derivative at that point

            # Calculates the index of the input variable
            j = 0
            while j < len(A[0]):
                if val == A[0][j]:
                    break
                else:
                    j = j+1

            # Checks if the index lies in the acceptable range
            
            # Left Endpoint is selected
            if end == 0:
                print("----------------- LEFT ENDPOINT -----------------")
                print("")
                # If the value is not in the input array then the following message is displayed
                if j > len(A[0])-1:
                    print("Invalid input!")
                    print("Input value out of the range of the input data")
            
                # If the value, with index j, and the point with index j+4 is not in the input data set then the following message is displayed
                elif j + 2 > len(A[0])-1:
                    print("Invalid input!")
                    print("Evaluation of the derivative using this method for the given value is not possible!")

                # If yes then it proceeds to evaluate the derivative
                else:
                    
                    # Iteratively calculate derivative using this method for increasing values of h (if possible) 
                    i = 1
                    while i > 0:
                        if j+2*i < len(A[0]):
                            h = A[0][j+i] - val
                            Deriv_f = (1/(2*h))*(-3*A[1][j] + 4*A[1][j+i] - A[1][j+2*i])
                            print("The derivative evaluated at",val,"using this method with h =",(sym.N(h,7)),"is:",sym.N(Deriv_f,8))
                            i = i+1
                        else:
                            break
            
            # Right Endpoint is selected
            elif end == 1:
                print("---------------- RIGHT ENDPOINT -----------------")
                print("")
                # If the value is not in the input array then the following message is displayed
                if j > len(A[0])-1:
                    print("Invalid input!")
                    print("Input value out of the range of the input data")
            
                # If the value, with index j, and the point with index j-4 is not in the input data set then the following message is displayed
                elif j - 2 < 0:
                    print("Invalid input!")
                    print("Evaluation of the derivative using this method for the given value is not possible!")

                # If yes then it proceeds to evaluate the derivative
                else:
                    
                    # Iteratively calculate derivative using this method for increasing values of h (if possible) 
                    i = 1
                    while i > 0:
                        if j-2*i >= 0:
                            h = A[0][j-i] - val
                            Deriv_f = (1/(2*h))*(-3*A[1][j] + 4*A[1][j-i] - A[1][j-2*i])
                            print("The derivative evaluated at",val,"using this method with h =",(sym.N(h,7)),"is:",sym.N(Deriv_f,8))
                            i = i+1
                        else:
                            break
                    
            else:
                print("Invalid endpoint selected!")
    
    else:
        print("Invalid method selected!")

# ----------------------------------------------------------------------------------------------------------------------------------------

def Second_Derivative_Midpoint_Formula(A,val):
    
    import numpy as np

    print("")
    print("-------------------------------------------------")
    print("        SECOND DERIVATIVE MIDPOINT FORMULA       ")
    print("-------------------------------------------------")
    print("")
    
    if len(A[0])<3:
        print("Invalid input!")
    
    else:
        # Check if the input value is valid and this method can be used to evaluate the Derivative at that point

        # Calculates the index of the input variable
        j = 0
        while j < len(A[0]):
            if val == A[0][j]:
                break
            else:
                j = j+1

        # Checks if the index lies in the acceptable range
            
        # If the value is not in the input array then the following message is displayed
        if j > len(A[0])-1:
            print("Invalid input!")
            print("Input value out of the range of the input data")

        # If the value lies in the region near the endpoints such that the midpoint formula cannot be used then the following message is displayed
        elif j < 1 or len(A[0])-j <2:
            print("Invalid input!")
            print("Evaluation of the derivative using this method for the given value is not possible!")
                
        # If yes then it proceeds to evaluate the derivative
        else:
                
            # Iteratively calculate derivative using this method for increasing values of h (if possible)
            
            # Calculating the maximum possible no.of times we can calculate the derivative iteratively
            i = 1
            while i > 0:
                if j+i < len(A[0]) and j-i>=0:
                    i = i+1
                else:
                    break
            
            # Maximum possible times the derviative can be computed for allowed values of h
            max_i = i -1
            
            k = 1
            Derivative = np.zeros(max_i)
            while k<=max_i:
                h = val - A[0][j-k]
                Deriv_f = (1/(h**2))*(A[1][j-k] - 2*A[1][j] + A[1][j+k])
                print("The derivative evaluated at",val,"using this method with h =",(sym.N(h,7)),"is:",sym.N(Deriv_f,8))
                Derivative[k-1] = Deriv_f
                k = k+1
            
            return Derivative