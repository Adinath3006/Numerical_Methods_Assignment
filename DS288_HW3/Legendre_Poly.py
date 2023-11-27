# Calculating the legendre polynomial "Bonnet’s recurrence relation"

#                                   n*P_n(x) = (2*n-1)*x*P_n-1(x) - (n-1)*P_n-2(x)                        

def P(n):
    import sympy as sym
    x = sym.symbols('x')
    if(n == 0):
        return 1 # P0 = 1
    elif(n == 1):
        return x # P1 = x
    else:
        return sym.simplify((((2 * n)-1)*x * P(n-1)-(n-1)*P(n-2))/float(n))