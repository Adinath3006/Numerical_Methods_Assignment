import numpy as np
def h(X,i):
    return X[0][i+1]-X[0][i]

def MatrixA(X):
    n = len(X[0])-1
    i = 0
    A = np.zeros((n+1, n+1), dtype=int)
    while i<n+1:
        j = 0
        while j<n+1:
            if abs(i-j) > 1:
                A[i][j] = 0
            else:
                if i > j:
                    if i == n:
                        A[i][j] = 0
                    else:
                        A[i][j] = h(X,i-1)
                elif j > i:
                    if i == 0:
                        A[i][j] = 0
                    else:
                        A[i][j] = h(X,j-1)
                else:
                    if i == 0 or i == n:
                        A[i][j] = 1
                    else:
                        A[i][j] = 2*(h(X,i-1)+h(X,i))
            j = j+1
        i =i+1
    return A

def Matrixb(X):
    n = len(X[0])-1
    i = 0
    b = np.zeros((n+1, 1))
    while i<n+1:
        if i == 0 or i == n:
            b[i][0] = 0
        else:
            b[i][0] = ((3*(X[1][i+1]-X[1][i]))/(h(X,i)))-((3*(X[1][i]-X[1][i-1]))/(h(X,i-1)))
        i = i+1
    return b

def Matrixc(X):
    return np.matmul(np.linalg.inv(MatrixA(X)),Matrixb(X))

def coefficent_b(X,i):
    return ((X[1][i+1]-X[1][i])/(h(X,i)))-(h(X,i)*((2*Matrixc(X)[i][0])+Matrixc(X)[i+1][0])/(3))

def coefficent_d(X,i):
    return (Matrixc(X)[i+1][0]-Matrixc(X)[i][0])/(3*h(X,i))

def Spline(X,k):
    import sympy as sym
    x = sym.Symbol('x')
    i = k
    S = X[1][i]+(coefficent_b(X,i)*(x-X[0][i]))+(Matrixc(X)[i][0]*(x-X[0][i])**2)+(coefficent_d(X,i)*(x-X[0][i])**3)
    return (sym.N(S,4))

def PolySpline(X,val):
    print("-------------------------------------------------")
    print("Cubic Spline")
    import sympy as sym
    x = sym.Symbol('x')
    i = 0
    n = len(X[0]) # Total length of the array A[0]
    while i<n-1:
        print("S(x) = ",Spline(X,i)," for [",X[0][i],",",X[0][i+1],"]")
        i = i+1
    i = 0
    while i<n-1:
        if val >= X[0][i] and val < X[0][i+1]:
            print(Spline(X,i).subs(x,val))
            break
        else:
            i = i+1