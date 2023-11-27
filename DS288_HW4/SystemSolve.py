def system_solver(fsys,a,b,m,N,init):

    import numpy as np
    import sympy as sym

    # Initializing the parameters
    h = (b-a)/N
    t = a
    W = sym.symbols(f'W(0:{m+1})')

    # Creating a array which stores the values w_i
    array_w = np.zeros(m+1)
    k = np.zeros([4,m+1])
    k[0][0] = h
    k[1][0] = h
    k[2][0] = h
    k[3][0] = h
    soln = np.zeros([N+1,m])
    array_w[0] = t

    for j in range(1,m+1):
        array_w[j] = init[j-1]
        soln[0][j-1] = array_w[j]

    print(array_w)

    for i in range(1,N+1):
        # Compute k1_j
        for j in range(1,m+1):
            k[0][j] = h*fsys[j-1].subs({W[i]:array_w[i] for i in range(0,m+1)})
                
        # Compute k2_j
        for j in range(1,m+1):
            k[1][j] = h*fsys[j-1].subs({W[i]:array_w[i] + (1/2)*k[0][i] for i in range(0,m+1)})
        
        # Compute k3_j
        for j in range(1,m+1):
            k[2][j] = h*fsys[j-1].subs({W[i]:array_w[i] + (1/2)*k[1][i] for i in range(0,m+1)})
        
        # Compute k4_j
        for j in range(1,m+1):
            k[3][j] = h*fsys[j-1].subs({W[i]:array_w[i] + k[2][i] for i in range(0,m+1)})

        # Compute w_j
        for j in range(1,m+1):
            array_w[j] = array_w[j] + ( k[0][j] + 2*k[1][j]+ 2*k[2][j]+ k[3][j] )/6
            soln[i][j-1] = array_w[j]
        t = a + i*h
        array_w[0] = t
        #print(array_w)

    print("")
    print("-------------------------------------------------")
    print("                 Runge-Kutta Method              ")
    print("-------------------------------------------------")
    print("")
    print(soln)
    return soln