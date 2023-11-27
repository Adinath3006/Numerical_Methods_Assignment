import cmath as cm

import numpy as np

def Muller(f,p0,p1,p2,N):
    print("-------------------------------------------------")
    print("Muller's Method")
    h_1 = p1-p0
    h_2 = p2-p1
    Delta1 = (f(p1)-f(p0))/h_1
    Delta2 = (f(p2)-f(p1))/h_2
    d = (Delta2 - Delta1)/(h_2+h_1) 
    i = 3
    while i <= N:
        b = Delta2 + (h_2*d)
        D = cm.sqrt(b**2 - 4*f(p2)*d)
        if np.sign(b)==np.sign(D):
            E = b+D
        else:
            E = b-D
        h = -2*((f(p2))/E)
        p = p2+h
        print("The", i ,"th iteration found the approximate root to be:", p)

        #Preparation for the next iteration

        p0 = p1
        p1 = p2
        p2 = p
        h_1 = p1-p0
        h_2 = p2-p1
        Delta1 = (f(p1)-f(p0))/h_1
        Delta2 = (f(p2)-f(p1))/h_2
        d = (Delta2 - Delta1)/(h_2+h_1)
        
        i = i+1

