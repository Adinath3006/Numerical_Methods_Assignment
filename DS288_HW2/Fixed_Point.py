def Fixed_Point(g,rTOL,p0,N):
    print("-------------------------------------------------")
    print("Fixed Point Iteration")
    i = 1
    while i <= N:
        p = g(p0)
        if abs(p-p0)/abs(p0) < rTOL:
            print("The approximate root with relative tolerance of ",rTOL," is found in ", i," iterations and its value is:",p)
            break
        else:
            print("The approximate root in the ",i, "th iteration is:",p)
            # Preparation for the next iteration
            p0 = p
            i = i+1
            if i == N:
                print("Maximum number of iterations reached, the desired in value with ",rTOL," relative tolerance was not found!!")
        