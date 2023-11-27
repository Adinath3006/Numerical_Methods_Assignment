

# Problem 1

print("-------------------------------------------------")
print("                    PROBLEM 1                    ")
print("-------------------------------------------------")

def f1(x,y):
    return (x-y)/2
a1 = 0
b1 = 3
y1_0 = 1

from scipy.integrate import odeint
import sympy as smp
import numpy as np
import matplotlib.pyplot as plt

x =smp.symbols('x')
g =smp.Function('g')
f_temp = smp.dsolve(smp.Derivative(g(x),x) + (g(x)-x)/2,ics = {g(0):1})
f_exact = smp.lambdify(x,f_temp.rhs)
x_d = np.linspace(0,3,1000)
y_d = f_exact(x_d)

print(f"The exact solution is {f_temp.rhs}")

plt.plot(x_d,y_d, c = 'm', label = "Exact solution")

# Part(a)
from Euler import *
dat11 = Euler_Method(f1,a1,b1,y1_0,1,1)
x11 = np.linspace(a1,b1,int((b1-a1)/1) + 1)
plt.plot(x11,dat11,label = 'h=1')

dat12 = Euler_Method(f1,a1,b1,y1_0,0.5,1)
x12 = np.linspace(a1,b1,int((b1-a1)/0.5) + 1)
plt.plot(x12,dat12,label = 'h=0.5')

dat13 = Euler_Method(f1,a1,b1,y1_0,0.25,1)
x13 = np.linspace(a1,b1,int((b1-a1)/0.25) + 1)
plt.plot(x13,dat13,label = 'h=0.25')

dat14 = Euler_Method(f1,a1,b1,y1_0,0.125,1)
x14 = np.linspace(a1,b1,int((b1-a1)/0.125) + 1)
plt.plot(x14,dat14,label = 'h=0.125')

plt.title('Comparison between Euler Method for different h')
plt.legend()

plt.show()

# Part (b)
plt.plot(x_d,y_d, c = 'm', label = "Exact solution")

dat21 = Modif_Euler_Method(f1,a1,b1,y1_0,1,1)
x11 = np.linspace(a1,b1,int((b1-a1)/1) + 1)
plt.plot(x11,dat21,label = 'h=1')

dat22 = Modif_Euler_Method(f1,a1,b1,y1_0,0.5,1)
x12 = np.linspace(a1,b1,int((b1-a1)/0.5) + 1)
plt.plot(x12,dat22,label = 'h=0.5')

dat23 = Modif_Euler_Method(f1,a1,b1,y1_0,0.25,1)
x13 = np.linspace(a1,b1,int((b1-a1)/0.25) + 1)
plt.plot(x13,dat23,label = 'h=0.25')

dat24 = Modif_Euler_Method(f1,a1,b1,y1_0,0.125,1)
x14 = np.linspace(a1,b1,int((b1-a1)/0.125) + 1)
plt.plot(x14,dat24,label = 'h=0.125')

plt.title('Comparison between Modified Euler Method for different h')
plt.legend()

plt.show()

# Part(c)
from Taylor import *

x,y =smp.symbols('x y')
f = (x-y)/2
X_1, Y_1 = Taylor_Method(0,3,1,f,1,4)
X_2, Y_2 = Taylor_Method(0,3,0.5,f,1,4)
X_3, Y_3 = Taylor_Method(0,3,0.25,f,1,4)
X_4, Y_4 = Taylor_Method(0,3,0.125,f,1,4)


plt.plot(x_d,y_d, label = "Exact solution")
plt.plot(X_1,Y_1, label = " h = 1")
plt.plot(X_2,Y_2, label = " h = 0.5")
plt.plot(X_3,Y_3, label = " h = 0.25")
plt.plot(X_4,Y_4, label = " h = 0.125")
plt.title("Comparison between Taylor's Method for different h")
plt.legend()
plt.show()

# Part (d)
plt.plot(x_d,y_d, c = 'm', label = "Exact solution")

from RK4 import *

dat41 = Runge_Kutta(f1,a1,b1,y1_0,1,1)
plt.plot(x11,dat41,label = 'h =1')

dat42 = Runge_Kutta(f1,a1,b1,y1_0,0.5,1)
plt.plot(x12,dat42,label = 'h =0.5')

dat43 = Runge_Kutta(f1,a1,b1,y1_0,0.25,1)
plt.plot(x13,dat43,label = 'h =0.25')

dat44 = Runge_Kutta(f1,a1,b1,y1_0,0.125,1)
plt.plot(x14,dat44,label = 'h =0.125')

plt.title('Comparison between Runge-Kutta for various h')
plt.legend()

plt.show()

# Part (e) and (f)
plt.plot(x_d,y_d, c = 'm', label = "Exact solution")

from Multistep import *
dat51 = MultistepMethod(f1,a1,b1,y1_0,0.125,0,1)
plt.scatter(x14,dat51,label = 'Adam-Bashforth h=0.125')

plt.title('Adam Bashforth')
plt.legend()

plt.show()

plt.plot(x_d,y_d, c = 'm', label = "Exact solution")

dat61 = MultistepMethod(f1,a1,b1,y1_0,0.125,1,1)
plt.scatter(x14,dat61,label = 'Adam-Moulton h=0.125')

plt.title('Comparison between Adam-Bashforth and Adam-Moulton')
plt.legend()

plt.show()



# Problem 2
print("-------------------------------------------------")
print("                    PROBLEM 2                    ")
print("-------------------------------------------------")
W = smp.symbols('W(0:3)')
fsys = [W[2], -(4*W[2]+5*W[1])]
init = [3,-5]

from SystemSolve import *

dat2 = system_solver(fsys,0,5,2,50,init)

def ty(x):
    return 3*np.exp(-2*x)*np.cos(x) + np.exp(-2*x)*np.sin(x)

x_t = np.linspace(0,3,1000)
y_t = ty(x_d)

plt.plot(x_t,y_t, c = 'm', label = "Exact solution")

xpoint = np.linspace(0,5,50)
ypoint = np.zeros(50)

for i in range(0,50):
    ypoint[i] = dat2[i][0]

plt.scatter(xpoint,ypoint,label = 'Runge-Kutta Solution')
plt.title('Comparison between Actual solution and Runge-Kutta method')
plt.legend()
plt.show()

# Problem 3
print("-------------------------------------------------")
print("                    PROBLEM 3                    ")
print("-------------------------------------------------")

#Part (a)
W = smp.symbols('W(0:3)')
fsys2 = [W[2], ((2*W[0])/(1+(W[0])**2))*W[2] - ((2)/(1+(W[0])**2))*W[1] + 1]

init1 = [1.25,0]
Soln1 = system_solver(fsys2,0,4,2,20,init1)
init2 = [0,1]
Soln2 = system_solver(fsys2,0,4,2,20,init2)

Xpt = np.linspace(0,4,20)
ysoln1 = np.zeros(20)
ysoln2 = np.zeros(20)

for i in range(0,20):
    ysoln1[i] = Soln1[i][0]
    ysoln2[i] = Soln2[i][0]

plt.plot(Xpt,ysoln1 + ((-0.95-ysoln1[19])/(ysoln2[19]))*(ysoln2), label = 'Shooting Method')

# Part (b)
from Finitediff import *

x = smp.symbols('x')

coeff = [ (2*x)/(1+(x**2)), (-2)/(1+(x**2)), x/x ]

soln1 = MatrixW(coeff,0,4,20,[1.25,-0.95])

x3pt = np.linspace(0,4,20)
y3pt = np.zeros(20)

for i in range(0,20):
    y3pt[i] = soln1[i][0]

plt.scatter(x3pt,y3pt,label='Finite-Difference')

plt.title('Comparison between Shooting Method and Finite-Difference method')


#Part(c)
def tty(x):
    return 1.25 + 0.486089652*x - 2.25*x**2 + 2*x*np.arctan(x) - (np.log(1 + x**2))/2 + (np.log(1+x**2))*(x**2)/2

x_tt = np.linspace(0,4,1000)
y_tt = tty(x_d)


plt.plot(x_tt,y_tt, c = 'm', label = "Exact solution")

plt.legend()
plt.show()