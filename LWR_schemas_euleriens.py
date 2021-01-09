#########################################################
#########################################################
############# LWR - schémas eulériens ###################
#########################################################
#########################################################
import numpy as np
import matplotlib.pyplot as plt
import math

def rho_0(x):
    return 1/(1+np.exp(-5*x))

def f(rho):
    return rho*(1-rho)

def schemas(type, dx, dt):
    #dx = 0.1
    #dt = 0.1
    road_length = 1
    simu_duration = 10 
    N = math.floor(road_length/dx)
    X = np.linspace(0,road_length,N+1)
    T = simu_duration/dt
    U = [0 for i in range(0,N+1)]
    for i in range (0,N+1):
        U[i] = rho_0(X[i])
    if type == "LF":
        t = 0
        while t < T :
            # il faut faire attention aux termes de bord, pour l'instant manque le premier et le dernier, conditions périodiques ?
            Uold = U
            print(Uold)
            for i in range(1,N):
                F_amont = (f(Uold[i])+f(Uold[i+1]))/2 - (dx/(2*dt))*(Uold[i+1]-Uold[i])
                F_aval = (f(Uold[i-1])+f(Uold[i]))/2 - (dx/(2*dt))*(Uold[i]-Uold[i-1])
                U[i] = Uold[i] + (dt/dx)*(F_amont-F_aval)
            t = t + dt
            #plt.plot(X,U)
            #plt.show()
        print(U)
        plt.plot(X,U)
        plt.show()

schemas("LF",0.1,1)