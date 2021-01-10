#########################################################
#########################################################
############# LWR - schémas eulériens ###################
#########################################################
#########################################################

import numpy as np
import matplotlib.pyplot as plt
import math

#########################################################
################## Paramétrisation ######################
#########################################################

def rho_0(x): # créneau sur [0;0.2]
    if x < 0.2 : 
        return 1
    else :
        return 0

def f(rho):
    return rho*(1-rho)

def f_prime(rho):
    return 1 - 2*rho

#########################################################
###################### Schémas ##########################
#########################################################

def schemas(type, dx, dt):
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
        plt.figure(1)
        plt.clf()
        plt.plot(X,U)
        while t < T :
            Uold = U
            for i in range(1,N):
                F_amont = (f(Uold[i])+f(Uold[i+1]))/2 - (dx/(2*dt))*(Uold[i+1]-Uold[i])
                F_aval = (f(Uold[i-1])+f(Uold[i]))/2 - (dx/(2*dt))*(Uold[i]-Uold[i-1])
                U[i] = Uold[i] + (dt/dx)*(F_amont-F_aval)
            U[0] = Uold[-1] # domaine périodique
            U[-1] = Uold[0] # domaine périodique
            t = t + dt
            plt.figure(1)
            plt.clf()
            plt.plot(X,U)
            plt.pause(0.1)
    if type == "LW":
        t = 0
        plt.figure(1)
        plt.clf()
        plt.plot(X,U)
        while t < T :
            Uold = U
            for i in range(1,N):
                F_amont = (f(Uold[i])+f(Uold[i+1]))/2 - (dx/(2*dt))*f_prime((Uold[i+1]+Uold[i])/2)*(f(Uold[i+1]-f(Uold[i])))
                F_aval = (f(Uold[i-1])+f(Uold[i]))/2 - (dx/(2*dt))*f_prime((Uold[i]+Uold[i-1])/2)*(f(Uold[i]-f(Uold[i-1])))
                U[i] = Uold[i] + (dt/dx)*(F_amont-F_aval)
            U[0] = Uold[-1] # domaine périodique
            U[-1] = Uold[0] # domaine périodique
            t = t + dt
            plt.figure(1)
            plt.clf()
            plt.plot(X,U)
            plt.pause(0.1)

schemas("LF",0.1,0.2)
#schemas("LW",0.1,0.2)