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

def rho_0(x,type): # créneau sur [0;0.2]
    if type == "creneau":
        if x < 0.2 : 
            return 1
        else :
            return 0
    if type == "constante":
        return 1
    if type == "parabole":
        if x<1:
            return 10*x*(1-x)
        else:
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
    pas = dx*np.ones(N)
    sommets = np.zeros(N+1)
    sommets[1:] = np.cumsum(pas) # coordonnées xi+1/2 des extrémités des mailles
    #print(sommets)
    mes = [] # les mesures des mailles mi
    for i in range(N) : 
        mes.append(sommets[i+1]-sommets[i])
    #print(mes)
    centres = 0.5*(sommets[:-1] + sommets[1:]) # les coordonnées xi des centres des mailles
    print(centres)
    T = simu_duration/dt
    U = [0 for i in range(0,N)]
    for i in range (0,N):
        U[i] = rho_0(centres[i],"parabole")
    if type == "LF":
        t = 0
        plt.figure(1)
        plt.clf()
        plt.plot(centres,U)
        plt.pause(2)
        while t < T :
            Uold = U
            for i in range(1,N-1):
                F_amont = (f(Uold[i])+f(Uold[i+1]))/2 - (dx/(2*dt))*(Uold[i+1]-Uold[i])
                F_aval = (f(Uold[i-1])+f(Uold[i]))/2 - (dx/(2*dt))*(Uold[i]-Uold[i-1])
                U[i] = Uold[i] + (dt/dx)*(F_amont-F_aval)
            # Calcul du premier terme 
            F_amont = (f(Uold[0])+f(Uold[1]))/2 - (dx/(2*dt))*(Uold[1]-Uold[0])
            F_aval = (f(Uold[-1])+f(Uold[0]))/2 - (dx/(2*dt))*(Uold[0]-Uold[-1])
            U[0] = Uold[0] + (dt/dx)*(F_amont-F_aval)
            # Calcul du dernier terme 
            F_amont = (f(Uold[-1])+f(Uold[0]))/2 - (dx/(2*dt))*(Uold[0]-Uold[-1])
            F_aval = (f(Uold[-2])+f(Uold[-1]))/2 - (dx/(2*dt))*(Uold[-1]-Uold[-2])
            U[-1] = Uold[-1] + (dt/dx)*(F_amont-F_aval)

            t = t + dt
            plt.figure(1)
            plt.clf()
            plt.plot(centres,U)
            plt.pause(2)
    if type == "LW":
        t = 0
        plt.figure(1)
        plt.clf()
        plt.plot(centres,U)
        plt.pause(2)
        while t < T :
            Uold = U
            for i in range(1,N-1):
                F_amont = (f(Uold[i])+f(Uold[i+1]))/2 - (dx/(2*dt))*f_prime((Uold[i+1]+Uold[i])/2)*(f(Uold[i+1]-f(Uold[i])))
                F_aval = (f(Uold[i-1])+f(Uold[i]))/2 - (dx/(2*dt))*f_prime((Uold[i]+Uold[i-1])/2)*(f(Uold[i]-f(Uold[i-1])))
                U[i] = Uold[i] + (dt/dx)*(F_amont-F_aval)
            # Calcul du premier terme 
            F_amont = (f(Uold[0])+f(Uold[1]))/2 - (dx/(2*dt))*f_prime((Uold[1]+Uold[0])/2)*(f(Uold[1]-f(Uold[0])))
            F_aval = (f(Uold[-1])+f(Uold[0]))/2 - (dx/(2*dt))*f_prime((Uold[0]+Uold[-1])/2)*(f(Uold[0]-f(Uold[-1])))
            U[0] = Uold[0] + (dt/dx)*(F_amont-F_aval)
            # Calcul du dernier terme 
            F_amont = (f(Uold[-1])+f(Uold[0]))/2 - (dx/(2*dt))*f_prime((Uold[0]+Uold[-1])/2)*(f(Uold[0]-f(Uold[-1])))
            F_aval = (f(Uold[-2])+f(Uold[-1]))/2 - (dx/(2*dt))*f_prime((Uold[-1]+Uold[-2])/2)*(f(Uold[-1]-f(Uold[-2])))
            U[-1] = Uold[-1] + (dt/dx)*(F_amont-F_aval)

            t = t + dt
            plt.figure(1)
            plt.clf()
            plt.plot(centres,U)
            plt.pause(2)

schemas("LF",0.1,0.04)
#schemas("LW",0.1,0.05)