#########################################################
#########################################################
################ LWR - mailles fixes ####################
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

def schemas(CI, type, dx, dt):
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
    #print(centres)
    T = simu_duration/dt
    U = [0 for i in range(0,N)]
    sigma = [0 for i in range(0,N+1)] # convention i + 1/2
    G = [0 for i in range(0,N+1)]
    for i in range (0,N):
        U[i] = rho_0(centres[i],CI)
    t=0
    plt.figure(1)
    plt.clf()
    plt.plot(centres,U)
    plt.pause(2)
    while t < T :
        for i in range (0,N):
            sommets[i] = sommets[i] + dt*sigma[i]
        centres = 0.5*(sommets[:-1] + sommets[1:])
        Uold = U
        if type == "upwind":
            for i in range (0,N-1): # convention i + 1/2
                if f_prime((Uold[i+1]+Uold[i])/2) >= 0:
                    G[i+1] = f(Uold[i]) - sigma[i+1]*Uold[i]
                else: 
                    G[i+1] = f(Uold[i+1]) - sigma[i+1]*Uold[i]
            if f_prime((Uold[0]+Uold[N-1])/2) >= 0:
                G[N] = f(Uold[N-1]) - sigma[N]*Uold[N-1]
            else: 
                G[N] = f(Uold[0]) - sigma[N]*Uold[N-1]
        for i in range (0,N):
            U[i] = ((sommets[i+1]-sommets[i])*Uold[i] - dt*(G[i+1]-G[i]))/(sommets[i+1]-sommets[i] + (sigma[i+1]-sigma[i])*dt)
        t = t + dt
        plt.figure(1)
        plt.clf()
        plt.plot(centres,U)
        plt.pause(2)

    

schemas("parabole","upwind",0.01,0.004)