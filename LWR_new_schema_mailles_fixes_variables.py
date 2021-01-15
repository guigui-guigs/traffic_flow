#########################################################
#########################################################
################ LWR - mailles fixes ####################
#########################################################
#########################################################

import numpy as np
import matplotlib.pyplot as plt
import math
from random import uniform

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

def V(Vmax,x):
    return Vmax*max([0,10-x])

#########################################################
###################### Schémas ##########################
#########################################################

def schemas(mailles, CI, type, dx, dt):
    l = dx/10 # taille d'une voiture
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
    if mailles == "fixes":
        sigma = [0 for i in range(0,N+1)] # convention i + 1/2
    if mailles == "alea": 
        sigma = [0 for i in range(0,N+1)]
        for i in range (1,N):
            sigma[i] = uniform(0,road_length/50)/dt # on prend un nombre aléatoire entre 0 et un 50e de la longeur de la route, divisé par dt, pour que ce soit visible à l'échelle de temps où on regarde
        sigma[0] = 0
        sigma[N] = 0
        dt = CFL(sommets,sigma,l) # on teste de prendre dx/10 comme taille de voiture
    if mailles == "FTL":
        sigma = [0 for i in range(0,N+1)]
        for i in range (1,N):
            sigma[i] = V(road_length/50,l/(centres[i]-centres[i-1]))
        sigma[0] = 0
        sigma[N] = 0
        dt = CFL(sommets,sigma,l) # on teste de prendre dx/10 comme taille de voiture

    G = [0 for i in range(0,N+1)]
    for i in range (0,N):
        U[i] = rho_0(centres[i],CI)
    t=0
    plt.figure(1)
    plt.clf()
    plt.plot(centres,U)
    plt.pause(1)
    while t < T :
        if mailles == "alea":
            for i in range (1,N):
                sigma[i] = uniform(0,road_length/50)/dt
            sigma[0] = 0
            sigma[N] = 0
            dt = CFL(sommets,sigma,dx/10)
        if mailles == "FTL":
            for i in range (1,N):
                sigma[i] = V(road_length/50,l/(centres[i]-centres[i-1]))
            sigma[0] = 0
            sigma[N] = 0
            dt = CFL(sommets,sigma,l)
        for i in range (0,N+1):
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
        plt.pause(1)

def CFL(sommets,sigma,l):
    # La condition CFL qui doit être vérifiée est : 
    # dt(sigma[i]-sigma[i+1]) < sommets[i+1] - sommets[i] - l où l est la taille d'un véhicule. On choisira de prendre dt = contrainte/2
    # cette contrainte n'a du sens que si sigma[i] > sigma[i+1], autrement il n'y a pas de contrainte sur dt causée par les sigma, et on prent dt = dx/2 

    # Cette fonction renvoie le dt à utiliser pour chaque nouvelle itération

    dt_list = []
    for i in range (0,len(sommets)-1):
        if sigma[i] > sigma[i+1]:
            dt_list.append((sommets[i+1]-sommets[i]-l)/(2*(sigma[i]-sigma[i+1]))) # on pourrait aussi ne pas diviser par 2, à tester
        else:
            dt_list.append((sommets[i+1]-sommets[i])/2)
    return min(dt_list)

##### INFORMATIONS #####

#schemas(mailles,CI,type,dx,dt)
# Pour mailles, "fixes", "alea" et "FTL" sont possibles
# Pour CI c'est-à-dire les conditions initiales, "creneau", "constante" et "parabole" sont possibles
# Pour type, c'est le type du schéma, "upwind" sont disponibles

##### EXEMPLES #####

#schemas("fixes","parabole","upwind",0.01,0.004)
#schemas("alea","parabole","upwind",0.01, 0.004) # le dx et le dt ne sont que pour les conditions initiales ici
#schemas("FTL","parabole","upwind",0.01, 0.004) # le dx et le dt ne sont que pour les conditions initiales ici