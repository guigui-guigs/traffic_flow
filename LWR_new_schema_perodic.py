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

def rho_0(x,type="creneau"): # créneau sur [0;0.2]
    if type == "creneau":
        if x < 0.6 and x> 0.3 : 
            return 0.5
        else :
            return 0
    if type == "constante":
        return 0.5
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

road_length = 1
simu_duration = 10 

#########################################################
###################### Schémas ##########################
#########################################################

def schemas_test(mailles, CI, type, N, dt, Vmax, pause):
    dx = road_length/(N+1)
    l = dx/2 # taille d'une voiture
    sommets = np.zeros(N+1)
    for i in range (0,N+1):
        sommets[i] = dx/2 + i*dx # coordonnées xi+1/2 des extrémités des mailles
    print(sommets)
    centres = 0.5*(sommets[:-1] + sommets[1:]) # les coordonnées xi des centres des mailles
    #centres.append((sommets[0] + 1 - sommets[N])/2)
    Y = [0 for i in range(N)]
    T = simu_duration/dt
    U = [0 for i in range(0,N)]
    if mailles == "fixes":
        sigma = [0 for i in range(0,N+1)] # convention i + 1/2
    if mailles == "alea": 
        sigma = [0 for i in range(0,N+1)]
        for i in range (1,N):
            sigma[i] = uniform(0,V(Vmax,U[i]))
        sigma[0] = 0
        sigma[N] = 0
        dt = CFL(U, sommets,sigma,l) 
    if mailles == "FTL":
        sigma = [0 for i in range(0,N+1)]
        for i in range (1,N):
            sigma[i] = V(Vmax,l/(centres[i]-centres[i-1]))
        sigma[0] = 0
        sigma[N] = 0
        dt = CFL(U, sommets,sigma,l) 

    G = [0 for i in range(0,N+1)]
    for i in range (0,N):
        U[i] = rho_0(centres[i],CI)
    t=0
    plt.figure(1)
    plt.clf()
    plt.plot(centres,U)
    plt.plot(centres, Y, "og")
    plt.pause(pause)
    while t < T :
        if mailles == "alea":
            for i in range (1,N):
                sigma[i] = uniform(0,V(Vmax,U[i]))
            sigma[0] = 0
            sigma[N] = 0
            dt = CFL(U, sommets,sigma,l)
            #print(dt)
        if mailles == "FTL":
            for i in range (1,N):
                sigma[i] = V(Vmax,l/(centres[i]-centres[i-1]))         
            if sommets[0] + sommets[N] > 1:
                centre_N = centres[0] + centres[N-1] - 1 # on crée artificiellement un centre en x0 et xN
                sigma[0] = V(Vmax,l/(centres[0] - centre_N))
                #print(sigma[0])
                sigma[N] = V(Vmax,l/(1 + centre_N - centres[N-1])) 
                #print(sigma[N])
            else : 
                centre_N = centres[0] + centres[N-1]
                sigma[0] = V(Vmax,l/(centres[0] - (1 - centre_N)))
                #print(sigma[0])
                sigma[N] = V(Vmax,l/(centre_N - centres[N-1])) 
                #print(sigma[N])
            dt = CFL(U, sommets,sigma,l)
            #print(dt)

        for i in range (0,N+1):
            sommets[i] = sommets[i] + dt*sigma[i]

        if sommets[N] > 1 : 
            to_insert = sommets[N]
            sommets = np.delete(sommets,N)
            sommets = np.insert(sommets, 0, 0)

        centres = 0.5*(sommets[:-1] + sommets[1:])
        Uold = U
        if type == "upwind":
            for i in range (1,N): # convention i + 1/2
                if (f_prime((Uold[i]+Uold[i-1])/2) - sigma[i]) >= 0:
                    G[i] = f(Uold[i-1]) - sigma[i-1]*Uold[i-1]
                else: 
                    G[i] = f(Uold[i]) - sigma[i]*Uold[i]
            if (f_prime((Uold[0]+Uold[N-1])/2) - sigma[0]) >= 0:
                G[0] = f(Uold[N-1]) - sigma[N-1]*Uold[N-1]
            else: 
                G[0] = f(Uold[0]) - sigma[0]*Uold[0]
            G[N] = G[0] # conditions périodiques
        for i in range (0,N):
            U[i] = ((sommets[i+1]-sommets[i])*Uold[i] - dt*(G[i+1]-G[i]))/(sommets[i+1]-sommets[i] + (sigma[i+1]-sigma[i])*dt)
        t = t + dt
        plt.figure(1)
        plt.clf()
        plt.plot(centres,U)
        plt.plot(centres, Y, "og")
        plt.pause(pause)

'''
def schemas_couplage(mailles, CI, type, l, N, sommets, U, t, go):

    if go == False: # attention ne prend plus en compte les dt et les sigma comme ça...
        if mailles == "alea":
            for i in range (1,N):
                sigma[i] = uniform(0,V(Vmax,U[i]))
            sigma[0] = 0
            sigma[N] = 0
            dt = CFL(sommets,sigma,l)
        if mailles == "FTL":
            for i in range (1,N):
                sigma[i] = V(Vmax,l/(centres[i]-centres[i-1]))
            sigma[0] = 0
            sigma[N] = 0
            dt = CFL(sommets,sigma,l)
        return dt

    else:
        G = [0 for i in range(0,N+1)]
        if t == 0 : 
            for i in range (0,N):
                U[i] = rho_0(centres[i],CI)

        for i in range (0,N+1):
            sommets[i] = sommets[i] + dt*sigma[i]
        if sommets[-1] - sommets[-2] < l:
            sommets = np.delete(sommets,-2)
            sommets = np.insert(sommets, 1, (sommets[0]+sommets[1])/4)

        centres = 0.5*(sommets[:-1] + sommets[1:])
        Uold = U
        if type == "upwind":
            for i in range (1,N): # convention i + 1/2
                if (f_prime((Uold[i]+Uold[i-1])/2) - sigma[i]) >= 0:
                    G[i] = f(Uold[i-1]) - sigma[i-1]*Uold[i-1]
                else: 
                    G[i] = f(Uold[i]) - sigma[i]*Uold[i]
            if (f_prime((Uold[0]+Uold[N-1])/2) - sigma[0]) >= 0:
                G[0] = f(Uold[N-1]) - sigma[N-1]*Uold[N-1]
            else: 
                G[0] = f(Uold[0]) - sigma[0]*Uold[0]
            G[N] = G[0] # conditions périodiques
        for i in range (0,N):
            U[i] = ((sommets[i+1]-sommets[i])*Uold[i] - dt*(G[i+1]-G[i]))/(sommets[i+1]-sommets[i] + (sigma[i+1]-sigma[i])*dt)
        return [sommets, centres, U]
'''

def CFL(U, sommets, sigma, l):
    # La condition CFL qui doit être vérifiée est : 
    # dt(sigma[i]-sigma[i+1]) < sommets[i+1] - sommets[i] - l où l est la taille d'un véhicule. On choisira de prendre dt = contrainte/2
    # cette contrainte n'a du sens que si sigma[i] > sigma[i+1], autrement il n'y a pas de contrainte sur dt causée par les sigma, et on prent dt = dx/2 

    # Cette fonction renvoie le dt à utiliser pour chaque nouvelle itération

    dt_list = []
    for i in range (0,len(sommets)-2):
        if sigma[i] > sigma[i+1] :
            dt_list.append((sommets[i+1]-sommets[i]-l)/(2*(sigma[i]-sigma[i+1]))) # on pourrait aussi ne pas diviser par 2, à tester
        else:
            dt_list.append((sommets[i+1]-sommets[i])/2)
        #dt_list.append((sommets[i+1]-sommets[i])*f_prime(U[i])) # ça fait des dt trop petit... vraiment ça la condition ?
    return min(dt_list)

def verif_entropic(x1,x2): # je ne comprends pas comment écrire la formule
    # x1 et x2 sont les coordonnées des chocs (dans le cas d'un créneau) où x2 > x1
    if (f(rho_0(x2,"creneau"))-f(rho_0(x1,"creneau")))/(rho_0(x2,"creneau")-rho_0(x1,"creneau")) < f_prime(rho_0(x2,"creneau")) and (f(rho_0(x2,"creneau"))-f(rho_0(x1,"creneau")))/(rho_0(x2,"creneau")-rho_0(x1,"creneau")) > f_prime(rho_0(x1,"creneau")) : 
        return "L'équation admet une solution entropique"
    else : 
        return "L'équation n'admet pas de solution entropique"
'''
def vitesse(centres,U):
    while i<len(centres):
'''


##### INFORMATIONS #####

#schemas(mailles,CI,type,dx,dt)
# Pour mailles, "fixes", "alea" et "FTL" sont possibles
# Pour CI c'est-à-dire les conditions initiales, "creneau", "constante" et "parabole" sont possibles
# Pour type, c'est le type du schéma, "upwind" sont disponibles

##### EXEMPLES #####

#schemas_test("fixes","creneau","upwind", 50,0.004, 0.01, 0.01)
#schemas_test("alea","creneau","upwind", 50, 0.004, 0.01, 0.1) 
schemas_test("FTL","creneau","upwind", 50, 0.004, 0.01, 0.1) 
#verif_entropic(0.3, 0.6)