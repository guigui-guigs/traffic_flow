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

def rho_0(x,type="creneau"): # créneau sur [0.1;0.2]
    if type == "creneau":
        if x < 0.3 and x> 0.1 : 
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

road_length = 0.5
simu_duration = 10 

#########################################################
###################### Schémas ##########################
#########################################################

def schemas_couplage(mailles, CI, type, N, dt, Vmax, pause):
    dx = road_length/(N+1)
    l = dx/2 # taille d'une voiture
    sommets = np.zeros(N+1)
    for i in range (0,N+1):
        sommets[i] = dx/2 + i*dx # coordonnées xi+1/2 des extrémités des mailles
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

    if t==0:
        old_choc = 0.3

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

        # Calcul de la vitesse du choc
        new_choc = centres[point_choc(centres, U)]
        #print(new_choc)
        print((new_choc-old_choc)/dt)
        old_choc = new_choc

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
        #print(U)
        plt.figure(1)
        plt.clf()
        plt.plot(centres,U)
        plt.plot(centres, Y, "og")
        plt.pause(pause)

# Construction d'un schéma qu'on appelle à chaque pas de temps
def schemas_couplage_iteratif(start_road, road_length, l, Vmax, mailles, type, sommets, centres, Uold, dt, insert, x1, x2):

    N = len(sommets) - 1 # signification de N ?
    # Retourne : densité U, les sommets, les centres et le pas de temps dt
    U = [0 for i in range(0,N)]
    if mailles == "fixes":
        sigma = [0 for i in range(0,N+1)] # convention i + 1/2
    if mailles == "FTL":
        sigma = [0 for i in range(0,N+1)]
        for i in range (1,N):
            sigma[i] = V(Vmax,l/(centres[i]-centres[i-1]))
        sigma[0] = 0
        sigma[N] = 0
        if len(sommets) > 2:
            dt = CFL(Uold, sommets,sigma,l) 

    G = [0 for i in range(0,N+1)]

    # Il faut voir ce qu'on fait, peut-être ne pas du tout réinsérer de mailles comme ce n'est plus périodique...
    for i in range (1,N):
        sommets[i] = sommets[i] + dt*sigma[i]
    to_insert = False
    if sommets[-1] > start_road + road_length :
        sommets = np.delete(sommets,-1)
        U = np.delete(U,-1)
        Uold = np.delete(Uold,-1)
        N = N-1
        if sommets[-1] > 1:
            to_insert = True

    if (insert == True):
        sommets = np.insert(sommets, 0, 0)

    centres = 0.5*(sommets[:-1] + sommets[1:])
    if type == "upwind":
        for i in range (1,N): # convention i + 1/2
            if (f_prime((Uold[i]+Uold[i-1])/2) - sigma[i]) >= 0:
                G[i] = f(Uold[i-1]) - sigma[i-1]*Uold[i-1]
            else: 
                G[i] = f(Uold[i]) - sigma[i]*Uold[i]
        G[0] = 0
    for i in range (0,N):
        U[i] = ((sommets[i+1]-sommets[i])*Uold[i] - dt*(G[i+1]-G[i]))/(sommets[i+1]-sommets[i] + (sigma[i+1]-sigma[i])*dt)
    return([U, sommets, centres, dt, to_insert])

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

def verif_entropic(x1,x2): 
    # x1 et x2 sont les coordonnées des chocs (dans le cas d'un créneau) où x2 > x1
    if (f(rho_0(x2,"creneau"))-f(rho_0(x1,"creneau")))/(rho_0(x2,"creneau")-rho_0(x1,"creneau")) < f_prime(rho_0(x2,"creneau")) and (f(rho_0(x2,"creneau"))-f(rho_0(x1,"creneau")))/(rho_0(x2,"creneau")-rho_0(x1,"creneau")) > f_prime(rho_0(x1,"creneau")) : 
        return "L'équation admet une solution entropique"
    else : 
        return "L'équation n'admet pas de solution entropique"

def point_choc(centres,U):
    # on parcourt la liste des U[i] pour voir où la dérivée est la plus grande
    # on suit la vitesse de ce point donné
    derive = []
    i = 0
    while i<len(centres)-1:
        derive.append((U[i+1]-U[i])/(centres[i+1]-centres[i]))
        i += 1 
    return np.argmax(derive)

##### INFORMATIONS #####

#schemas(mailles,CI,type,dx,dt)
# Pour mailles, "fixes", "alea" et "FTL" sont possibles
# Pour CI c'est-à-dire les conditions initiales, "creneau", "constante" et "parabole" sont possibles
# Pour type, c'est le type du schéma, "upwind" sont disponibles

##### EXEMPLES #####

#schemas_couplage("fixes","creneau","upwind", 50,0.004, 0.01, 1)
#schemas_couplage("alea","creneau","upwind", 50, 0.004, 0.01, 0.1) 
#schemas_couplage("FTL","creneau","upwind", 20, 0.004, 0.01, 0.1) 
#verif_entropic(0.2, 0.4)