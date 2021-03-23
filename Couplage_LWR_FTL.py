# Brainstorming

# A chaque pas de temps on fait tourner le schéma puis le coupleur 

from LWR_new_schema_perodic import schemas_couplage_iteratif
from LWR_new_schema_perodic import rho_0
from Coupleur_LWR_FTL import coupleur_LWRversFTL
from Coupleur_LWR_FTL import ajout_voitures_FTL
from Coupleur_LWR_FTL import coupleur_FTLversLWR

import numpy as np
import matplotlib.pyplot as plt
import math

def main():

    #########################################################
    ################### Paramétrisation #####################
    #########################################################

    dx = 0.01
    dt = 0.004 
    # Ces données évolueront mais sont ici pour paramétrer 
    simu_duration = 5
    T = simu_duration/dt
    taille_voiture = dx/2

    Vmax_LWR1 = 0.01
    Vmax_FTL = 0.01
    Vmax_LWR2 = 0.01

    M1 = 100 # nombre de voitures dans la zone LWR 1

    # Définition des différentes zones de circulation entre [0;1]
    x1 = 0.4 # début de la zone FTL
    x2 = 0.8 # fin de la zone FTL 

    N1 = math.floor(x1/dx) 
    N3 = math.floor((1-x2)/dx) 

    dx_1 = x1/(N1+1)
    dx_3 = (1-x2)/(N3+1)

    sommets_1 = np.zeros(N1+1)
    for i in range (0,N1+1):
        sommets_1[i] = dx_1/2 + i*dx_1 # coordonnées xi+1/2 des extrémités des mailles
    centres_1 = 0.5*(sommets_1[:-1] + sommets_1[1:]) # les coordonnées xi des centres des mailles

    sommets_2 = []
    centres_2 = []

    sommets_3 = np.zeros(N3+1)
    for i in range (0,N3+1):
        sommets_3[i] = x2 + dx_3/2 + i*dx_3 # coordonnées xi+1/2 des extrémités des mailles
    centres_3 = 0.5*(sommets_3[:-1] + sommets_3[1:]) # les coordonnées xi des centres des mailles

    U_1 = [0 for i in range(0,N1)]
    for i in range (0,N1):
        U_1[i] = rho_0(centres_1[i],"creneau") # voitures initialement présentes dans la zone LWR1
    U_2 = []
    U_3 = [0 for i in range(0,N3)]

    surface_LWR1 = calcul_aire(rho_0,0,1,1/1000)
    check_surface = surface_LWR1/M1 # surface correspondant à une voiture
    surface_tampon_1 = 0

    insert_LWR1 = False # on insère de la densité enf fonction du flux G mais lequel ?
    insert_LWR2 = False

    t = 0
    centres = [*centres_1, *centres_2, *centres_3]
    U = [*U_1, *U_2, *U_3]
    plt.figure(1)
    plt.clf()
    plt.plot(centres,U)
    Y = [0 for i in range(len(centres))]
    plt.plot(centres, Y, "og")
    plt.axis([0, 1, -0.1, 1])
    plt.pause(0.5)

    while t<T: 

        result_LWR1 = schemas_couplage_iteratif(0, x1, taille_voiture, Vmax_LWR1, "fixes", "upwind", sommets_1, centres_1, U_1, dt, insert_LWR1)
        result_FTL = schemas_couplage_iteratif(x1, x2-x1, taille_voiture, Vmax_FTL, "FTL", "upwind", sommets_2, centres_2, U_2, dt, False)
        result_LWR2 = schemas_couplage_iteratif(x2, 1-x2, taille_voiture, Vmax_LWR2, "fixes", "upwind", sommets_3, centres_3, U_3, dt, insert_LWR2)

        U_1 = result_LWR1[0]
        sommets_1 = result_LWR1[1]
        centres_1 = result_LWR1[2]
        dt_1 = result_LWR1[3]

        U_2 = result_FTL[0]
        sommets_2 = result_FTL[1]
        centres_2 = result_FTL[2]
        dt_2 = result_FTL[3]

        U_3 = result_LWR2[0]
        sommets_3 = result_LWR2[1]
        centres_3 = result_LWR2[2]
        dt_3 = result_LWR2[3]
        insert_LWR1 = result_LWR2[4]

        # Calcul du dt 
        dt = min(dt_1, dt_2, dt_3)
            
        result_couplage_1 = coupleur_LWRversFTL(U_1, sommets_1, check_surface, surface_tampon_1, x1)
        transfert_voitures_FTL = result_couplage_1[0]
        surface_tampon_1 = result_couplage_1[1]
        U_transfert_versFTL = result_couplage_1[2]
        
        #result_2 = schemas_couplage()
        #sommets_2 = result_2[0]
        #centres_2 = result_2[1]
        #U_2 = result_2[2]

        centres = [*centres_1, *centres_2, *centres_3]
        U = [*U_1, *U_2, *U_3]
        plt.figure(1)
        plt.clf()
        plt.plot(centres,U)
        Y = [0 for i in range(len(centres))]
        plt.plot(centres, Y, "og")
        plt.axis([0, 1, -0.1, 1])
        plt.pause(0.01)

        # Insertion des voitures dans la zone FTL 
        # Si la zone est saturée on ne fait rien et on attend l'itération suivante pour voir si l'insertion est possible : à implémenter
        if transfert_voitures_FTL > 0:
            result = ajout_voitures_FTL(x1, sommets_2, transfert_voitures_FTL, taille_voiture, U_2, U_transfert_versFTL)
            sommets_2 = result[0]
            centres_2 = 0.5*(sommets_2[:-1] + sommets_2[1:])
            U_2 = result[1]
        t=t+dt

def calcul_aire(f, x1, x2, dx):
    x = x1
    aire = 0
    while x <= x2:
        aire += dx*f(x)
        x += dx
    return aire

main()