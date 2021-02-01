# Brainstorming

# A chaque pas de temps on fait tourner le schéma puis le coupleur 

from LWR_new_schema_mailles_fixes_variables import schemas_couplage
from LWR_new_schema_mailles_fixes_variables import rho_0
from Coupleur_LWR_FTL import coupleur_LWRversFTL
from Coupleur_LWR_FTL import coupleur_FTLversLWR

import numpy as np
import matplotlib.pyplot as plt
import math

#########################################################
################### Paramétrisation #####################
#########################################################

dx = 0.01 
dt = 0.004 
simu_duration = 5
T = simu_duration/dt

taille_voiture = dx/10

M1 = 30 # nombre de voitures dans la zone LWR 1
M2 = 15 # zone FTL
M3 = 20 # zone LWR 2

# Définition des différentes zones de circulation entre [0;1]
x1 = 0.3 # début de la zone FTL
x2 = 0.4 # fin de la zone FTL 

N1 = math.floor(x1/dx) 
N2 = math.floor((x2-x1)/dx) 
N3 = math.floor((1-x2)/dx) 

pas_1 = dx*np.ones(N1)
pas_1[0] = 0
pas_2 = dx*np.ones(N2)
pas_2[0] = x1
pas_3 = dx*np.ones(N3+1)
pas_3[0] = x2
sommets_1 = np.cumsum(pas_1) # coordonnées xi+1/2 des extrémités des mailles
sommets_2 = np.cumsum(pas_2)
sommets_3 = np.cumsum(pas_3)
centres_1 = 0.5*(sommets_1[:-1] + sommets_1[1:]) # les coordonnées xi des centres des mailles
centres_2 = 0.5*(sommets_2[:-1] + sommets_2[1:])
centres_3 = 0.5*(sommets_3[:-1] + sommets_3[1:])
U_1 = [0 for i in range(0,N1-1)] 
U_2 = [0 for i in range(0,N2-1)]
U_3 = [0 for i in range(0,N3)]

nb_voitures_LWR1 = 400
surface_LWR1 = calcul_aire(rho_0,0,1,1/1000)
check_surface_1 = surface_LWR1/nb_voitures_LWR1 # surface correspondant à une voiture
surface_tampon_1 = 0

t = 0
centres = [*centres_1, *centres_2, *centres_3]
U = [*U_1, *U_2, *U_3]
plt.figure(1)
plt.clf()
plt.plot(centres,U)
plt.pause(1)

while t<T: 

    # Calcul du dt 
    dt_1 = schemas_couplage(null, null, null, null, N1, sommets_1, U_1, null, True)
    dt_2 = schemas_couplage(null, null, null, null, N2, sommets_2, U_2, null, True)
    dt_3 = schemas_couplage(null, null, null, null, N3, sommets_3, U_3, null, True)
    dt = min(dt_1, dt_2, dt_3)

    result_1 = schemas_couplage("fixes", "parabole", "upwind", taille_voiture, N1, sommets_1, U_1, t, True):
    sommets_1 = result_1[0]
    centres_1 = result_1[1]
    U_1 = result_1[2]
    result_couplage_1 = coupleur_LWRversFTL(U_1, sommets_1, nb_voitures_LWR1, surface_LWR1, surface_tampon, x1)
    transfert_voitures_FTL = result_couplage_1[0]

    result_2 = schemas_couplage()
    sommets_2 = result_2[0]
    centres_2 = result_2[1]
    U_2 = result_2[2]

    result_3 = schemas_couplage()
    sommets_3 = result_3[0]
    centres_3 = 0.5*(sommets_3[:-1] + sommets_3[1:])
    U_3 = result_3[1]

    t=t+dt

def calcul_aire(f, x1, x2, dx):
    x = x1
    aire = 0
    while x <= x2:
        aire += dx*f(x)
        x += dx
    return aire