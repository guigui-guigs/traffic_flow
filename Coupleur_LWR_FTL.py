# Brainstorming 

# Route constituée des zones LWR - FTL - LWR 
# Passage LWR --> FTL 
# Entre chaque zone, on met une "zone de test" dans laquelle on test si la densité atteint S/N
# Le cas échéant, on rajoute une voiture dans le modèle FTL 

# Quelle taille pour la zone de test ? 
# Si la densité n'atteint pas voire jamais la valeur S/N on fait quoi ?
# --> on met la densité de cette zone à 0 et on stocke la valeur.
# --> à chaque pas de temps on réitère le processus en sommant jusqu'à atteindre S/N et ajouter un véhicule dans la partie FTL 
# On peut prendre comme zone de test la taille d'une maille, celle "au contact" de la zone péage

# Architecture du coupleur 
# Une fonction
# Appelée à chaque pas de temps
# Prend la liste des rho_i
# Calcule l'intégrale sur la zone de test 
# Retourne un nombre de véhicules à insérer si besoin 
# Stocke le reste de la surface pour la boucle suivante
# En fait on va retourner une liste à chaque fois [nb de voitures à insérer, intégrale restante, résidu à ajouter à la surface tampon)]

import numpy as np
import matplotlib.pyplot as plt
import math

def coupleur_LWRversFTL(U, sommets, check_surface, surface_tampon, x1):
    # les surfaces initiales sont calculées en dehors de la fonction coupleur, uniquement avec la donnée de rho_0
    L1 = x1 # début de la zone du péage
    #### Couplage LWR --> FTL ####
    i = 0
    length = 0
    while length < 0.3 and i < len(sommets)-1:
        length = sommets[i]
        i+=1
    surface = U[i]*(sommets[i]-sommets[i-1])
    voitures_versFTL = math.floor(surface/check_surface) 
    voitures_surface_tampon = math.floor(surface_tampon/check_surface) 
    surface_tampon = surface_tampon - (voitures_versFTL + voitures_surface_tampon)*check_surface
    U_transfert = U[i]

    return [voitures_versFTL + voitures_surface_tampon, surface_tampon, U_transfert]

def ajout_voitures_FTL(debut_FTL, sommets, nb_voitures, taille_voiture, U, U_transfert):
    new_sommets = sommets[:]
    new_U = U[:]
    saturation = False
    reste_U_transfert = U_transfert
    reste_voitures = nb_voitures
    if len(new_sommets) > 2:
        derniere_voiture = new_sommets[1]
    else :
        derniere_voiture = debut_FTL 
    
    for i in range (1,nb_voitures + 1):
        if len(new_sommets) > 2:
            if (derniere_voiture - 2*i*taille_voiture < debut_FTL):
                saturation = True
            else : 
                new_sommets = np.insert(new_sommets, 1, derniere_voiture - i*taille_voiture)
                reste_voitures = reste_voitures - 1
        else : 
            new_sommets = np.insert(new_sommets, i, debut_FTL + i*taille_voiture)
            reste_voitures = reste_voitures - 1
        if saturation == False :
            new_U.insert(0, U_transfert/nb_voitures) # on répartit la densité sur le nombre de mailles qu'on crée
    reste_U_transfert = reste_U_transfert - (U_transfert/nb_voitures)*(nb_voitures-reste_voitures)
    return [new_sommets, new_U, reste_voitures, reste_U_transfert]

def coupleur_FTLversLWR(U_to_add, U_3): 
    U_3[0] = U_3[0] + U_to_add
    return U_3