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

def coupleur_LWRversFTL(densite, sommets, surface_LWR1, check_surface_1, surface_tampon, x1):
    # les surfaces initiales sont calculées en dehors de la fonction coupleur, uniquement avec la donnée de rho_0
    L1 = 0.3 # début de la zone du péage
    #### Couplage LWR --> FTL ####
    i = 0
    length = 0
    while length < 0.3 and i < len(sommets)-1:
        length = sommets[i]
        i+=1
    surface = densite[i]*(sommets[i]-sommets[i-1])
    voitures_versFTL = math.floor(surface/check_surface_1) 
    voitures_versFTL_tampon = math.floor(surface_tampon/check_surface_1)

    return [voitures_versFTL + voitures_versFTL_tampon, surface_LWR1-surface, surface-voitures_versFTL*check_surface_1]

def coupleur_FTLversLWR(x2): 
    L2 = 0.4 # fin de la zone du péage
    return null