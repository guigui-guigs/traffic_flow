#########################################################
#########################################################
############# LWR - schémas eulériens ###################
#########################################################
#########################################################
import numpy as np
import matplotlib.pyplot as plt
import math

def rho(x):
    return x

def f(rho):
    return rho(1-rho)

N = 50
Z = []
X = np.linspace(0,1,N+1)

for i in range(0,N+1):
    Z.append(f(rho(x)))

plt.plot()