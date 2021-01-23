import numpy as np
import matplotlib.pyplot as plt


##################### test ####################






###############################################



l = 0.01
Mtot = 0.5*(2-l+np.sin(np.pi*(1-l))/np.pi)
N = int(Mtot/l)
Vm = 1
dt = 0.75*l/Vm
Tfin = 0.9

Z = []

for j in range(0, N+1):
    if j == 0:
        z = -1+l
    else:
        z = Z[-1] + l

    err = 10
    limit = 10**(-8)
    while err > limit:
        z_old = z
        z = z - (z + np.sin(np.pi*z)/np.pi - 2*j*l + 1)/(1 + np.cos(np.pi*z))
        err = abs(z - z_old)/abs(z)

    Z += [z]



plt.figure(1)
X = np.linspace(-1, 1, N+1)
Densite = (np.cos(np.pi*X)+1)/2
plt.scatter(Z,np.zeros(len(Z)))
plt.plot(X, Densite)
#plt.show()

### Etape 2 : recurrence
t = 0
Y = np.array([(Z[i+1] - Z[i])/l for i in range(0, N)])
Y = np.append(Y, [(2 - Z[-1] + Z[0])/l])

while t < Tfin:
    t = t+dt
    Yold = Y
    Vy = np.array([Vm * max(0, 1-1/Y[i]) for i in range(0, len(Y))])
    Z = Z+dt*Vy
    Y = np.array([(Z[i + 1] - Z[i]) / l for i in range(0, N)])
    Y = np.append(Y, [(2 - Z[-1] + Z[0])/l])

    #print(len(X), len(Y), len(Z), len(Vy), len(Densite))

    Densite = (Densite*Yold)/Y

    if Z[-1] >= 1:
        Z = np.concatenate(([Z[-1]-2], Z[0:-1]))
        Y = np.concatenate(([Y[-1]], Y[0:-1]))
        Densite = np.concatenate(([Densite[-1]], Densite[0:-1]))

    #print(Z)

def findIndexDensite(x, Z):
    for i in range(0, len(Z)-1):
        if Z[i] <= x < Z[i + 1]:
            return i


Zp = np.concatenate((Z[1:], [Z[0]+2]))
coord = (Z + Zp)/2
Np = 2000
Dens = np.zeros((Np, 1))
Discretisation = np.linspace(Z[0], Z[-1], Np)
for n in range(0, Np-1):
    i = findIndexDensite(Discretisation[n], Z)
    Dens[n] = Densite[i]
Dens[-1] = Densite[-1]

plt.figure(2)
plt.plot(coord, Densite, 'x')
plt.plot(Discretisation, Dens, 'r')
plt.xlabel("x")
plt.ylabel("densite(t,x)")

JJ = [i for i in range(0, len(Z))]
plt.figure(3)
plt.scatter(JJ, Z)
plt.xlabel("indice Lagrange point")
plt.ylabel("position Lagrange point")

plt.figure(4)
plt.scatter(Z, t*np.ones(len(Z)))
plt.xlabel("position")
plt.ylabel("temps")

plt.show()








print("Well done")












