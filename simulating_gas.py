from re import escape
from warnings import warn

import numpy as np 
from numba import njit
import time


K = 1.38*10**(-23)
T = 3*10**3 
N = 10**5
L = 10**(-6)
m = 3.35*10**(-24)
dt = 10**(-12)
sd = np.sqrt(K*T/m)
timeinit = time.time()

Npos = np.zeros((N, 3), dtype = float)
Npos = Npos - L/2
Nvelocity = np.random.normal(0, sd, (N, 3))
Npos = np.random.uniform(-L/2, L/2, (N, 3))
Npos[0] = [-L/2, -L/2, -L/2]

print("velocities generated at "+ str(time.time() - timeinit))
#print("-----------------")
#print("sum = " + str(np.sum(np.sum(Nvelocity))/10**5))
#Nvelocitysquared = Nvelocity**2
#totalv = []
#for i in Nvelocitysquared:
#totalv.append(np.sqrt(i[0] + i[1] + i[2]))
#print(np.average(totalv))
collisioncount = np.zeros((N,3), dtype = int)
momentum = np.zeros((N, 3), dtype = float) 
escaped = np.zeros(N) 
'''
for i in range(1000):
    Npos = Npos + Nvelocity * dt
    #print(len(Npos))
    for a, b in enumerate(Npos):
        for j, k in enumerate(b):
            if abs(k) > abs(L/2):
                if (Npos[a][2] < -L/2) and (abs(Npos[a][0]) < 0.25*L) and (abs(Npos[a][1]) < 0.25*L):
                    #print("bye fuckers")
                    momentum += Nvelocity[a] * m
                    Npos[a] = [0, 0, L/2]
                    Nvelocity[a] = np.random.normal(0, sd, (1, 3))
                    escapecount += 1
                else:
                    Npos[a][j] = k/abs(k) * L/2 - k/(abs(k)) * (abs(k) - abs(L/2)) #Teleporter den over veggen til riktig sted dersom den gikk over 
                    Nvelocity[a][j] = Nvelocity[a][j] * -1
                    collisioncount += 1
'''
momentumx = np.zeros(N, dtype = float)
momentumy = np.zeros(N, dtype = float) 
momentumz = np.zeros(N, dtype = float)
for i in range(1, 1001):
    Npos = Npos + Nvelocity * dt 
    Noutside = abs(Npos) > L/2
    #Npos = Npos - Noutside*2*(Npos - L/2)
    Nvelocity = Nvelocity - 2*Nvelocity * Noutside
    Nsenteredx = abs(Npos[:,0]) < 0.25*L 
    Nsenteredy = abs(Npos[:,1]) < 0.25*L 
    Centerednoth = np.logical_and(Nsenteredx, Nsenteredy)
    Escape = np.logical_and((Centerednoth), (Npos[:,2] < -L/2)) 
    momentumx += Nvelocity[:,0] * m * Escape
    momentumy += Nvelocity[:,1] * m * Escape
    momentumz += Nvelocity[:,2] * m * Escape
    escaped += Escape 
    collisioncount += Noutside
    if i%20 == 0:
        print("#"*int((i/20)) + "-"*(50-int(i/20)))

print("Simulasjon ferdig på " + str(round(time.time()-timeinit, 2)) + " sekunder \n")
print(f"momentumsum fra unslippende partikler: [{sum(momentumx)}, {sum(momentumy)}, {sum(momentumz)}]")
print("Partikler som slapp ut: " + str(round(sum(escaped)))) 
#print("Antall kollisjoner med vegger: " + str(sum(collisioncount)))
