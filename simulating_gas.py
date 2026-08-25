import numpy as np 
from numba import njit
import time


K = 1.38*10**-23
T = 3*10**3 
N = 10**5
L = 10*10**(-6)
m = 3.35*10**(-24)
dt = 10**(-12)
sd = np.sqrt(K*T/m)
timeinit = time.time()
def randomv():
    return(np.random.default_rng().normal(0, sd))

Npos = np.zeros((N, 3), dtype = float)
Npos = Npos - L/2
Nvelocity = np.zeros((N, 3), dtype = float)

for i in range(len(Nvelocity)):
    for j in range(3):
        Nvelocity[i][j] = randomv()
 
Npos[0] = [-L/2, -L/2, -L/2]
for i in range (1, len(Npos)):
    if Npos[i-1][0] + L/(N**(1/3)) <= L/2:
        Npos[i][0] = Npos[i-1][0] + L/(N**(1/3))
        Npos[i][1] = Npos[i-1][1]
        Npos[i][2] = Npos[i-1][2]
    elif Npos[i-1][1] + L/(N**(1/3)) <= L/2:
        Npos[i][0] = -L/2
        Npos[i][1] = Npos[i-1][1] + L/(N**(1/3))
        Npos[i][2] = Npos[i-1][2]
    else:
        Npos[i][0] = -L/2
        Npos[i][1] = -L/2
        Npos[i][2] = Npos[i-1][2] + L/(N**(1/3))

print("velocities generated at "+ str(time.time() - timeinit))
#print(Npos)
#print("-----------------")
#print(Nvelocity)
#print("sum = " + str(np.sum(np.sum(Nvelocity))/10**5))

#Nvelocitysquared = Nvelocity**2

#totalv = []
#for i in Nvelocitysquared:
#totalv.append(np.sqrt(i[0] + i[1] + i[2]))
#print(np.average(totalv))

#@njit 
#def update(timestep):
collisioncount = 0
for i in range(1000):
    Npos = Npos + Nvelocity * dt
    #print(len(Npos))
    for a, b in enumerate(Npos):
        for j, k in enumerate(b):
            if abs(k) > abs(L/2):
                Npos[a][j] = k/abs(k) * L/2 - k/(abs(k)) * (abs(k) - abs(L/2)) #Teleporter den over veggen til riktig sted dersom den gikk over 
                Nvelocity[a][j] = Nvelocity[a][j] * k/abs(k)
                collisioncount += 1
    print("stage " + str(i) + " complete")
print(collisioncount)
