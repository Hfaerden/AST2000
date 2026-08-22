import numpy as np 


K = 1.38*10**-23
T = 3*10**3 
N = 100
L = 10**(-6)
error = 10^(-10)
dt = 10^(-12)

Npos = np.array[x, y, z]
Nvelocity = np.array[x, y, z] 


rng = np.random.default_rng()

sample = rng.normal()

for i in range(1000):
    Npos = Npos + Nvelocity * dt 
    for i in range(3):
        for n in Npos:
            if abs(Npos[n][i]) - 5*10**(-7) < error: 
                #flip angle 
