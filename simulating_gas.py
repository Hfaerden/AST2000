import numpy as np 


K = 1.38*10**-23
T = 3*10**3 
N = 100
L = 10*10**(-6)
error = 10^(-10)
dt = 10^(-12)

Npos = np.zeros((100, 3), dtype = float)
Npos = Npos - L/2
Nvelocity = []

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

print(Npos)

'''
for i in range(5):
    for j in range(5):
        for k in range(4):
            Npos.append([(2*10**(-6)-10*10**(-6))/2], [2*10**(-6)-10*10**(-6))/2], 2.5*10*10**(-6)])

rng = np.random.default_rng()

sample = rng.normal()

for i in range(1000):
    Npos = Npos + Nvelocity * dt 
    for i in range(3):
        for n in Npos:
            if abs(Npos[n][i]) - 5*10**(-7) < error: 
                '''
   #flip angle 
