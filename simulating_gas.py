import numpy as np 
import time
import math

K = 1.38*10**(-23)
T = 3*10**3 
N = 10**5
L = 10**(-6)
m = 3.35*10**(-24)
dt = 10**(-12)
total_time=10**-9
timesteps = int(total_time/dt)
rocket_width = 2

sd = np.sqrt(K*T/m)
timeinit = time.time()
nchambers = (rocket_width/L)**2

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
#collisioncount = np.zeros((N,3), dtype = int)
momentum = np.zeros((N, 3), dtype = float) 
escaped = np.zeros(N) 
momentumx = np.zeros(N, dtype = float)
momentumy = np.zeros(N, dtype = float) 
momentumz = np.zeros(N, dtype = float)



for i in range(1, timesteps+1):
    
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
    #collisioncount += Noutside
    if i%20 == 0:
        print("#"*int((i/20)) + "-"*(50-int(i/20)))

Forcex = nchambers*(sum(momentumx)/total_time)
Forcey = nchambers*(sum(momentumy)/total_time)
Forcez = nchambers*(sum(momentumz)/total_time)

kg_pr_s = nchambers*((round(sum(escaped))*m)/total_time)

print(f"drivstoff brukt i kg/s er {kg_pr_s}")
print(f"Kraften i x retning er {Forcex}")
print(f"Kraften i y retning er {Forcey}")
print(f"Kraften i z retning er {Forcez}")

print("Simulasjon ferdig på " + str(round(time.time()-timeinit, 2)) + " sekunder \n")
print(f"momentumsum fra unslippende partikler: [{sum(momentumx)}, {sum(momentumy)}, {sum(momentumz)}]")
print("Partikler som slapp ut: " + str(round(sum(escaped))))
print(f"gjennomsnittlig hastighet for unsluppne partikler i z-retning er {(sum(momentumz)/m)/round(sum(escaped))}")





