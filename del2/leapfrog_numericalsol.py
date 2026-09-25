#UTEN KODEMAL

import numpy as np
import matplotlib.pyplot as plt
from ast2000tools.space_mission import SpaceMission as SM
from ast2000tools.solar_system import SolarSystem as SS
import ast2000tools.constants as const
import ast2000tools.utils as utils
from numba import jit 


seed = utils.get_seed('natanies')

mission = SM(seed)  #setter opp mission-raketten
system = SS(seed)   #setter opp solsystemet vårt

pos = system.initial_positions #Henter posisjonene ved t = 0
vel = system.initial_velocities #Henter hastighetene ved t = 0
starmass = system.star_mass #Henter solmassen vår 


yearlen = system.semi_major_axes[0]**(3/2) #Regner ut lengden på et år i jordår, ved Kepler's tredje
print(yearlen)
timesteps_per_year = 10000 #Setter hvor mange tidssteg per år 
dt = yearlen/(timesteps_per_year) #Regner ut tidsintervalet vi skal ved dt og årlengden
timesteps = round(200 * timesteps_per_year) #Regner ut hvor mange tidssteg vi trenger totalt for 20 år 

@jit 
def solve(index):
    x1 = pos[0][index] #Henter posisjonene ved t = 0
    y1 = pos[1][index]
    ax1 = -(4*np.pi**2 * starmass / (x1**2+y1**2)) #Regner ut aksellerasjonen ved t = 0, merk at y = 0, og dermed ser uttrykket annerledes ut
    x = [x1]
    y = [y1] #Lager arrays for poisjonene våre 
    ax = [ax1]
    ay = [0]
    vx = [vel[0][index]] #lager arrays for hastighetene våre 
    vy = [vel[1][index]]
    for i in range(timesteps - 1): #Leapfrog loop for å regne ut aksellerasjon, posisjon, og hastighet
        r = np.sqrt( x[i]**2 + y[i]**2 ) #Regner ut distansen fra sola
        ax.append(abs( x[i]/(r)) * 4*np.pi**2 * -np.sign(x[i]) * starmass / (r**2) ) #Regner ut aksellerasjon i x retning og legger til i array
        ay.append(abs( y[i]/(r)) * 4*np.pi**2 * -np.sign(y[i]) * starmass / (r**2) ) #Regner ut aksellerasjon i y retning og legger til i array
        x.append( x[i] + dt*vx[i] + 0.5*ax[i]*dt**2 ) #Regner ut posisjonen i x retning og legger til i array 
        y.append( y[i] + dt*vy[i] + 0.5*ay[i]*dt**2)  #Regner it posisjonen i y retning og legger til i array
        vx.append( vx[i] + 0.5*dt*(ax[i+1]+ax[i]) ) #Regner ut hastighet i x retning og legger til i array
        vy.append( vy[i] + 0.5*dt*(ay[i+1]+ay[i])) #Regner ut hastighet i y retning og legger til i array

    return([x, y]) #Returnerer x og y posisjonene 

results = []

for i in range(len(system.radii)):
    results.append(solve(i))
for j in range(len(results)):
    plt.plot(results[j][0], results[j][1])
#plt.plot(result[0], result[1])
print(len(system.radii))
plt.show()
