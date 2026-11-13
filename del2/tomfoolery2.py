
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
starmass = system.star_mass
yearlen = (system.semi_major_axes[0]**(3)/starmass)**(0.5) #Regner ut lengden på et år i jordår, ved Kepler's tredje
pos = system.initial_positions #Henter posisjonene ved t = 0
vel = system.initial_velocities #Henter hastighetene ved t = 0
starmass = system.star_mass #Henter solmassen vår 


yearlen = (system.semi_major_axes[0]**(3)/(system.masses[0]+system.star_mass) )**0.5 #Regner ut lengden på et år i jordår, ved Kepler's tredje
print(yearlen)
timesteps_per_year = 30_000 #Setter hvor mange tidssteg per år 
dt = yearlen/(timesteps_per_year) #Regner ut tidsintervalet vi skal ved dt og årlengden
years = 0 #Setter antall år fra planeten vår vi kjører simulasjonen i
timesteps = round(years * timesteps_per_year) #Regner ut hvor mange tidssteg vi trenger totalt for 20 år 

@jit 
def solve(index):
    x1 = pos[0][index] #Henter posisjonene ved t = 0
    y1 = pos[1][index]
    r = np.sqrt(x1**2 + y1**2) 
    ax1 = (abs(x1/(r)) * 4*np.pi**2 * -np.sign(x1) * starmass / (r**2)) #Regner ut aksellerasjonen ved t = 0, merk at y = 0, og dermed ser uttrykket annerledes ut
    ax2 = (abs(y1/(r)) * 4*np.pi**2 * -np.sign(y1) * starmass / (r**2)) 
    x = [x1]
    y = [y1] #Lager arrays for poisjonene våre 
    ax = [ax1]
    ay = [0]
    vx = [vel[0][index]] #lager arrays for hastighetene våre 
    vy = [vel[1][index]]
    for i in range(timesteps): #Euler-cromer loop for å regne ut aksellerasjon, posisjon, og hastighet
        r = np.sqrt(x[i]**2 + y[i]**2) #Regner ut distansen fra sola 
        ax.append(abs(x[i]/(r)) * 4*np.pi**2 * -np.sign(x[i]) * starmass / (r**2)) #Regner ut aksellerasjon i x retning og legger til i array
        ay.append(abs(y[i]/(r)) * 4*np.pi**2 * -np.sign(y[i]) * starmass / (r**2)) #Regner ut aksellerasjon i y retning og legger til i array
        vx.append(vx[i] + dt*ax[i+1]) #Regner ut hastighet i x retning og legger til i array
        vy.append(vy[i] + dt*ay[i+1]) #Regner ut hastighet i y retning og legger til i array
        x.append(x[i] + dt*vx[i+1]) #Regner ut posisjonen i x retning og legger til i array 
        y.append(y[i] + dt*vy[i+1]) #Regner it posisjonen i y retning og legger til i array 
    return(x, y) #Returnerer x og y posisjonene 

results = []
planet = []
xlist = []
ylist = []

for i in range(len(system.radii)):
    (x, y) = solve(i)
    xlist.append(x)
    ylist.append(y)

positions = [xlist, ylist]

system.verify_planet_positions(timesteps*dt, positions, filename = 'orbits.npz')
