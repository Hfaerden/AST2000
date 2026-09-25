#UTEN KODEMAL!!

import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from ast2000tools.space_mission import SpaceMission as SM
from ast2000tools.solar_system import SolarSystem as SS
import ast2000tools.constants as const
import ast2000tools.utils as utils
from numba import jit 

seed = utils.get_seed('natanies')

mission = SM(seed)  #setter opp mission-raketten
system = SS(seed)   #setter opp solsystemet vårt

k = 5 #planet nr 5
pos = system.initial_positions #Henter posisjonene ved t = 0
vel = system.initial_velocities #Henter hastighetene ved t = 0
starmass = system.star_mass #Henter solmassen vår 
planetmass = system.masses
p_pos = system.initial_positions[k]

yearlen = system.semi_major_axes[0]**(3/2) #Regner ut lengden på et år i jordår, ved Kepler's tredje
#print(yearlen)
timesteps_per_year = 10000 #Setter hvor mange tidssteg per år 
dt = yearlen/(timesteps_per_year) #Regner ut tidsintervalet vi skal ved dt og årlengden
#print(dt)
t_end = 100
timesteps = round(t_end * timesteps_per_year) #Regner ut hvor mange tidssteg vi trenger totalt for t år
CM = 

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
    #r_meter = np.sqrt( (x[0]**2*const.AU) + (y[0]*const.AU)**2)
    #E_tot = [ 0.5*planetmass[index]*const.m_sun*((vx[0]*const.AU/const.yr)**2+(vy[0]*const.AU/const.yr)**2)**2 - const.G*starmass*planetmass[index]*const.m_sun**2/r_meter ]
       
    for i in range(timesteps - 1): #Euler-cromer loop for å regne ut akselerasjon, posisjon, og hastighet
        r = np.sqrt(x[i]**2 + y[i]**2) #Regner ut distansen fra sola 
        ax.append(abs(x[i]/(r)) * 4*np.pi**2 * -np.sign(x[i]) * starmass / (r**2)) #Regner ut aksellerasjon i x retning og legger til i array
        ay.append(abs(y[i]/(r)) * 4*np.pi**2 * -np.sign(y[i]) * starmass / (r**2)) #Regner ut aksellerasjon i y retning og legger til i array
        vx.append(vx[i] + dt*ax[i+1]) #Regner ut hastighet i x retning og legger til i array
        vy.append(vy[i] + dt*ay[i+1]) #Regner ut hastighet i y retning og legger til i array
        x.append(x[i] + dt*vx[i+1]) #Regner ut posisjonen i x retning og legger til i array 
        y.append(y[i] + dt*vy[i+1]) #Regner it posisjonen i y retning og legger til i array 
        
        #r_meter = np.sqrt( (x[i]**2*const.AU) + (y[i]*const.AU)**2)
        #E_tot.append( 0.5*planetmass[index]*const.m_sun*((vx[i]*const.AU/const.yr)**2+(vy[i]*const.AU/const.yr)**2)**2 - const.G*starmass*planetmass[index]*const.m_sun**2/(r_meter) )
        
        
    return([x, y], [vx, vy])#, E_tot) #Returnerer x og y posisjonene , samt vx og vy
