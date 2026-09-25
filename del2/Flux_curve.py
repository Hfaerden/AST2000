#UTEN KODEMAL!!


import numpy as np
import matplotlib.pyplot as plt
from ast2000tools.space_mission import SpaceMission as SM
from ast2000tools.solar_system import SolarSystem as SS
import ast2000tools.constants as const
import ast2000tools.utils as utils


seed = utils.get_seed('natanies')

mission = SM(seed)  #setter opp mission-raketten
system = SS(seed)   #setter opp solsystemet vårt

dt = 100/const.yr
pos = system.initial_positions #Henter posisjonene ved t = 0
yearlen = system.semi_major_axes[0]**(3/2) #Regner ut lengden på et år i jordår, ved Kepler's tredje
i = 5 #planet-index 
radius = system.radii[i]*1000/const.AU #radius til planeten i AU
v_p = np.linalg.norm([system.initial_velocities[0][i], system.initial_velocities[0][i]]) #beholder alt i yr og au, får overflow på liten dt ellers
mass = system.masses[i]
r_s =system.star_radius*1000/const.AU #radius på stjerne i AU
print(system.radii)
#print(mass)
#print(system.semi_major_axes)

F_eclipse = (r_s**2-radius**2)/r_s**2   #finner fluksen som sett av en observatør lang unna i ratio der 1 er fluksen til den udekkede stjernen
values = np.ones(round( (2*r_s-4*radius)/(v_p*dt) ))     #lager en stund der solen er udekket
values = values*np.random.normal( 1, 0.001, len(values) )
values = np.append(values, np.linspace( 1, F_eclipse, round((2*radius)/(v_p*dt))) ) #finner tiden planeten tar å krysse fra delvis formørkelse til total
values = np.append(values, np.ones(round((2*r_s-4*radius)/(v_p*dt)))*F_eclipse )    #Legger til verdiene fluksen får i antall dt det tar å krysse fra ene siden av sole til den andre
values = np.append(values, np.linspace( F_eclipse, 1, round((2*radius)/(v_p*dt))) )  #samme som to linjer over, bare omvendt
values = np.append(values, np.ones(round( (2*r_s-4*radius)/(v_p*dt) )))  #lager en stund der solen er udekket
values = values * np.random.normal( 1, 0.001, len(values) )     #legger til gausisk støy

plt.plot( np.linspace(0, len(values), len(values)) , values )
plt.show()
