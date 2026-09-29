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
pos = system.initial_positions #Henter posisjonen ved t = 0
vel = system.initial_velocities #Henter hastigheten ved t = 0
starmass = system.star_mass #Henter solmassen vår 
planetmass = system.masses[k]
p_pos = np.array([system.initial_positions[0][k], system.initial_positions[1][k]])
#print(p_pos)

yearlen = system.semi_major_axes[0]**(3/2) #Regner ut lengden på et år i jordår, ved Kepler's tredje
#print(yearlen)
timesteps_per_year = 10000 #Setter hvor mange tidssteg per år 
dt = yearlen/(timesteps_per_year) #Regner ut tidsintervalet vi skal ved dt og årlengden
#print(d++t)
t_end = 40
timesteps = round(t_end * timesteps_per_year) #Regner ut hvor mange tidssteg vi trenger totalt for t år

CM = planetmass*p_pos/(starmass+planetmass)     #Finner posisjon og hastighet til CM
CMvx = float(planetmass*vel[0][k]/(starmass+planetmass))
CMvy = float(planetmass*vel[1][k]/(starmass+planetmass))

v_s = [float(-CMvx), float(-CMvy)]      #Konverterer til CM-ref systemet
v_p = [ [float(vel[0][k]-CMvx)], [float(vel[1][k]-CMvy)] ]
p_pos -= CM
s_pos = -CM
print(CMvx, CMvy)
print(v_s)


#print(CM)



@jit
def solve(pos_px, pos_py, pos_sx, pos_sy, v_px, v_py, v_sx, v_sy):
    x1_p, y1_p = pos_px, pos_py    #Henter posisjonene ved t = 0
    x1_s, y1_s = pos_sx, pos_sy
    #CMx, CMy = [mu*(x1_p-x1_s)/starmass], [mu*(y1_p-y1_s)/planetmass]
    r = np.sqrt((-x1_p+x1_s)**2 + (-y1_p+y1_s)**2)    #relative posisjonen 
    ax_p, ay_p = [-((x1_p-x1_s) * 4*np.pi**2 * starmass / (r**3) )], [-((y1_p-y1_s) * 4*np.pi**2 * -np.sign(y1_p) * starmass / (r**3) )]    #Regner ut aksellerasjonen ved t = 0
    ax_s, ay_s = [-((x1_s-x1_p) * 4*np.pi**2 * starmass / (r**3) )], [-((y1_s-y1_p) * 4*np.pi**2 * -np.sign(y1_s) * starmass / (r**3) )]
    x_p, y_p = [x1_p], [y1_p]   #Lager arrays for poisjonene våre  
    x_s, y_s = [x1_s], [y1_s]
    vx_p, vy_p = [v_px], [v_py]  #lager arrays for hastighetene 
    vx_s, vy_s = [v_sx], [v_sy] 
    
     
    for i in range(timesteps - 1): #Euler-cromer loop for å regne ut akselerasjon, posisjon, og hastighet
    
        r = np.sqrt((x_p[i]-x_s[i])**2 + (y_p[i]-y_s[i])**2)
        rx_p2s = (x_p[i]-x_s[i])
        ry_p2s = (y_p[i]-y_s[i])
        
        ax_p.append(-( rx_p2s * 4*np.pi**2 * starmass / (r**3) ))
        ay_p.append(-( ry_p2s * 4*np.pi**2  * starmass / (r**3) ))
        vx_p.append(vx_p[i]+ax_p[i+1]*dt)
        vy_p.append(vy_p[i]+ay_p[i+1]*dt)
        x_p.append(x_p[i]+vx_p[i+1]*dt)
        y_p.append(y_p[i]+vy_p[i+1]*dt)
        
        ax_s.append(-(-rx_p2s * 4*np.pi**2  * planetmass / (r**3) ))
        ay_s.append(-(-ry_p2s * 4*np.pi**2 * planetmass / (r**3) ))
        vx_s.append(vx_s[i]+ax_s[i+1]*dt)
        vy_s.append(vy_s[i]+ay_s[i+1]*dt)
        x_s.append(x_s[i]+vx_s[i+1]*dt)
        y_s.append(y_s[i]+vy_s[i+1]*dt)
        
        #CMx.append( mu*(x1_p-x1_s)/starmass )
        #CMy.append( planetmass)
        
        
    return([[x_p, y_p], [vx_p, vy_p]], [[x_s, y_s], [vx_s, vy_s]]) #, [CMx, CMy])#, E_tot) #Returnerer x og y posisjonene , samt vx og vy


print(v_p[0])
print(v_p[1])
results = solve(p_pos[0], p_pos[1], s_pos[0], s_pos[1], v_p[0], v_p[1],  v_s[0], v_s[1])
plt.plot(results[0][0][0], results[0][0][1])
plt.plot(results[1][0][0], results[1][0][1])
plt.show()

