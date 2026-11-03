# UTEN KODEMAL

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
pos = system.initial_positions
vel = system.initial_velocities
yearlen = (system.semi_major_axes[0]**(3)/starmass)**(0.5) #Regner ut lengden på et år i jordår, ved Kepler's tredje
timesteps_per_year = 10_000 #Setter hvor mange tidssteg per år 
dt = yearlen/(timesteps_per_year) #Regner ut tidsintervalet vi skal ved dt og årlengden


@jit 
def runsim(time, index):
    timesteps = round(time * timesteps_per_year)
    x1 = pos[0][index] #Henter posisjonene ved t = 0
    y1 = pos[1][index]
    ax1 = -(4*np.pi**2 * starmass / (x1**2+y1**2)) #Regner ut aksellerasjonen ved t = 0, merk at y = 0, og dermed ser uttrykket annerledes ut
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

    return((x[-1], y[-1])) #Returnerer x og y posisjonene 

def solver(A, B, C, D, E, F):
    #regner ut x og y verdiene fra et system med to ligninger 
    matrix = np.array([A, B], [D, E])
    vec = np.array([C, F])
    return(np.linalg.solve(matrix, vec))

def locator(time, measurements):
    #Setter x1, y1 i origo, da vi bruker sola for dette
    x1 = 0 
    y1 = 0 
    #Regner ut posisjonene til planet 1 og 2 (indeksert 0 og 1)
    x2, y2 = runsim(time, 0)
    x3, y3 = runsim(time, 1)
    #Setter A, B, C, D, E, og F verdier etter formellen vi fant
    A = 2*(x2-x1)
    B = 2*(y2-y1)
    C = measurements[-1]**2 - measurements[1]**2+x2**2+y2**2
    D = 2*(x3-x1)
    E = 2*(y2-y1)
    F = measurements[-1]**2 - measurements[2]**2 + x3**2 + y3**2
    (x, y) = solver(A, B, C, D, E, F)
    return(x, y)


def plot_circle (time, measurements):
    #her er measurements en 1d array med tre verdier
    #finner posisjonen for planetene vi vil ha for den tiden:
    x = [ 0, runsim(time, 0)[0], runsim(time, 1)[0] ]
    y = [ 0, runsim(time, 0)[1], runsim(time, 1)[1] ]
    n = 1000 #oppløsning på sirkelen
    meas = [measurements[-1], measurements[0], measurements[1]]
    
    print(meas)
    
    for i in range(len(x)):
        #lager sirkler med sentrum og radius til korresponderende planeter
        R_x = meas[i]*np.cos( np.linspace(0, 2*np.pi, n) ) + x[i]
        R_y = meas[i]*np.sin( np.linspace(0, 2*np.pi, n) ) + y[i]
        plt.plot(R_x, R_y)
        plt.legend(f'sirkel til objekt {i}')
        plt.grid(True)
    
    for i in range(len(x)):
        plt.plot(x[i], y[i], 'ro')
        plt.legend(f'objekt {i}')
        plt.axis('square')
        
    
    plt.show()
        
    
    
    
    




