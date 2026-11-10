# UTEN KODEMAL

import numpy as np
import matplotlib.pyplot as plt
from ast2000tools.space_mission import SpaceMission as SM
from ast2000tools.solar_system import SolarSystem as SS
import ast2000tools.constants as const
import ast2000tools.utils as utils

seed = utils.get_seed('natanies')

mission = SM(seed)  #setter opp mission-raketten
system = SS(seed)   #setter opp solsystemet vårt
starmass = system.star_mass
pos = system.initial_positions
vel = system.initial_velocities
yearlen = (system.semi_major_axes[0]**(3)/starmass)**(0.5) #Regner ut lengden på et år i jordår, ved Kepler's tredje
timesteps_per_year = 10_000 #Setter hvor mange tidssteg per år 
dt = yearlen/(timesteps_per_year) #Regner ut tidsintervalet vi skal ved dt og årlengden


def getpositions(time, index):
    #henter planetbanene fra ast2000tools
    with np.load("orbits.npz") as f:
        times = f['times']
        positions = f['planet_positions']
        time = np.argmin(abs(times - time))
        xval = positions[0][index][time]
        yval = positions[1][index][time]
        return(xval, yval)


def getpositions_time(time, index):
    #henter planetbanene fra ast2000tools
    with np.load("orbits.npz") as f:
        times = f['times']
        positions = f['planet_positions']
        time = np.argmin(abs(times - time))
        xval = positions[0][index][time]
        yval = positions[1][index][time]
        return (xval, yval), times[time]


def solver(A, B, C, D, E, F):
    #regner ut x og y verdiene fra et system med to ligninger 
    matrix = np.array([ [A, B], [D, E] ])
    vec = np.array([C, F])
    print(f'matrix: {matrix}')
    print(vec)
    return(np.linalg.solve(matrix, vec))

def locator(time, measurements):
    #Setter x1, y1 i origo, da vi bruker sola for dette
    x1 = 0 
    y1 = 0 
    #Regner ut posisjonene til planet 1 og 2 (indeksert 0 og 1)
    x2, y2 = getpositions(time, 0)
    x3, y3 = getpositions(time, 1)
    #Setter A, B, C, D, E, og F verdier etter formellen vi fant
    A = 2*(x2-x1)
    B = 2*(y2-y1)
    C = measurements[-1]**2 - measurements[0]**2+x2**2+y2**2
    D = 2*(x3-x1)
    E = 2*(y3-y1)
    F = measurements[-1]**2 - measurements[1]**2 + x3**2 + y3**2
    (x, y) = solver(A, B, C, D, E, F)
    print (f'pos = {x, y}')
    return ([x, y])


def plot_circle (time, measurements, rocketpos):
    #her er measurements en 1d array med tre verdier
    #finner posisjonen for planetene vi vil ha for den tiden:
    x = [ 0, getpositions(time, 0)[0], getpositions(time, 1)[0] ]
    y = [ 0, getpositions(time, 0)[1], getpositions(time, 1)[1] ]
    n = 10000000 #oppløsning på sirkelen
    meas = [measurements[-1], measurements[0], measurements[1]]
    calc_pos = locator(time, measurements)
    #print(meas)
    
    for i in range(len(x)):
        #lager sirkler med sentrum og radius til korresponderende planeter
        R_x = meas[i]*np.cos( np.linspace(0, 2*np.pi, n) ) + x[i]
        R_y = meas[i]*np.sin( np.linspace(0, 2*np.pi, n) ) + y[i]
        plt.plot(R_x, R_y)
        #plt.legend(f'sirkel til objekt {i}')
        plt.grid(True)
    
    for i in range(len(x)):
        plt.plot(x[i], y[i], 'ro')
        #plt.legend(f'objekt {i}')
        
        
    plt.plot(rocketpos[0], rocketpos[1], 'bo')
    plt.plot(calc_pos[0], calc_pos[1],'ko')
    plt.axis('square')
    plt.xlabel('avstand i AU')
    plt.ylabel('avstand i AU')
    plt.show()
        
    
    
    
    




