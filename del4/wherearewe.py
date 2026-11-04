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

def getpositions(index, time):
    with np.load("orbits.npz") as f:
        times = f['times']
        positions = f['planet_positions']
        index = np.argmin(times - time)
        xval = positions[0][index][time]
        yval = positions[1][index][time]
        return(xval, yval)

def solver(A, B, C, D, E, F):
    #regner ut x og y verdiene fra et system med to ligninger 
    matrix = np.array([A, B], [D, E])
    vec = np.array([C, F])
    return(np.linalg.solve(matrix, vec))

def locator(time, measurements):
    #Setter x1, y1 i origo, da vi bruker sola for dette
    x1 = 0 
    y1 = 0 
    #Henter ut posisjonene til planet 1 og 2 (indeksert 0 og 1)
    x2, y2 = getpositions(time, 0)
    x3, y3 = getpositions(time, 1)
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
    x = [ 0, getpositions(time, 0)[0], getpositions(time, 1)[0] ]
    y = [ 0, getpositions(time, 0)[1], getpositions(time, 1)[1] ]
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

print(getpositions(0, 0))
