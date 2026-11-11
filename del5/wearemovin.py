#Uten kodemal :)
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

masses = system.masses

def interpolation(time):
    with np.load("orbits.npz") as f:
        times = f["times"]
        positions = f["planet_positions"]
        closest = np.argmin(abs(times-time))
        if times[closest] > time:
            t2in = closest 
            t2 = times[closest]
            t1in = t2in - 1
            t1 = times[t1in]
        else:
            t1in = closest 
            t1 = times[closest]
            t2in = t1in +1 
            t2 = times[t2in]

        p1 = positions[:,:,t1in] 
        p2 = positions[:,:,t2in] 
        
        interpolatedpositions = (p1 + (p2-p1)/(t2-t1) * (time-t1))
        return(interpolatedpositions)

def webemovin(starttime, startpos, startv, timestepsN, dt):
    x = np.zeros((timestepsN, 2))
    x[0] = np.array([startpos])
    v = np.zeros((timestepsN, 2))
    v[0] = np.array([startv])
    initpos = interpolation(starttime)
    r0 = np.sqrt(initpos[0]**2 + initpos[1]**2)
    a0comp = [(-initpos[0]/r0) * (2*np.pi)**2 * masses/(r0**2), (-initpos[1]/r0) * (2*np.pi)**2 * masses/(r0**2)]
    a0 = [np.sum(a0comp[0]), np.sum(a0comp[1])]
    a = np.zeros((timestepsN, 2))
    a[0] = np.array(a0)
    t = [starttime]

    for i in range(timestepsN-1):
        positions = interpolation(t[i])
        r = np.sqrt((positions[0]-x[i][0])**2 + (positions[1]-x[i][1])**2)
        acomponents = [(-(positions[0]-x[i][0])/r) * (2*np.pi)**2 * masses/(r**2), -((positions[1]-x[i][1])/r) * (2*np.pi)**2 * masses/(r**2)]
        a[i+1] = [np.sum(acomponents[0]), np.sum(acomponents[1])]
        v[i+1] = (a[i+1] * dt) + v[i]
        x[i+1] = (v[i+1] * dt) + x[i]
        t.append(t[i] + dt)

    return((starttime + timestepsN * dt), x[-1], v[-1])

print(webemovin(0.3, [-0.2438, 0.1800], [-6.99, -7.39], 10_000, 1/(10_000)))
