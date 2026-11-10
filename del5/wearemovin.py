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

def interpolation(time):
    with np.load("orbits.nz") as f:
        times = f["times"]
        positions = f["planet_positions"]
        closest = np.argmin(abs(times-time))
        if closest > time:
            t2 = closest 
            t2in = times.index(closest)
            t1in = times.index(t2) - 1
            t1 = times[t1in]
        else:
            t1 = closest 
            t1in = times.index(closest)
            t2in = times.index(t1) + 1
            t2 = times[times[t2in]]

        p1 = positions[:,:,t1in] 
        p2 = positions[:,:,t2in] 
        
        interpolatedpositions = (p1 + (p2-p1)/(t2-t1) * (time-t1))
        return(interpolatedpositions)
