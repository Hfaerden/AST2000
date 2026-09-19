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

pos = system.initial_positions
vel = system.initial_velocities 
starmass = system.star_mass
timesteps = 20 * 10_000
dt = 1/10_000

@jit
def solve(index):
    x1 = pos[index][0]
    y1 = pos[index][1]
    ax1 = -(4*np.pi**2 * starmass / x1**2)
    #ay1 = -(y1 + 4*np.pi**2 * starmass / y1**2)
    x = [x1]
    y = [y1]
    ax = [ax1]
    ay = [0]
    vx = [vel[index][0]]
    vy = [vel[index][1]]
    for i in range(timesteps - 1):
        ax.append((-abs(x[i])/x[i]) * 4*np.pi**2 * starmass / x[i]**2) if x[i] != 0 else 0
        ay.append((-abs(y[i])/y[i]) * 4*np.pi**2 * starmass  / y[i]**2) if y[i] != 0 else 0 
        vx.append(vx[i] + dt*ax[i+1])
        vy.append(vy[i] + dt*ay[i+1])
        x.append(x[i] + dt*vx[i+1])
        y.append(y[i] + dt*vy[i+1])
    return(x)

plt.plot(range(timesteps), solve(0))
plt.show()
