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
    x1 = pos[0][index]
    y1 = pos[1][index]
    ax1 = -(4*np.pi**2 * starmass / (x1**2+y1**2))
    #ay1 = -np.sin(np.arctan(y1/x1))* (4*np.pi**2 * starmass / (x1**2+y1**2))
    x = [x1]
    y = [y1]
    ax = [ax1]
    ay = [0]
    vx = [vel[index][0]]
    vy = [vel[index][1]]
    for i in range(timesteps - 1):
        r = np.sqrt(x[i]**2 + y[i]**2)
        ax.append((x[i]/(r)) * 4*np.pi**2 * ((x[i])+10**(-11))/abs(x[i]+10**(-11)) * starmass / (r**2))
        ay.append((y[i]/(r)) * 4*np.pi**2 * ((y[i])+10**(-11))/abs(y[i]+10**(-11)) * starmass / (r**2))
        vx.append(vx[i] + dt*ax[i+1])
        vy.append(vy[i] + dt*ay[i+1])
        x.append(x[i] + dt*vx[i+1])
        y.append(y[i] + dt*vy[i+1])
    return(x)

plt.plot(range(timesteps), solve(0))
plt.show()
