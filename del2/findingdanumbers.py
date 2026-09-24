#IKKE BRUKT KODEMAL

import numpy as np
import matplotlib.pyplot as plt
from ast2000tools.space_mission import SpaceMission as SM
from ast2000tools.solar_system import SolarSystem as SS
import ast2000tools.constants as const
import ast2000tools.utils as utils


seed = utils.get_seed('natanies')

mission = SM(seed)  #setter opp mission-raketten
system = SS(seed)   #setter opp solsystemet vårt
f = np.linspace(0, 2*np.pi, 100000)
r = np.array(list(np.zeros(len(f)) for i in range(len(system.radii))))
e = system.eccentricities
a = system.semi_major_axes
p = a*(1-e**2)
x_analytisk =  []
y_analytisk = []
print(system.radii)
for i in range(len(system.radii)):
    f_i = f + (np.pi-system.aphelion_angles[i])
    r[i] = (p[i]/(1+e[i]*np.cos(f_i)))

plt.axes(projection="polar")
for j in range(len(r)):
    plt.polar(f, r[j], label = (f"planet {j}"))
plt.show()

