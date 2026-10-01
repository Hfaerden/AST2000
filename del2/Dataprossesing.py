
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from ast2000tools.space_mission import SpaceMission as SM
from ast2000tools.solar_system import SolarSystem as SS
import ast2000tools.constants as const
import ast2000tools.utils as utils
from numba import jit 


flux_data = np.load("\home\nataniel\ast\kode\AST2000\del2\Flux_curve_data")

radial_vel = flux_data["arr_0"]
t = flux_data["arr_1"]

plt.plot(radial_vel, t)
plt.show()

