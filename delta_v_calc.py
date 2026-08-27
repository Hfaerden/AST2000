

import numpy as np
from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
seed = utils.get_seed('user_str')



mission = SpaceMission(seed)

print(mission.spacecraft_mass)
print(mission.spacecraft_area)

rocket_thrust_F = 1

fuel_cons = 1
speed_boost = 1
fuel_kg = 1
wet_mass = mission.spacecraft_mass + fuel_kg



def speed_boost (fuel, dv, )


