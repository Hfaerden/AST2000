# UTEN KODEMAL

import numpy as np
import matplotlib.pyplot as plt
from ast2000tools.space_mission import SpaceMission as SM
from ast2000tools.solar_system import SolarSystem as SS
import ast2000tools.constants as const
import ast2000tools.utils as utils
import wearemovin as yargh



seed = utils.get_seed('natanies')

mission = SM(seed)  #setter opp mission-raketten
system = SS(seed)   #setter opp solsystemet vårt

print(system.radii[1])
print( ((0.001)/(system.rotational_periods[0]/const.day)) )#/(60*60) )
print(system.rotational_periods[0])

print(f'esc_vel: {np.sqrt((2*const.*system.masses[1]*system.star_mass/(system.radii[1]*1000))}')

'''
t = 0.18

pos = yargh.interpolation(t)

p1 = np.array([ pos[0][0], pos[1][0] ])
p2 = np.array([ pos[0][1], pos[1][1] ])


print(p2)
print(np.linalg.norm(p1-p2))

pos2 = yargh.interpolation(t+0.1)
p1 = np.array([ pos2[0][0], pos2[1][0] ])
p2 = np.array([ pos2[0][1], pos2[1][1] ])
print(f'p1 {p1} p2: {p2}')


'''



