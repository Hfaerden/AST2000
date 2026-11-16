#UTEN KODEMAL!!

import numpy as np
import matplotlib.pyplot as plt
from ast2000tools.space_mission import SpaceMission as SM
from ast2000tools.solar_system import SolarSystem as SS
from ast2000tools.shortcuts import SpaceMissionShortcuts
import ast2000tools.constants as const
import ast2000tools.utils as utils
import wearemovin as yargh


seed = utils.get_seed('natanies')

mission = SM(seed)  #setter opp mission-raketten
system = SS(seed)   #setter opp solsystemet vårt
t = 0.18
flight_time = 0.06
timesteps = 50000
dt = flight_time/timesteps



def getpositions_time(index, time):
    with np.load("orbits.npz") as f:
        times = f['times']
        positions = f['planet_positions']
        time = np.argmin(abs(times - time))
        xval = positions[0][index][time]
        yval = positions[1][index][time]
        print(times[time])
        return(xval, yval), times[time]

pos = yargh.interpolation(t)
p1 = np.array([ pos[0][0], pos[1][0] ])
p2 = np.array([ pos[0][1], pos[1][1] ])

theta_offset = 0.35
theta = np.arctan((p1[1]-p2[1])/(p1[0]-p2[0]))+theta_offset
    
print(theta)

def launch_r (theta):
    dt = 0.005
    F = 128060.01625967553 #henter launch-verdiene fra del 1
    fuel_cons = 20.09
    fuel = 11000
    R = system.radii[0]*1000/const.AU * np.array([np.cos(theta), np.sin(theta)]) #finer hvor på planeten vi er
    pos = yargh.interpolation(t)
    launch_position = np.array([pos[0][0], pos[1][0]]) #finner posisjonen iforhold til solen ved t 
    print(launch_position)
    launch_position += R
    print(launch_position)
    mission.set_launch_parameters(F, fuel_cons, fuel, 500, launch_position, t)
    mission.launch_rocket(dt) 
    shortcuts = SpaceMissionShortcuts(mission, [33759])
    
    launch = shortcuts.get_launch_results()
    mission.verify_launch_result(launch[2])
    
    return launch

pee = launch_r(theta)
time_used_during_launch = pee[1]/const.yr

trajectory = yargh.webemovin(t+time_used_during_launch, pee[2], pee[3], timesteps, dt)
planet_positions_start = yargh.interpolation(t)
planet_positions_end = yargh.interpolation(t+flight_time)

closest_point = trajectory[4]



plt.plot(trajectory[1][:, 0], trajectory[1][:, 1])
plt.plot(planet_positions_start[0][0], planet_positions_start[1][0], 'ro')
plt.plot(planet_positions_start[0][1], planet_positions_start[1][1], 'ro')
plt.plot(planet_positions_end[0][0], planet_positions_end[1][0], 'ro', color = 'blue')
plt.plot(planet_positions_end[0][1], planet_positions_end[1][1], 'ro', color = 'blue')
plt.legend(['rakett-bane','posisjon før' ,'','posisjon_etter'])
plt.axis('square')
plt.xlabel('avstand i AU')
plt.ylabel('avstand i AU')
plt.title(f'rakettutskytning ved t = {t} år, med lengde på flyving = {flight_time}år, og sikter {theta_offset} radianer foran')
plt.show()



