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

print(seed)
mission = SM(seed)  #setter opp mission-raketten
system = SS(seed)   #setter opp solsystemet vårt
t = 0.18
flight_time = 0.07
timesteps = 400000
dt = flight_time/timesteps
dtheta = 0.0001   #setter endringen i vinkelen for hver launch
n_thetaoffsets = 30  #setter antall launch-er vi skal ha
closest_point_str = []






def launch_r (theta): #samme launch-kode som del4
    dt = 0.005
    F = 128060.01625967553 #henter launch-verdiene fra del 1
    fuel_cons = 20.09
    fuel = 11000
    R = system.radii[0]*1000/const.AU * np.array([np.cos(theta), np.sin(theta)]) #finer hvor på planeten vi er
    pos = yargh.interpolation(t)
    launch_position = np.array([pos[0][0], pos[1][0]]) #finner posisjonen for planeten iforhold til solen ved t 
    print(launch_position)
    launch_position += R    #legger til radiusen til planeten i retningen bestemt av theta
    print(launch_position)  
    mission.set_launch_parameters(F, fuel_cons, fuel, 500, launch_position, t) #setter launch-parameterene med verdier fra del1
    mission.launch_rocket(dt) 
    shortcuts = SpaceMissionShortcuts(mission, [33759])
    
    launch = shortcuts.get_launch_results()
    mission.verify_launch_result(launch[2]) #kjører ast2000tools sin launch-kode
    print(launch[2], launch[3])
    
    return launch




def interpolation_plan(time, planned_traj):
    upper = np.searchsorted(times, time) #Finner nærmeste indeks til input tiden, rundet opp, som en upper bound
    lower = upper - 1 #Tar forrige indeks som lower bound 
    t2, t1 = times[upper], times[lower] #Finner tidene ved indeksene til bounds 
    p1 = positions[:,:,lower] #Henter posisjonene ved tidene til bounds 
    p2 = positions[:,:,upper] 
    interpolatedpositions = (p1 + (p2-p1)/(t2-t1) * (time-t1)) #Regner ut lineær interpolerte posisjoner
    return(interpolatedpositions) #Returnerer posisjonene





    
pos = yargh.interpolation(t)
p1 = np.array([ pos[0][0], pos[1][0] ])
p2 = np.array([ pos[0][1], pos[1][1] ])
    
base_theta_offset = 0.645 #
theta = np.arctan((p2[1]-p1[1])/(p2[0]-p1[0])) + base_theta_offset











