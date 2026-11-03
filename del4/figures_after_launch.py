#UTEN KODEMAL!!

import numpy as np
import matplotlib.pyplot as plt
from ast2000tools.space_mission import SpaceMission as SM
from ast2000tools.solar_system import SolarSystem as SS
from ast2000tools.shortcuts import SpaceMissionShortcuts
import ast2000tools.constants as const
import ast2000tools.utils as utils
from numba import jit 
from wherearewe import plot_circle
from wherearewe import runsim
import referencevelocity 

seed = utils.get_seed('natanies')

mission = SM(seed)  #setter opp mission-raketten
system = SS(seed)   #setter opp solsystemet vårt
t = 0


def launch ():
    dt = 0.005
    F = 128060.01625967553 #henter launch-verdiene fra del 1
    fuel_cons = 20.09
    fuel = 11000
    theta = np.pi
    R = system.radii[0]*1000/const.AU * np.array([np.cos(theta), np.sin(theta)]) #finer hvor på planeten vi er
    launch_position = runsim(t, 0) + R #finner posisjonen iforhold til solen
    print(launch_position)
    
    mission.set_launch_parameters(F, fuel_cons, fuel, 500, launch_position,  t)
    mission.launch_rocket(dt) 
    shortcuts = SpaceMissionShortcuts(mission, [33759])

    launch = shortcuts.get_launch_results()
    mission.verify_launch_result(launch[2])


def visualize_trilateration():  
    d = mission.measure_distances()
    plot_circle(t, d)

def test_reference_vel():
    dop = np.array(mission.measure_star_doppler_shifts())
    angles = mission.star_direction_angles
    vel = referencevelocity.referencevelocitytosun(dop[0], dop[1])
    tol = 0.01
    
    if abs(np.linalg.norm(vel - launch[3])) < tol :
           print(f'passed! difference in velocity is {abs(np.linalg.norm(vel - launch[3]))}')
    else:
        print(f'failed :( poopoohead! difference was {abs(np.linalg.norm(vel - launch[3]))}')
        print(f'individual values: from doppler {vel}, from ast2000 {launch[3]} ')

launch()
test_reference_vel()





