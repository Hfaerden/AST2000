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
t = 1000


def launch ():
    dt = 0.005
    F = 128060.01625967553
    fuel_cons = 20.09
    fuel = 11000
    theta = np.pi
    R = system.radii[0]*1000/const.AU * np.array([np.cos(theta), np.sin(theta)])
    launch_position = runsim(t, 0) + R
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
    






