#UTEN KODEMAL!!

import numpy as np
import matplotlib.pyplot as plt
from ast2000tools.space_mission import SpaceMission as SM
from ast2000tools.solar_system import SolarSystem as SS
from ast2000tools.shortcuts import SpaceMissionShortcuts
import ast2000tools.constants as const
import ast2000tools.utils as utils
import wherearewe
import referencevelocity 

seed = utils.get_seed('natanies')

mission = SM(seed)  #setter opp mission-raketten
system = SS(seed)   #setter opp solsystemet vårt
t = 0.19
theta = np.pi


def launch_r ():
    dt = 0.005
    F = 128060.01625967553 #henter launch-verdiene fra del 1
    fuel_cons = 20.09
    fuel = 11000
    R = system.radii[0]*1000/const.AU * np.array([np.cos(theta), np.sin(theta)]) #finer hvor på planeten vi er
    launch_position, l_time = wherearewe.getpositions_time(t, 0) #finner posisjonen iforhold til solen ved t 
    launch_position += R
    print(launch_position)
    mission.set_launch_parameters(F, fuel_cons, fuel, 500, launch_position, l_time)
    mission.launch_rocket(dt) 
    shortcuts = SpaceMissionShortcuts(mission, [33759])
    
    launch = shortcuts.get_launch_results()
    mission.verify_launch_result(launch[2])
    
    return launch


def visualize_trilateration(info):  
    d = mission.measure_distances()
    wherearewe.plot_circle(t, d, info[2])
    


def test_reference_vel(info):
    #lager en testfunksjon som sjekker om hastighetene vi måler er de samme som vi får fra AST 2000
    dop = np.array(mission.measure_star_doppler_shifts())   #henter dopplerverdiene
    #dop = dopplershifts = mission.star_doppler_shifts_at_sun
    vel = referencevelocity.referencevelocitytosun(-dop[0], -dop[1]) #kallerløsningsfunksjonen vår
    tol = 0.001
    vel = vel*const.yr/const.AU
    
    if abs(np.linalg.norm(vel - np.array(info[3])) ) < tol :
           print(f'passed! difference in velocity is {abs(np.linalg.norm(vel - info[3]))} au/yr')
    else:
        print(f'failed :( poopoohead! difference was {abs(np.linalg.norm(vel - info[3]))}')
        print(f'individual values: from doppler {vel}, from ast2000 {info[3]} ')




info = launch_r()
t = t + 492/const.yr #time after launch when we want our position
visualize_trilateration(info)
test_reference_vel(info)
d = mission.measure_distances()
print(info[2], info[3])
print(f'diff launch-snarvei avstand til planet og measure_distance i AU: {np.linalg.norm(info[2]-wherearewe.getpositions(t, 0) ) -d[0]}')
mission.verify_manual_orientation(wherearewe.locator(t, d), info[3], np.pi)



