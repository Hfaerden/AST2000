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
tidligere_beregnet_launch_varighet = 0.22621271691949735
flight_time =  0.07#tidligere_beregnet_launch_varighet - t
timesteps = 400000
dt = flight_time/timesteps
dtheta = 0   #setter endringen i vinkelen for hver launch
n_thetaoffsets = 1 #setter antall launch-er vi skal ha
closest_point_str = []


def getpositions_time(index, time): #Henter posisjonen for planetene fra npz-filen til bruk senere til intrepoleringen
    with np.load("orbits.npz") as f:
        times = f['times']
        positions = f['planet_positions']
        time = np.argmin(abs(times - time))
        xval = positions[0][index][time]
        yval = positions[1][index][time]
        print(times[time])
        return(xval, yval), times[time]



def spaceflight(theta_offset): #kjører launch- og interplanetarisk- simulering for oppskytning ved vinkel theta 
    
     
    print(theta)
    
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
    
    pee = launch_r(theta)   #kjører launch for vinkelen theta
    time_used_during_launch = pee[1]  
    
    print(pee[3])
    #kjører løseren fra wherearewe.py med verdiene vi har regnet ut
    trajectory = yargh.webemovin(time_used_during_launch, pee[2], pee[3], timesteps, dt)
    planet_positions_start = yargh.interpolation(pee[1])
    planet_positions_end = yargh.interpolation(t+flight_time)
    
    planet_traj = trajectory[3]
    
    distance_str = []
    for i in range(len(planet_traj)):
        current_p_pos = ([planet_traj[i][0][1],planet_traj[i][1][1] ]) 
        distance_str.append(np.linalg.norm(current_p_pos - trajectory[1][i]) )
    
    min_distance = np.min(distance_str)
    print(min_distance)
    print(f'argmin: {np.argmin(distance_str)}')
    print(min_distance)
    closest_point_time = np.argmin(distance_str)*dt + time_used_during_launch
    print(closest_point_time)
    abs_r = np.linalg.norm(trajectory[1][np.argmin(distance_str)])
    l = abs_r * np.sqrt(system.masses[1]/(system.star_mass*10))
    print (f'max distance for orbit insertion for our closest approach is {l}, while we are {min_distance} at closest')
    
    return min_distance, closest_point_time, trajectory, planet_positions_start, planet_positions_end


    
pos = yargh.interpolation(t)
p1 = np.array([ pos[0][0], pos[1][0] ])
p2 = np.array([ pos[0][1], pos[1][1] ])
    
base_theta_offset = 0.645
theta = np.arctan((p2[1]-p1[1])/(p2[0]-p1[0])) + base_theta_offset


closest_point = []
flights = []

for i in range(n_thetaoffsets):
    theta = i*dtheta + base_theta_offset
    current_flight = spaceflight(theta)
    flights.append(current_flight)
    closest_point.append(current_flight[0])


closest_point = np.array(closest_point)
closest_point = np.where( closest_point < (system.radii[1]+500)*1000/const.AU, closest_point+10, closest_point  )
closest_flight_index = np.argmin(closest_point)  
  
distance, close_time, trajectory, planet_positions_start, planet_positions_end = flights[closest_flight_index]

print(f'closest viable point: {np.min(closest_point)}')
best_theta = closest_flight_index*dtheta + base_theta_offset
print(f'best theta = {best_theta} , end theta = {base_theta_offset + dtheta*n_thetaoffsets}')
print(f'time at closest approach: {close_time}')
print(f'closest point list {closest_point}')
print(f'speed at closest approach: {trajectory[2][int( (close_time-t)/dt)] }')
print(f'point min {trajectory[1][int( (close_time-t)/dt)]}')


distance, close_time, trajectory, planet_positions_start, planet_positions_end = flights[closest_flight_index]

plt.plot(trajectory[1][:, 0], trajectory[1][:, 1])
plt.plot(planet_positions_start[0][0], planet_positions_start[1][0], 'ro')
plt.plot(planet_positions_start[0][1], planet_positions_start[1][1], 'ro')
plt.plot(planet_positions_end[0][0], planet_positions_end[1][0], 'o', color = 'blue')
plt.plot(planet_positions_end[0][1], planet_positions_end[1][1], 'o', color = 'blue')
plt.plot(0, 0, 'o', color='black' )
plt.legend(['rakett-bane','posisjon før' ,'','posisjon_etter'])
plt.axis('square')
plt.xlabel('avstand i AU')
plt.ylabel('avstand i AU')
plt.title(f'rakettutskytning ved t = {t} år, med lengde på flyving = {flight_time}år, og sikter {best_theta} radianer foran')
plt.show()



