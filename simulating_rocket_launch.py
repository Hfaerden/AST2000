


import numpy as np


from ast2000tools.space_mission import SpaceMission as SM
from ast2000tools.solar_system import SolarSystem as SS
import ast2000tools.utils as utils


seed = utils.get_seed('natanies')

mission = SM(seed)  #setter opp mission-raketten
system = SS(seed)   #setter opp solsystemet vårt

t = 0
dt = 1

F = 82144
fuel_cons = 18.6
fuel = 8000
wet_mass = mission.spacecraft_mass + fuel    #massen til HELE raketten, inkl brennstoff
planet_mass = system.masses[0]*1.988*10**30  #masses gitt i solmasser. konverterer til kilo
planet_radius = system.radii[0]*1000         #radii spytter ut i kilometer. konverterer til meter

G_konst = 6.67 * 10**(-11)
print(wet_mass)


p_position = np.zeros(2)
r_position = np.array([planet_radius, 0]) + p_position #SM.launch_position
rel_position = r_position-p_position
rotation_speed = 0

v_p = np.zeros(2)
v_rocket = v_p + np.array([0, rotation_speed])



rel_v = v_rocket-v_p
v_rad = np.dot(v_rocket, rel_position) / np.sqrt(rel_position[0]**2+rel_position[1]**2)
v_rad_rel = np.dot(rel_v, rel_position) / np.sqrt(rel_position[0]**2+rel_position[1]**2)

print(np.sqrt(rel_position[0]**2+rel_position[1]**2))
print(rel_position[1])


def integrator (r_mass, r_pos, r_vec, p_pos, v_r_abs, v_rel, v_plan, t, dt, F):
    2

    def gravity_a (r_m, rvec):
        G = (G_konst*(r_m*planet_mass)*rvec)/( np.dot(rvec, rvec)**3 )
        return G

    
    def motor_a (r_m, rvec):
        a = ( rvec / np.linalg.norm(r_vec) )*F/r_m
        
        return a

    
    a_i = motor_a(r_mass, r_vec) - gravity_a(r_mass, r_vec)
    v_rad_rel = np.dot(v_rel, r_vec) / np.sqrt(np.dot(r_vec, r_vec))
    
    tot_fuel_cons = 0
    v_esc = np.sqrt(2*planet_mass*G_konst/planet_radius)
    plot = []
    
    while v_rad_rel < v_esc :
        
        t += dt
        
        r_pos += v_rocket*dt + 0.5*a_i*(dt**2)
        a_ip1 = motor_a(r_mass, r_vec) - gravity_a(r_mass, r_vec)
        v_r_abs += 0.5*(a_i+a_ip1)*dt
        a_i=a_ip1
        print(v_rad_rel)
        
        #a_i += motor_a(r_mass, r_vec) - gravity_a(r_mass, r_vec)
        #v += a_i*dt
        
        
        r_mass = r_mass -fuel_cons*dt
        tot_fuel_cons += fuel_cons*dt
        
        r_vec = r_pos-p_pos
        
        
        v_rel = v_r_abs - v_p
        v_rad_rel = np.dot(v_rel, r_vec) / np.sqrt(np.dot(r_vec, r_vec))
        v_esc = np.sqrt(2*planet_mass*G_konst/planet_radius)
        
        
        
        if r_mass < mission.spacecraft_mass :
            print(f'break due to out of fuel. r_mass = {r_mass}. ran out at t = {t}. position from the planet is {np.linalg.norm(r_vec)-planet_radius}')
            print(f'missing {v_rad_rel - v_esc} delta v to escape')
            print(v_esc)
            print(v_rad_rel)
            break
        #print(np.linalg.norm(a_i))
        #print(np.linalg.norm(r_vec))
        #print(v_rad_rel)
        #if np.absolute(a_i) <= 0:
         #   print('not enough force to liftoff')
        
        
        
    return r_pos, r_mass, tot_fuel_cons, v_r_abs, t, plot
        


print(f'position at t=0 {r_position}, ')
print()

x = integrator(wet_mass, r_position, rel_position, p_position, v_rocket, rel_v, v_p, t, dt, F) 

if x[1] < mission.spacecraft_mass :
    print(f'not enough fuel. missing {x[1]-mission.spacecraft_mass}')



print(f'{x}' )






