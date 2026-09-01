


import numpy as np
from ast2000tools.space_mission import SpaceMission as SM
import ast2000tools.utils as utils
seed = utils.get_seed('natanies')

mission = SM(seed)


t = 0
dt = 60

F = 500000
fuel_cons = 1
fuel = 1
wet_mass = mission.spacecraft_mass + fuel
planet_mass = 1
planet_radius = 1
G_konst = 6.67 * 10**(-11)



p_position = np.zeros(2)
r_position = np.array([planet_radius, 0]) + p_position #SM.launch_position
rel_position = r_position-p_position
rotation_speed = 0

v_p = np.zeros(2)
v_esc = np.sqrt(2*planet_mass*G_konst/planet_radius)
v_rocket = v_p + np.array([0, rotation_speed])



rel_v = v_rocket-v_p
v_rad = np.dot(v_rocket, rel_position) / np.sqrt(rel_position[0]**2+rel_position[1]**2)
v_rad_rel = np.dot(rel_v, rel_position) / np.sqrt(rel_position[0]**2+rel_position[1]**2)

print(np.sqrt(rel_position[0]**2+rel_position[1]**2))
print(rel_position[1])


def integrator (r_mass, r_pos, r_vec, p_pos, v_r_abs, v_rel, v_plan, t, dt, F):
    2

    def gravity_a (r_m, rvec):
        G = G_konst*(r_m*planet_mass)*rvec/np.sqrt(np.dot(rvec, rvec))**3
        return G

    
    def motor_a (r_m, rvec):
        print(np.sqrt(rvec[0]**2+rvec[1]**2))
        a = ( rvec / np.sqrt(np.dot(rvec, rvec)) )*F/r_m
        
        return a

    print(np.sqrt(r_vec[0]**2+r_vec[1]**2))
    a_i = motor_a(r_mass, r_vec) - gravity_a(r_mass, r_vec)
    v_rad_rel = np.dot(v_rel, r_vec) / np.sqrt(np.dot(r_vec, r_vec))
    
    while v_rad_rel < v_esc :
        
        t += dt
        r_pos += v_rocket*dt + 0.5*a_i*(dt**2)
        a_ip1 = motor_a(r_mass, r_vec) - gravity_a(r_mass, r_vec)
        v_r_abs = 0.5*(a_i+a_ip1)*dt
        a_i=a_ip1
        
        r_mass = r_mass -fuel_cons*dt
        r_vec = r_pos-p_pos
        
        v_rel = v_r_abs - v_p
        v_rad_rel = np.dot(v_rel, r_vec) / np.sqrt(np.dot(r_vec, r_vec))
        
        #print(t)
        
        
    return r_pos, v_r_abs, t
        

print(f'position at t=0 {r_position}, ')
print()

x = integrator(wet_mass, r_position, rel_position, p_position, v_rocket, rel_v, v_p, t, dt, F) 

print(f'{x}' )






