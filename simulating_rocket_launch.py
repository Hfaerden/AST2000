


import numpy as np
from ast2000tools.SpaceMission import SpaceMission as SM
import ast2000tools.utils as utils
seed = utils.get_seed('user_str')

mission = SM(seed)



thrust = 1
fuel_cons = 1
fuel = 1
wet_mass = mission.spacecraft_mass + fuel
planet_mass = 1
planet_radius = 1
G_konst = 6.67 * 10**(-11)
planet_posisjon = np.zeros(2)
rel_v = np.zeros(N, 3)


rocket_position = np.array([planet_radius, 0])


def gravity (r_mass, r_pos, p_pos):
    r_vec = r_pos-p_pos
    R = np.linalg.norm(r_vec, axis=1)
    G = G_konst*(r_mass*planet_mass)*r_vec/R**3
    return G



while rel_v < v_esc :
    




