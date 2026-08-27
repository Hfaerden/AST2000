

from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
seed = utils.get_seed('user_str')



mission = SpaceMission(seed)

print(mission.spacecraft_mass)
print(mission.spacecraft_area)

F = 1
delta_v = 1
v = mission.speed
fuel_cons = 1
speed_boost = v + delta_v
fuel_kg = 1
wet_mass = mission.spacecraft_mass + fuel_kg
dt = 10**(-4)
time = 0



while v < speed_boost :
    
    v += (F*dt)/wet_mass
    fuel_kg -= fuel_cons*dt
    wet_mass = mission.spacecraft_mass + fuel_kg
    
    time += dt

print(f"time needed to achive delta_v {time}")
print()
