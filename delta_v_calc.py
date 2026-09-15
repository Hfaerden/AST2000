from ast2000tools.space_mission import SpaceMission
from simulating_gas import Gassimulation
import ast2000tools.constants as const
import ast2000tools.utils as utils
seed = utils.get_seed('user_str')



mission = SpaceMission(seed)

print(mission.spacecraft_mass)
print(mission.spacecraft_area)

model = Gassimulation(round(10**5.1), 3.5*10**3, const.m_H2, 10**(-9), 10**(-12), 10**(-7))
model.runsim()
f_box = model.Forcez
n_box = mission.spacecraft_area/(model.len)**2

F = n_box*f_box
print(F)
delta_v = 1
v: float = 1
fuel_cons = (sum(model.escaped)*const.m_H2*n_box)
speed_boost = v + delta_v
fuel_start: float = 1
fuel_kg = fuel_start
wet_mass = mission.spacecraft_mass + fuel_kg
dt = 10**(-4)
time = 0


while v < speed_boost :
    print(v)
    v += (F*dt)/wet_mass
    fuel_kg -= fuel_cons*dt
    wet_mass = mission.spacecraft_mass + fuel_kg

    time += dt
    
    if fuel_kg < 0:
        break


if fuel_kg > 0:
   print(f"time needed to achive delta_v: {time} amout of fuel needed {fuel_start - fuel_kg} ")

else:
    print(f"ran out of fuel. it would take {-fuel_kg}kg of more fuel to do that manouver")


