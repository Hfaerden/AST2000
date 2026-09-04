


import numpy as np


from ast2000tools.space_mission import SpaceMission as SM
from ast2000tools.solar_system import SolarSystem as SS
import ast2000tools.constants as const
import ast2000tools.utils as utils


seed = utils.get_seed('natanies')

mission = SM(seed)  #setter opp mission-raketten
system = SS(seed)   #setter opp solsystemet vårt

t = 0
dt = 0.1

F = 82144           #kraften som motoren vår gir
fuel_cons = 18.6    #hvor mye drivstoff raketten bruker pr sekund
fuel = 8000         #hvor mye drivstoff vi har med oss
wet_mass = mission.spacecraft_mass + fuel    #massen til HELE raketten, inkl brennstoff
planet_mass = system.masses[0]*1.988*10**30  #masses gitt i solmasser. konverterer til kilo
planet_radius = system.radii[0]*1000         #radii spytter ut i kilometer. konverterer til meter

G_konst = const.G

p_position = np.array([system.initial_positions[0][0], system.initial_positions[1][0]]) * const.AU     #spytter ut posisjonen i AU, konverterer til meter
r_position = np.array([planet_radius, 0]) + p_position  #posisjonen til rakketten vår ved launch
rel_position = r_position-p_position                    #posisjonen til raketten sett fra planeten

rotation_speed = 2*np.pi*planet_radius/( system.rotational_periods[0]*const.day )     #rotasjonshastigheten
v_p = np.array([system.initial_velocities[0][0],system.initial_velocities[1][0] ]) * const.yr/const.AU  #konverterer hastigheten til planeten til SI-enheter
v_rocket = v_p + np.array([0, rotation_speed])       #hastigheten til raketten sett fra solen (inkludert rotasjonshastigheten til planeten)

rel_v = v_rocket-v_p    #hastigheten som sett fra planeten
v_rad_rel = np.dot(rel_v, rel_position) / np.linalg.norm(rel_position)  #hastigheten til raketten radielt utover som sett fra planeten




def integrator (r_mass, r_pos, r_vec, p_pos, v_r_abs, v_rel, v_plan, t, dt, F):

    def gravity_a (r_m, rvec):
        g = (G_konst*(r_m*planet_mass)*rvec)/( np.dot(rvec, rvec)**3 )  #regner ut akselerasjonen fra gravitasjonen
        return g

    
    def motor_a (r_m, rvec):
        a_m = ( rvec / np.linalg.norm(r_vec) )*F/r_m        #regner ut akselerasjonen fra motoren, denne endres med massen til skipet
        
        return a_m  #vi gjør også at akselerasjonsretnigen alltid peker direkte ut av planeten ved hjelp av retningsvektoren til R

    
    a_i = motor_a(r_mass, r_vec) - gravity_a(r_mass, r_vec)          #setter opp første runde av ODE-løseren
    v_rad_rel = np.dot(v_rel, r_vec) / np.sqrt(np.dot(r_vec, r_vec)) #finner hastigheten radielt utover
    
    tot_fuel_cons = 0 
    v_esc = np.sqrt(2*planet_mass*G_konst/planet_radius)        #regner ut unnslipnings-hastigheten

    
    while v_rad_rel < v_esc :   #for å unnslippe må vi at den radielle farten er større enn unnslipningsfarten
        
        t += dt
        
        r_pos += v_rocket*dt + 0.5*a_i*(dt**2)                          #denne blokken er leapfrog-algoritmen (ODE-løseren)
        a_ip1 = motor_a(r_mass, r_vec) - gravity_a(r_mass, r_vec)
        v_r_abs += 0.5*(a_i+a_ip1)*dt
        a_i=a_ip1
        
        
        r_mass = r_mass -fuel_cons*dt   #oppdaterer massen
        tot_fuel_cons += fuel_cons*dt   #oppdaterer menden drivstoff brukt
        
        r_vec = r_pos-p_pos             #oppdaterer vekotren fra planetens midtpunkt til raketten (altså R radielt utover)
        
        
        v_rel = v_r_abs - v_p           #oppdaterer den relative hastigheten
        v_rad_rel = np.dot(v_rel, r_vec) / np.sqrt(np.dot(r_vec, r_vec))    #oppdaterer den radielle hastigheten
 
        
        
        
        if r_mass < mission.spacecraft_mass :
            print(f'break due to out of fuel. r_mass = {r_mass}. ran out at t = {t}. position from the planet is {np.linalg.norm(r_vec)-planet_radius}')
            print(f'missing {v_rad_rel - v_esc} delta v to escape')
            print(v_esc)
            print(v_rad_rel)
            break

        
        
        
    return r_pos, r_mass, tot_fuel_cons, v_r_abs, t
        


print(f'position at t=0 {r_position}, ')

x = integrator(wet_mass, r_position, rel_position, p_position, v_rocket, rel_v, v_p, t, dt, F) 

if x[1] < mission.spacecraft_mass :
    print(f'not enough fuel. missing {x[1]-mission.spacecraft_mass}')



print(f'{x}' )
print(f'fuel used = {x[2]}')






