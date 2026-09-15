


import numpy as np
import matplotlib.pyplot as plt
from simulating_gas import Gassimulation
from ast2000tools.space_mission import SpaceMission as SM
from ast2000tools.solar_system import SolarSystem as SS
import ast2000tools.constants as const
import ast2000tools.utils as utils


seed = utils.get_seed('natanies')

mission = SM(seed)  #setter opp mission-raketten
system = SS(seed)   #setter opp solsystemet vårt


def integrator (r_mass, r_pos, r_vec, p_pos, v_r_abs, v_rel, v_plan, t, dt, F):     #integrator-funksjonen som inneholder både gravitasjon, akselerasjon, og ODE-løseren

    def gravity_a (rvec):
        g = (G_konst*planet_mass*rvec)/( np.linalg.norm(r_vec)**3 )  #regner ut akselerasjonen fra gravitasjonen
        return g

    
    def motor_a (r_m, rvec):
        a_m = ( rvec / np.linalg.norm(r_vec) )*F/r_m        #regner ut akselerasjonen fra motoren, denne endres med massen til skipet
        
        return a_m  #vi gjør også at akselerasjonsretnigen alltid peker direkte ut av planeten ved hjelp av retningsvektoren til R

    
    #a_i = motor_a(r_mass, r_vec) - gravity_a(r_vec)          #setter opp første runde av ODE-løseren
    v_rad_rel = np.dot(v_rel, r_vec) / np.linalg.norm(r_vec)         #finner hastigheten radielt utover
    
    tot_fuel_cons = 0 
    v_esc = np.sqrt(2*planet_mass*G_konst/planet_radius)        #regner ut unnslipnings-hastigheten

    
    while np.linalg.norm(v_rel) < v_esc :   #for å unnslippe må vi at den radielle farten er større enn unnslipningsfarten
        
        
        #r_pos += v_rocket*dt + 0.5*a_i*(dt**2)                          #denne blokken er leapfrog-algoritmen (ODE-løseren)
        #a_ip1 = motor_a(r_mass, r_vec) - gravity_a(r_vec)
        #v_r_abs += 0.5*(a_i+a_ip1)*dt
        #a_i=a_ip1
        
        t += dt
        a_ip1 = motor_a(r_mass, r_vec) - gravity_a(r_vec)              #Euler-cromer gir ca samme svar
        v_r_abs += a_ip1*dt
        r_pos += v_r_abs*dt
        
        p_pos += v_p*dt     #oppdaterer planetens posisjon
        
        r_mass = r_mass -fuel_cons*dt   #oppdaterer massen
        tot_fuel_cons += fuel_cons*dt   #oppdaterer mengden drivstoff brukt
        
        r_vec = r_pos-p_pos             #oppdaterer vekotren fra planetens midtpunkt til raketten (altså R radielt utover)
        
        
        v_rel = v_r_abs - v_p           #oppdaterer den relative hastigheten
        v_esc = np.sqrt(2*planet_mass*G_konst/np.linalg.norm(r_vec))
        
 
        global r_array
        r_array = np.append(r_array, [r_vec], axis=0)   #appender posisjonen i forhold til planeten for å plotte
        
        if np.linalg.norm(r_vec) < planet_radius :
            print(f'underground :D needs { r_mass*np.linalg.norm(motor_a(r_mass, r_vec) - gravity_a(r_vec)) } Newtons of thrust')
            break
        
        
        if r_mass < mission.spacecraft_mass :   #sjekker om vi går tom for drivstoff (negativ masse gir rare resultater)
            print(f'break due to out of fuel. r_mass = {r_mass}. ran out at t = {t}. position from the planet is {np.linalg.norm(r_vec)-planet_radius}')
            print(f'missing {v_rad_rel - v_esc} delta v to escape')
            print(v_esc)
            print(v_rad_rel)
            break


    return r_pos, r_pos-p_pos, r_mass, tot_fuel_cons, v_r_abs, t
        
        #resten printer info om rakettens forsjellige paramentre





t = 0
dt = 0.005


box_area = ( 10**(-6) )**2
n_box = mission.spacecraft_area/box_area
model = Gassimulation(round(10**5), 3.5*10**3, const.m_H2, 10**(-9), 10**(-12), 10**(-6))
model.runsim()
f_box = model.Forcez
partikkel_masse = 3.36*10**(-27)

F = n_box*f_box        #kraften som motoren vår gir
fuel_cons = (sum(model.escaped)*partikkel_masse/10**(-9) )*n_box   #hvor mye drivstoff raketten bruker pr sekund
print(fuel_cons)
fuel = 11000        #hvor mye drivstoff vi har med oss
wet_mass = mission.spacecraft_mass + fuel    #massen til HELE raketten, inkl brennstoff
planet_mass = system.masses[0]*1.988*10**30  #masses gitt i solmasser. konverterer til kilo
planet_radius = system.radii[0]*1000         #radii spytter ut i kilometer. konverterer til meter

print(f'planet mass = {planet_mass}')
print(f' planet radius = {planet_radius}')

G_konst = const.G

p_position = np.array( [system.initial_positions[0][0], system.initial_positions[1][0]] ) * const.AU     #spytter ut posisjonen i AU, konverterer til meter
r_position = np.array([planet_radius, 0]) + p_position  #posisjonen til rakketten vår ved launch
rel_position = r_position-p_position                    #posisjonen til raketten sett fra planeten


rotation_speed = 2*np.pi*planet_radius/( system.rotational_periods[0]*const.day )     #rotasjonshastigheten
v_p = np.array([system.initial_velocities[0][0],system.initial_velocities[1][0] ]) * const.yr/const.AU  #konverterer hastigheten til planeten til SI-enheter
v_rocket = v_p + np.array([0, rotation_speed])       #hastigheten til raketten sett fra solen (inkludert rotasjonshastigheten til planeten)
print(v_rocket)

rel_v = v_rocket-v_p    #hastigheten som sett fra planeten
v_rad_rel = np.dot(rel_v, rel_position) / np.linalg.norm(rel_position)  #hastigheten til raketten radielt utover som sett fra planeten

r_array = np.array([rel_position])  #gjør klar til plotting av posisonen




 # finner ikke bug-en, har lett lenge. Posisjonen er rett, tiden er ca rett, men koden min gir feil plassering av raketten
 # og den sier at det er for lite drivstoff til å komme til orbit


print(f'position at t=0 {r_position}, ')

print(wet_mass)
x = integrator(wet_mass, r_position, rel_position, p_position, v_rocket, rel_v, v_p, t, dt, F) 

print(f'array of values = {x}' )
print(f'fuel used = {wet_mass - x[2]}')

plt.plot(r_array[:,0], r_array[:,1])    #sjekker for plotter som gir fysisk mening
plt.show()

mission.set_launch_parameters(F, fuel_cons, fuel, x[-1], r_position/const.AU,  t)
mission.launch_rocket(dt)
mission.verify_launch_result(x[0])
