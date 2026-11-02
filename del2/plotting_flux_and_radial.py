#UTEN KODEMAL!!



import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from ast2000tools.space_mission import SpaceMission as SM
from ast2000tools.solar_system import SolarSystem as SS
import ast2000tools.constants as const
import ast2000tools.utils as utils 


def sinus(x, v, p, t0):
    return v*np.sin(2*np.pi/p*x-t0)
    

star_mass = 3.3389946609527805 #fått fra samarbeidspartner
G = const.G

fluxdata = np.load('Relative_flux_curve.npz') #henter ut verdiene for flux vi fikk
flux = fluxdata["relative_flux"]
t_flux = fluxdata["t"]


plt.plot(t_flux, flux)
plt.xlabel('tid i år')
plt.ylabel('Flux relativt til normal output')
plt.grid(True)
plt.ticklabel_format(useOffset=False, axis='y')
plt.show()


rad_vel = np.load('radial_vel_curve 6.npz') #henter ut verdiene for radiell hastighet
vel = rad_vel['radial_vel']
t_rad = rad_vel['t']
vrel_npec = vel-np.mean(vel) #fjerner relativ hastighet for curve_fit


param, p_co = sc.optimize.curve_fit(sinus, t_rad, vrel_npec, p0 = [0.001, 170, 50] )    #finner gode parametere for modellen vår
ans = sinus(t_rad, param[0], param[1], param[2]) + np.mean(vel)     #lager array for plotting

print(param[0])
print(param[1])
print(param[2])
print(f'peculiar = {np.mean(vel)}')

print(f'massen til planeten er: {(param[1]/(2*np.pi*G))**(1/3)}')

plt.plot(t_rad, vel, label = 'samlet data')
plt.plot(t_rad, ans, label = 'modell')
plt.xlabel('tid i år')
plt.ylabel('hastighet i AU/yr')
plt.grid(True)
plt.legend()
plt.ticklabel_format(useOffset=False, axis='y')
plt.show()