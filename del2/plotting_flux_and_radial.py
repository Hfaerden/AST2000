#UTEN KODEMAL!!



import numpy as np
import scipy as sc
import matplotlib.pyplot as plt 



fluxdata = np.load('Relative_flux_curve.npz')
flux = fluxdata["relative_flux"]
t_flux = fluxdata["t"]


plt.plot(t_flux, flux)
plt.ticklabel_format(useOffset=False, axis='y')
plt.show()


rad_vel = np.load('radial_vel_curve2.npz')
vel = rad_vel['radial_vel']
t_rad = rad_vel['t']

plt.plot(t_rad, vel)
plt.ticklabel_format(useOffset=False, axis='y')
plt.show()