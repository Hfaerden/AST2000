#UTEN KODEMAL!!



import numpy as np
import scipy as sc
import matplotlib.pyplot as plt 



fluxdata = np.load('Flux_curve_data.npz')
flux = fluxdata["flux"]
t_flux = fluxdata["time"]


plt.plot(t_flux, flux)
plt.ticklabel_format(useOffset=False, axis='y')
plt.show()


rad_vel = np.load('')

