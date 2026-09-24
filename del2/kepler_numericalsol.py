
#KODEMAL ER IKKE BRUKT

import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from ast2000tools.space_mission import SpaceMission as SM
from ast2000tools.solar_system import SolarSystem as SS
import ast2000tools.constants as const
import ast2000tools.utils as utils
from numba import jit 


seed = utils.get_seed('natanies')

mission = SM(seed)  #setter opp mission-raketten
system = SS(seed)   #setter opp solsystemet vårt

pos = system.initial_positions #Henter posisjonene ved t = 0
vel = system.initial_velocities #Henter hastighetene ved t = 0
starmass = system.star_mass #Henter solmassen vår 
planetmass = system.masses

yearlen = system.semi_major_axes[0]**(3/2) #Regner ut lengden på et år i jordår, ved Kepler's tredje
#print(yearlen)
timesteps_per_year = 10000 #Setter hvor mange tidssteg per år 
dt = yearlen/(timesteps_per_year) #Regner ut tidsintervalet vi skal ved dt og årlengden
#print(dt)
t_end = 100
timesteps = round(t_end * timesteps_per_year) #Regner ut hvor mange tidssteg vi trenger totalt for t år

@jit 
def solve(index):
    x1 = pos[0][index] #Henter posisjonene ved t = 0
    y1 = pos[1][index]
    ax1 = -(4*np.pi**2 * starmass / (x1**2+y1**2)) #Regner ut aksellerasjonen ved t = 0, merk at y = 0, og dermed ser uttrykket annerledes ut
    x = [x1]
    y = [y1] #Lager arrays for poisjonene våre 
    ax = [ax1]
    ay = [0]
    vx = [vel[0][index]] #lager arrays for hastighetene våre 
    vy = [vel[1][index]]
    #r_meter = np.sqrt( (x[0]**2*const.AU) + (y[0]*const.AU)**2)
    #E_tot = [ 0.5*planetmass[index]*const.m_sun*((vx[0]*const.AU/const.yr)**2+(vy[0]*const.AU/const.yr)**2)**2 - const.G*starmass*planetmass[index]*const.m_sun**2/r_meter ]
       
    for i in range(timesteps - 1): #Euler-cromer loop for å regne ut akselerasjon, posisjon, og hastighet
        r = np.sqrt(x[i]**2 + y[i]**2) #Regner ut distansen fra sola 
        ax.append(abs(x[i]/(r)) * 4*np.pi**2 * -np.sign(x[i]) * starmass / (r**2)) #Regner ut aksellerasjon i x retning og legger til i array
        ay.append(abs(y[i]/(r)) * 4*np.pi**2 * -np.sign(y[i]) * starmass / (r**2)) #Regner ut aksellerasjon i y retning og legger til i array
        vx.append(vx[i] + dt*ax[i+1]) #Regner ut hastighet i x retning og legger til i array
        vy.append(vy[i] + dt*ay[i+1]) #Regner ut hastighet i y retning og legger til i array
        x.append(x[i] + dt*vx[i+1]) #Regner ut posisjonen i x retning og legger til i array 
        y.append(y[i] + dt*vy[i+1]) #Regner it posisjonen i y retning og legger til i array 
        
        #r_meter = np.sqrt( (x[i]**2*const.AU) + (y[i]*const.AU)**2)
        #E_tot.append( 0.5*planetmass[index]*const.m_sun*((vx[i]*const.AU/const.yr)**2+(vy[i]*const.AU/const.yr)**2)**2 - const.G*starmass*planetmass[index]*const.m_sun**2/(r_meter) )
        
        
    return([x, y], [vx, vy])#, E_tot) #Returnerer x og y posisjonene , samt vx og vy




def kepler(pos_array, v_array, index):  
    x = pos_array[0]
    y = pos_array[1]
    l_ap = 0       #lengden planetene har bevegd seg ila tiden
    l_per = 0
    v_ap = []
    v_per = []
    A_ap = []       #Lagrer arealene, de skal alle være like store, da dt er like stor
    A_per = []
    n_t_ap = 0      #antall timesteps på periapsis og apoapsis, for å regne gjennomsnittshastighet
    n_t_per = 0
    r = []
    print(index)
    
    for i in range(len(x)):
        r.append(np.sqrt(x[i]**2 + y[i]**2))    #konverterer til polarkoordinater
        
    ap_i = r.index(min(r))   #finner tidspunktet når planeten er nærmest solen
    for j in range(ap_i-100, ap_i+100):   #beregner for posisjonene inærheten av periapsen, spesifikt dt*100 på hver side
        A_ap.append( 0.5*np.linalg.norm( np.cross((x[j], y[j], 0), (x[j+1], y[j+1], 0)) ) ) #tar kryssproduktet for å regne arealet til trekanten ved apoapsen
        v_ap.append(np.linalg.norm([ v_array[0][j], v_array[1][j] ])) #sjekker hastigheten i dette punktet
        l_ap = v_ap[n_t_ap]*dt
        n_t_ap +=1
        
    per_i = r.index(max(r))  #finner tidspunktet når planeten er nærmest solen
    for k in range(per_i-100, per_i+100) :   #beregner for posisjonene inærheten av periapsen, spesifikt dt*100 på hver side
        A_per.append( 0.5*np.linalg.norm( np.cross((x[k], y[k], 0), (x[k+1], y[k+1], 0)) ) ) #tar kryssproduktet for å regne arealet til trekanten ved periapsen
        v_per.append(np.linalg.norm([ v_array[0][k], v_array[1][k] ]))  #sjekker hastigheten i dette punktet
        l_per += v_per[n_t_per]*dt
        n_t_per +=1
    
    #finner snitthastigheten
    v_ap = sum(v_ap)/n_t_ap      
    v_per = sum(v_per)/n_t_per
    
    #finner snittratioen av A ved apoapsen og periapsen
    dA = np.sum( np.array(A_ap) ) / np.sum(np.array(A_per)) 
    
    return(v_ap, v_per, A_ap, A_per, dA, l_ap, l_per)


#kjører analytisk-plot for å sammenligne med numerisk
f = np.linspace(0, 2*np.pi, 100000)     
r = np.array(list(np.zeros(len(f)) for i in range(len(system.radii))))
e = system.eccentricities
a = system.semi_major_axes
p = a*(1-e**2)
x_analytisk =  []
y_analytisk = []
print(system.radii)
for i in range(len(system.radii)):
    f_i = f + (np.pi-system.aphelion_angles[i])
    r[i] = (p[i]/(1+e[i]*np.cos(f_i)))

#for i in range(len(system.radii)):
#    r[i] += system.aphelion_angles[i]

for i in range(len(system.radii)):
    x_analytisk.append( list(r[i]*np.cos(f)) )
    y_analytisk.append( list(r[i]*np.sin(f)) )
#print(x_analytisk)
#print(y_analytisk)

#plt.axes(projection="polar")
#for j in range(len(r)):
#    plt.polar(f, r[j], label = (f"planet {j}"))
#plt.show()



results = []
kepler_results = []

for i in range(len(system.radii)):
    results.append(solve(i))

p_analytisk = []
for i in range(len(system.radii)):  #finner analytisk omløpstid ved å se på når avstanden fra planeten til solen er størst
    x = np.array(results[i][0][0])
    y = np.array(results[i][0][1])
    r = np.sqrt(x**2 + y**2)     #konverterer x og y posisjonene til avstand
    p_analytisk.append( sc.signal.find_peaks(r) )
print(p_analytisk)

deviation_yr1 = []
deviation_tend = []
for i in range(len(system.radii)):      #Sjekker for forskjeller i periodene mellom Kepler og Newton
    print(system.semi_major_axes)
    deviation_yr1.append( abs((p_analytisk[i][0][0]*dt-p_analytisk[i][0][0]*dt)**2 - system.semi_major_axes[i]**3)/system.semi_major_axes[i]**3 )
    deviation_tend.append( abs((p_analytisk[i][0][0]*dt-p_analytisk[i][0][0]*dt)**2 - system.semi_major_axes[i]**3)/system.semi_major_axes[i]**3 )

for k in range(len(system.radii)):  #Løser for arealet sveipet ut av planetene
    kepler_results.append( kepler(results[k][0], results[k][1], k) )

labels = ['planet 1','planet 2','planet 3','planet 4','planet 5','planet 6','planet 7', 'planet 8' ]
for j in range(len(system.radii)):
    plt.plot(results[j][0][0], results[j][0][1], label=labels[j])
    plt.plot(x_analytisk[j], y_analytisk[j], color = "b", linestyle="dashed")

print(len(system.radii))
plt.title("Planetenes baner, analytisk og numerisk")
plt.xlabel('x-posisjon i AU')
plt.ylabel('y-posisjon i AU')
plt.axis('square')
plt.grid(True)
plt.legend()
plt.show()

#for h in range(len(results)):
#    plt.plot( np.linspace(0, t_end, timesteps), results[h][2])

for l in range(len(results)):
#    plt.plot(np.linspace(0 , 1000/dt, 1000) , np.array(kepler_results[l][3][0:1000:1])-np.array(kepler_results[l][2][0:1000:1]))
    print(f'forskjell i omløpstid fra forventet i prosent for planet{l} er ved år 1 {deviation_yr1[l]} og ved slutten {deviation_tend[l]}')
    print(f'ratio på arealet = {kepler_results[l][4]}')
    print(f'planet med index {l} har gjennomsnittlig fart på {kepler_results[l][0]}AU/yr ved apoapsen og {kepler_results[l][1]}AU/yr ved periapsen')
    print(f'planet med index {l} dro {kepler_results[l][5]}AU for å spanne ut arealet ved apoapsen, og {kepler_results[l][5]}AU for samme arealet ved periapsen')
#plt.show()






