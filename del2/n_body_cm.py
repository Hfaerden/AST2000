# UTEN KODEMAL!!

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

k = [3, 4, 5] #planet nr 3, 4, og 5
pos = system.initial_positions #Henter posisjonen ved t = 0
vel = system.initial_velocities #Henter hastigheten ved t = 0
starmass = system.star_mass #Henter solmassen vår 
planetmass = []
p_pos = []
p_posx, p_posy = [], []

for i in range(len(k)):
    planetmass.append(system.masses[k[i]])
    p_pos.append( [system.initial_positions[0][k[i]], system.initial_positions[1][k[i]]] )  #henter planetposisjonene
    p_posx.append( system.initial_positions[0][k[i]]) #bruker disse senere for å få rett datastruktur
    p_posy.append( system.initial_positions[1][k[i]])

p_pos = np.array(p_pos)

yearlen = system.semi_major_axes[0]**(3/2) #Regner ut lengden på et år i jordår, ved Kepler's tredje
#print(yearlen)
timesteps_per_year = 10000 #Setter hvor mange tidssteg per år 
dt = yearlen/(timesteps_per_year) #Regner ut tidsintervalet vi skal ved dt og årlengden
#print(d++t)
t_end = 100
timesteps = round(t_end * timesteps_per_year) #Regner ut hvor mange tidssteg vi trenger totalt for t år
M=0

for i in range(len(k)):
    M +=float(planetmass[i])

M += float(starmass)
mu = (planetmass[0]*planetmass[1]*planetmass[2]*starmass)/(M)
G = 4*np.pi**2

#print(p_pos)
CM = np.array([(planetmass[0]*p_pos[0] + planetmass[1]*p_pos[1] + planetmass[2]*p_pos[2]) / M ])    #Finner posisjon og hastighet til CM
CMvx = float( (planetmass[0]*vel[0][0]+planetmass[1]*vel[0][1]+planetmass[2]*vel[0][2]) / M)
CMvy = float( (planetmass[0]*vel[1][0]+planetmass[1]*vel[1][1]+planetmass[2]*vel[1][2]) / M)

v_s = [float(-CMvx), float(-CMvy)]      #Konverterer til CM-ref systemet
v_p = [ [[],[],[]], [[],[],[]] ]
for i in range(len(k)):
    v_p[0][i].append( float(vel[0][k[i]]-CMvx) )
    v_p[1][i].append( float(vel[1][k[i]]-CMvy) )

#print(v_p)

CM = CM[0]
#print(CM)
s_pos = -CM
#print(CMvx, CMvy)
#print(v_s)


#print(CM)



#@jit
def solve(pos_px, pos_py, pos_sx, pos_sy, v_px, v_py, v_sx, v_sy):
    x1_p, y1_p = pos_px, pos_py    #Henter posisjonene ved t = 0 som en liste
    x1_s, y1_s = pos_sx, pos_sy
    ax_p, ay_p = [[],[],[]], [[],[],[]]     #setter opp nøstede lister for akselerasjonene våre
    r, rx_p2s, ry_p2s = 0, 0, 0
    asx1, asy1 = 0, 0
    
    
    for i in range(len(k)):
        r = float(np.sqrt((x1_p[i][0]-x1_s)**2 + (y1_p[i][0]-y1_s)**2))    #relative posisjonen 
        ax_p[i].append(float( -((x1_p[i][0]-x1_s) * 4*np.pi**2 * starmass / (r**3) ) ))
        print(ax_p)
        ay_p[i].append(float( -((y1_p[i][0]-y1_s) * 4*np.pi**2 * starmass / (r**3) ) ))    #Regner ut aksellerasjonen ved t = 0
        
        asx1 += ( -(-rx_p2s * 4*np.pi**2 * planetmass[i] / (r**3) ) )
        asy1 += ( -(-ry_p2s * 4*np.pi**2 * planetmass[i] / (r**3) ) )
    
    ax_s, ay_s = [asx1], [asy1]
    x_p, y_p = x1_p, y1_p   #Lager arrays for poisjonene våre  
    x_s, y_s = [x1_s], [y1_s]
    vx_p, vy_p = v_px, v_py  #lager arrays for hastighetene 
    #print(vx_p)
    vx_s, vy_s = [v_sx], [v_sy] 
    #U = G*M*mu/r
    #K = 0.5*mu*(x_p[0]*vx_p[0] + y_p[0]*vy_p[0] - x_s[0]*vx_s[0] + y_s[0]*vy_s[0])
    #(x_p[0]*vx_p[0] + y_p[0]*vy_p[0] - x_s[0]*vx_s[0] + y_s[0]*vy_s[0])
    #np.dot([x_p[0], y_p[0]], [vx_p[0], vy_p[0]]) - np.dot([x_s[0], y_s[0]] , [vx_s[0], vy_s[0]])
    #E = [K-U]
    
     
    for i in range(timesteps - 1): #Euler-cromer loop for å regne ut akselerasjon, posisjon, og hastighet
    
        asx = 0
        asy = 0

        for j in range(len(k)):     #løper gjennom for hver planet i simulasjonen
            #print(x_p)
            r = float(np.sqrt((x_p[j][i]-x_s[i])**2 + (y_p[j][i]-y_s[i])**2))
            rx_p2s = (x_p[j][i]-x_s[i])
            ry_p2s = (y_p[j][i]-y_s[i])
            #print(rx_p2s, ry_p2s)
            #break
            #print(r)
            ax_p[j].append(float( -( rx_p2s * 4*np.pi**2 * starmass / (r**3)) ))
            ay_p[j].append(float( -( ry_p2s * 4*np.pi**2 * starmass / (r**3)) ))
            vx_p[j].append(vx_p[j][i]+ax_p[j][i+1]*dt)
            vy_p[j].append(vy_p[j][i]+ay_p[j][i+1]*dt)
            x_p[j].append(x_p[j][i]+vx_p[j][i+1]*dt)
            y_p[j].append(y_p[j][i]+vy_p[j][i+1]*dt)
            #print(ax_p, ay_p)
            
            asx += ( -(-rx_p2s * 4*np.pi**2  * planetmass[j] / (r**3) ) )
            asy += ( -(-ry_p2s * 4*np.pi**2  * planetmass[j] / (r**3) ) )
            #print(asx, asy)
            
        
        ax_s.append(asx)
        ay_s.append(asy)
        vx_s.append(vx_s[i]+ax_s[i+1]*dt)
        vy_s.append(vy_s[i]+ay_s[i+1]*dt)
        x_s.append(x_s[i]+vx_s[i+1]*dt)
        y_s.append(y_s[i]+vy_s[i+1]*dt)
        #break
        #U = G*M*mu/r
        #K = 0.5*mu*(x_p[0]*vx_p[0] + y_p[0]*vy_p[0] - x_s[0]*vx_s[0] + y_s[0]*vy_s[0])
        #E.append(K-U)
        
        
        
    return([[x_p, y_p], [vx_p, vy_p]], [[x_s, y_s], [vx_s, vy_s]])#, E) #Returnerer x og y posisjonene , samt vx og vy, og totalenergien


print(v_p)
print(v_p[1])
px_format, py_format = [ [p_posx[0]], [p_posx[1]], [p_posx[2]] ], [[p_posy[0]], [p_posy[1]], [p_posy[2]]]
p_pos = [ px_format , py_format]
vx_format, vy_format = 1, 1

#print(p_pos)
results = solve(p_pos[0], p_pos[1], s_pos[0], s_pos[1], v_p[0], v_p[1],  v_s[0], v_s[1])

for i in range(len(k)):
    plt.plot(results[0][0][0][i], results[0][0][1][i], label='planet')
    
plt.plot(results[1][0][0], results[1][0][1], label='Stjernen')
plt.title('n-legemet systemet, med planet nr 4, 5 og 6')
plt.xlabel('x-posisjon i AU')
plt.ylabel('y-posisjon i AU')
plt.axis('square')
plt.grid(True)
plt.legend()
plt.show()


'''E = np.array(results[-1])
E_mean = np.mean(E)
E_deviation = (abs(E/E_mean - 1))*100
peaks = sc.signal.find_peaks(E)
peak_diff = abs(peaks[0][0]-peaks[0][-1])/(peaks[0][0]+peaks[0][-1])
print(peak_diff)

plt.plot(np.linspace(0, t_end, timesteps), E_deviation )
plt.show()

vel_curve = np.array(results[1][1][0])
vel_curve += np.random.normal(0, 0.2*np.max(vel_curve), len(vel_curve))    #legger til støy lik 1/5 av den største verdien
peculiar_vel = 0.1
vel_curve += peculiar_vel #legger til hastigheten til CM som sett fra observatøren
'''
'''
plt.plot( np.linspace(0, t_end, timesteps), vel_curve )
plt.title('radiell hastighet')
plt.xlabel('tid i år')
plt.ylabel('hastighet i AU/yr')
plt.grid(True)
plt.show()
'''
