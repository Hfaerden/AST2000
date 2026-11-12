#Uten kodemal :)
import numpy as np
from ast2000tools.space_mission import SpaceMission as SM
from ast2000tools.solar_system import SolarSystem as SS
import ast2000tools.utils as utils

seed = utils.get_seed('natanies')

mission = SM(seed)  #setter opp mission-raketten
system = SS(seed)   #setter opp solsystemet vårt

masses = system.masses #henter massene til planetene våre 
starmass = system.star_mass #henter massene til stjernen vår 

#Henter planetposisjonene våre fra fil og lagrer de i variabel navn
with np.load("orbits.npz") as f: 
    times = f["times"]
    positions = f["planet_positions"]

def interpolation(time):
    upper = np.searchsorted(times, time) #Finner nærmeste indeks til input tiden, rundet opp, som en upper bound
    lower = upper - 1 #Tar forrige indeks som lower bound 
    t2, t1 = times[upper], times[lower] #Finner tidene ved indeksene til bounds 
    p1 = positions[:,:,lower] #Henter posisjonene ved tidene til bounds 
    p2 = positions[:,:,upper] 
        
    interpolatedpositions = (p1 + (p2-p1)/(t2-t1) * (time-t1)) #Regner ut lineær interpolerte posisjoner
    return(interpolatedpositions) #Returnerer posisjonene

def webemovin(starttime, startpos, startv, timestepsN, dt):
    #Definerer arrays for å lagre variabler
    x = np.zeros((timestepsN+1, 2))
    v = np.zeros((timestepsN+1, 2))
    a = np.zeros((timestepsN+1, 2))

    #Definerer startposisjon og hastighet i første element i arrays
    x[0] = np.array([startpos])
    v[0] = np.array([startv])

    #Finner aksellerasjon ved starten for x of y komponenter for hver planet 
    initpos = interpolation(starttime)
    r0 = np.linalg.norm(initpos)
    a0comp = [(-initpos[0]/r0) * (2*np.pi)**2 * masses/(r0**2), (-initpos[1]/r0) * (2*np.pi)**2 * masses/(r0**2)]
    a0 = [np.sum(a0comp[0]), np.sum(a0comp[1])]

    #Setter aksellerasjons startverdien i array 
    a[0] = np.array(a0)

    #definerer tidsarray
    t = np.zeros(timestepsN+1)
    t[0] = starttime

    for i in range(timestepsN):
        #Henter posisjonene ved tiden
        positions = interpolation(t[i])
        #Regner ut r vektoren fra oss til hver planet 
        r = np.sqrt((positions[0]-x[i][0])**2 + (positions[1]-x[i][1])**2)

        #Regner ut aksellerasjonen med posisjonene og r vektorer
        acomponents = [(-(positions[0]-x[i][0])/r) * (2*np.pi)**2 * masses/(r**2), -((positions[1]-x[i][1])/r) * (2*np.pi)**2 * masses/(r**2)]
        #Finner aksellerasjonen mot solen 
        rsun = np.sqrt(x[i][0]**2 + x[i][1]**2)
        asun = [(-x[i][0]/rsun) * (2*np.pi)**2 * starmass / (rsun**2), (-x[i][1]/rsun) * (2*np.pi)**2 * starmass / (rsun**2)]
        #Regner ut total aksellerasjon
        totalac = np.array([np.sum(acomponents[0]), np.sum(acomponents[1])]) + asun

        #Oppdaterer verdier etter hverandre
        a[i+1] = totalac
        v[i+1] = (a[i+1] * dt) + v[i]
        x[i+1] = (v[i+1] * dt) + x[i]
        t[i+1] = t[i] + dt

    return((starttime + timestepsN * dt), x[-1], v[-1])

print(webemovin(0.3, [-0.2438, 0.1800], [-6.99, -7.39], 100_000, 1/(100_000)))
