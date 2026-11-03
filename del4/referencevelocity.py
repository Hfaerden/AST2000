import numpy as np
import matplotlib.pyplot as plt
from ast2000tools.space_mission import SpaceMission as SM
from ast2000tools.solar_system import SolarSystem as SS
import ast2000tools.constants as const
import ast2000tools.utils as utils


seed = utils.get_seed('natanies')

mission = SM(seed)  #setter opp mission-raketten
system = SS(seed)   #setter opp solsystemet vårt

Hspectral = 656.3
c = const.c

angles = mission.star_direction_angles
dopplershifts = mission.star_doppler_shifts_at_sun

#Regner ut radiell hastighet til noe gitt en dopllershift 
def velocityfind(dopplershift):
    fractionsquared = ((dopplershift/Hspectral)+1)**2
    return(c * (fractionsquared-1)/(fractionsquared+1))

#Finner hastigheten av stjernen vår i forhold til refereansestjernene
print("Radial velocity for our sun: \n")
for i in (dopplershifts):
    print(f"phi1: {velocityfind(i)}")


def referencevelocity(lam1, lam2, phi1, phi2):
    v1 = velocityfind(lam1)
    v2 = velocityfind(lam2)
    matrix = np.array([np.sin(phi2), -np.sin(phi1)], [-np.cos(phi2), np.cos(phi1)])
    return(1/(np.sin(phi2-phi1)) * (matrix @ [v1,v2]))

def referencevelocitytosun(lam1, lam2):
    #Finner vår hastighet i x, y i forhold til stjernene
    ourrefv = referencevelocity(lam1, lam2, angles[0], angles[1])
    #Finner solen sin hastighet i x, y i forhold til stjernene
    sunrefv = referencevelocity(dopplershifts[0], dopplershifts[1], angles[0], angles[1])
    return(ourrefv - sunrefv)


