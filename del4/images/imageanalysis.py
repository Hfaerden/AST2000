import numpy as np
import matplotlib.pyplot as plt
from ast2000tools.space_mission import SpaceMission as SM
from ast2000tools.solar_system import SolarSystem as SS
import ast2000tools.constants as const
import ast2000tools.utils as utils
from PIL import Image 


seed = utils.get_seed('natanies')
mission = SM(seed)  #setter opp mission-raketten
system = SS(seed)   #setter opp solsystemet vårt

#Oppgave 1

img = Image.open("sample0000.png") #Åpner bildet
width, height = img.size #Henter størrelsen
print(str(width) + "\n" + str(height)) #Printer ut størrelsen

#Oppgave 2:
print(f"x range = {width/2} \n y range = {height/2}") #Printer ut rekkevidden på X og Y rundt origo i midten av bildet

#Oppgave 3:
#Funksjon som regner ut theta og phi verdier for en gitt x og y verdi, og gitt en theta0 og phi0
def coords(xcor, ycor, theta0, phi0):
    p = np.sqrt(xcor**2+ycor**2)
    B = 2 * np.arctan(p/2)
    theta = theta0 - np.arcsin(np.cos(B) * np.cos(theta0) + (ycor/p) * np.sin(B) * np.sin(theta0))
    phi = phi0 + np.arctan((xcor * np.sin(B))/(p*np.sin(theta0) * np.cos(B) - ycor * np.cos(theta0) * np.sin(B)))
    return(theta, phi)

fov = 70 * np.pi/180 #Setter field of view i radianer

#Regner ut max og min verdier
xmax = 2*np.sin(fov/2)/(1+np.cos(fov/2))
ymax = 2*np.sin(fov/2)/(1+np.cos(fov/2))
xmin = -1 * xmax 
ymin = -1 * ymax

x = np.linspace(xmin, xmax, width)
y = np.linspace(ymin, ymax, height)

X, Y = np.meshgrid(x, y) #Setter opp meshgrid for X, Y

theta, phi = coords(X, Y, np.pi/2, 0) #Bruker X, Y for å lage et theta, phi meshgrid med funksjonen definert over

#Oppgave 4:
img = Image.new(mode = "RGB", size = (len(theta[0]), len(theta))) #Lager et bildet.
#  Merk at theta og phi har samme dimensjoner, så kun en brukes. innerste liste blir brukt først grunnet motsatt X, Y for matriser og bilder

file = np.load("himmelkule.npy") #Henter himmmelkule filen
for i, j in enumerate (phi): #Går over antall y verdier i listen 
    for k, l in enumerate(j): #Går over antall x verdier i listen 
        index = mission.get_sky_image_pixel(theta[i][k], phi[i][k]) #Henter indeks for RGB verdien for gitte kordinat
        RGB = file[index] #Henter RGB verdien fra filen
        img.putpixel((k, (height - 1) - i), tuple(int(x) for x in (RGB[2:]))) #Setter RGB verdien i bildet. K brukt først siden x og y er reversert. minuser høyden siden PIL begynner øverst

img.save("output.png") #Largrer bildet


def generalcase(degree):
    theta, phi = coords(X, Y, np.pi/2, degree * np.pi / 180) #Bruker X, Y for å lage et theta, phi meshgrid med funksjonen definert over
    for i, j in enumerate (phi): #Går over antall y verdier i listen 
        for k, l in enumerate(j): #Går over antall x verdier i listen 
            index = mission.get_sky_image_pixel(theta[i][k], phi[i][k]) #Henter indeks for RGB verdien for gitte kordinat
            RGB = file[index] #Henter RGB verdien fra filen
            img.putpixel((k, (height - 1) - i), tuple(int(x) for x in (RGB[2:]))) #Setter RGB verdien i bildet. K brukt først siden x og y er reversert. minuser høyden siden PIL begynner øverst

    img.save(f"output {degree}.png")

input("ADVARSEL FOR MEG SELV / NATANIEL I FREMTIDEN: Denne filen kommer til å legge 360 bilder i mappen du er i. Klikk enter for å fortsette, ctrl C for å avbryte \n -Hennie :)")
for i in range(360):
    generalcase(i)
