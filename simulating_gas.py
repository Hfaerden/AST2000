import numpy as np 
import time
import math
import ast2000tools.constants as const

K = const.k_B

class Gassimulation:
    def __init__(self, N, T, m, totaltime, dt, len):
        #Denne tar in paramterene vi kan redigere for å justere motoren vår
        self.sd = np.sqrt(K*T/m) #Regner ut standardavvik for hastighet
        self.N = N
        self.len = len
        self.m = m
        self.dt = dt
        self.T = T
        self.Npos = np.random.uniform(-self.len/2, self.len/2, (self.N, 3)) #Genererer posisjoner med uniform fordeling
        self.Nvelocity = np.random.normal(0, self.sd, (self.N, 3)) #Genererer hastighet i hver dimensjon med standardavvik fra tidligere og gjennomsnitt = 0
        self.momentumz = np.zeros(self.N, dtype = float) #Lager tomme lister for senere bruk 
        self.escaped = np.zeros(self.N)
        self.timesteps = int(totaltime / dt) + 1 #legger til 1, som blir fjerna i loopen under, dette er for å lage en fin output lettere 
        self.totaltime = totaltime
    def runsim(self):
        #Denne funksjonen kjører simulasjonen på de initialiserte verdiene fra over
        timeinit = time.time()
        for i in range(1, self.timesteps):
            self.Npos = self.Npos + self.Nvelocity * self.dt #Oppdaterer posisjonene, ved å legge til hastighet ganger dt 
            self.Noutside = abs(self.Npos) > self.len/2 #Setter en array med lengde N, hvor vært element har 3 verdier basert på om den er utenfor boksen i x, y, og z dimensjonene
            self.Nvelocity = self.Nvelocity - 2*self.Nvelocity * self.Noutside #Dersom Noutside = false endrer denne ingen av tallene, dersom Noutside = true får vi Nvelocity = -Nvelocity for den partikellen i den dimensjonen
            self.Nsenteredx = abs(self.Npos[:,0]) < 0.25*self.len #denne og den under sjekker om x og y verdiene er innenfor åpningen
            self.Nsenteredy = abs(self.Npos[:,1]) < 0.25*self.len
            self.Centerednoth = np.logical_and(self.Nsenteredx, self.Nsenteredy)
            self.Escape = np.logical_and((self.Centerednoth), (self.Npos[:,2] < -self.len/2)) #denne sjekker om den er utenfor boksen akkurat i hullet
            self.momentumz += self.Nvelocity[:,2] * self.m * self.Escape #Legger til momentum fra unslippende partikeller
            self.escaped += self.Escape #Teller hvor mange partikeller har unsluppet. 
            if i%20 == 0: #lager en fancy progress bar for simulasjonen
                print("#"*int((i/20)) + "-"*(50-int(i/20)))
        Forcez = (sum(self.momentumz)/self.totaltime) #regner ut total kraft 
        print(f"Kraften i z retning er {Forcez}") #Denne og under printer relevant informasjon
        print("Simulasjon ferdig på " + str(round(time.time()-timeinit, 2)) + " sekunder \n")
        print(f"momentumsum fra unslippende partikler: {sum(self.momentumz)}")
        print("Partikler som slapp ut: " + str(round(sum(self.escaped))))
        print(f"gjennomsnittlig hastighet for unsluppne partikler i z-retning er {(sum(self.momentumz)/self.m)/round(sum(self.escaped))}")
#Lager en klasse med paramterene oppgitt i oppgaven og kjører med den. Hvis filen er importert kjører ikke det under, og vi kan heller lage en annen boks med andre parametere. 
if __name__ == "__main__":
    engine = Gassimulation(10**5, 3.5*10**3, const.m_H2, 10**(-9), 10**(-12), 10**(-7))
    engine.runsim()
    Vsquared = (engine.Nvelocity)**2 
    print(f"Average velocity = {np.average(engine.Nvelocity)}")
    print(f"Average velocity squared = {np.average(Vsquared)}")
    print(f"Average velocity squared expected = {3 * K*engine.T/engine.m}")
