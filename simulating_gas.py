import numpy as np 
import time
import math

K = 1.38*10**(-23)

class Gassimulation:
    def __init__(self, number, temperature, mass, totaltime, dt, length):
        self.sd = np.sqrt(K*temperature/mass)
        self.N = number 
        self.len = length
        self.m = mass
        self.dt = dt
        self.Npos = np.random.uniform(-self.len/2, self.len/2, (self.N, 3))
        self.Nvelocity = np.random.normal(0, self.sd, (self.N, 3))
        self.momentumz = np.zeros(self.N, dtype = float)
        self.escaped = np.zeros(self.N)
        self.timesteps = int(totaltime / dt) + 1
        self.totaltime = totaltime
    def runsim(self):
        timeinit = time.time()
        for i in range(1, self.timesteps):
            
            self.Npos = self.Npos + self.Nvelocity * self.dt 
            self.Noutside = abs(self.Npos) > self.len/2
            self.Nvelocity = self.Nvelocity - 2*self.Nvelocity * self.Noutside
            
            self.Nsenteredx = abs(self.Npos[:,0]) < 0.25*self.len
            self.Nsenteredy = abs(self.Npos[:,1]) < 0.25*self.len
            self.Centerednoth = np.logical_and(self.Nsenteredx, self.Nsenteredy)
            self.Escape = np.logical_and((self.Centerednoth), (self.Npos[:,2] < -self.len/2)) 
            
            self.momentumz += self.Nvelocity[:,2] * self.m * self.Escape
            self.escaped += self.Escape
            #collisioncount += Noutside
            if i%20 == 0:
                print("#"*int((i/20)) + "-"*(50-int(i/20)))
        Forcez = (sum(self.momentumz)/self.totaltime)
        print(f"Kraften i z retning er {Forcez}")
        print("Simulasjon ferdig på " + str(round(time.time()-timeinit, 2)) + " sekunder \n")
        print(f"momentumsum fra unslippende partikler: {sum(self.momentumz)}")
        print("Partikler som slapp ut: " + str(round(sum(self.escaped))))
        print(f"gjennomsnittlig hastighet for unsluppne partikler i z-retning er {(sum(self.momentumz)/self.m)/round(sum(self.escaped))}")

if __name__ == "__main__":
    engine = Gassimulation(10**5, 3*10**3, 3.35*10**(-24), 10**(-9), 10**(-12), 10**(-6))
    engine.runsim()
