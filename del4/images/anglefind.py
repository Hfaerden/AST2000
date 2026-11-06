#Ingen kodemal :)
from PIL import Image 
import numpy as np

def anglefind(imagename):
    img = Image.open(imagename)
    imgpixel = np.array(img)
    imgpixel = imgpixel.astype(int)
    best = 10**(10)
    bestangle = 0 
    for i in range(360):
        compare = Image.open(f"output {i}.png")
        refarray = np.array(compare)
        refarray = refarray.astype(int)
        d = (np.sum((imgpixel - refarray)**2))
        if d < best:
            best = d
            bestangle = i
    return(bestangle)

print(anglefind("output 127.png"))
#Følgende tester ble også kjørt ved å kjøre generalcase() med vinkelene nedenfor fra imagegeneration
print(anglefind("output 338.6.png"))
print(anglefind("output 90.2.png"))
print(anglefind("output 180.2.png"))
