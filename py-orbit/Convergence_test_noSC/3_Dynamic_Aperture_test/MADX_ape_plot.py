import numpy as np
from math import *
import matplotlib.pylab as plt
plt.style.use("IceCube")


alphax = 5.378281542e-09
betax = 0.7905949642
alphay = 5.27404231e-09
betay = 1.163203072
emitlimx = 4.016e-6
emitlimy = 4.016e-6
gx = (1 + alphax**2)/betax
gy = (1 + alphay**2)/betay

sx = np.sqrt(betax * emitlimx)
spx = np.sqrt(gx * emitlimx)
sy = np.sqrt(betay * emitlimy)
spy = np.sqrt(gy * emitlimy)


def getData():

    p = open("./6D_survivor_1p.txt", 'r')

    survivor = []

    while 1:
        temp = p.readline().split()

        if len(temp) < 1:
            break

        if temp[0].isdigit():
            survivor.append(int(temp[0]))

    return survivor

x_list = []
y_list = []

survivor = getData()

m = 0
count = 1
while (m < 50):
    l = 0
    while (l < 100):
        l = l + 1
        x = float(l)/100 * sx * 10 * cos(float(m)/50 * pi/2)
        y = float(l)/100 * sx * 10 * sin(float(m)/50 * pi/2)


        if count in survivor:
            #print(count, x, y)
            x_list.append(x/sx)
            y_list.append(y/sy)

        count += 1


    m = m + 1


f = plt.figure()

plt.scatter(x_list, y_list)

plt.xlabel("x in meter")
plt.ylabel("y in meter")
plt.xlim((0, 5))
plt.ylim((0, 5))
plt.legend()
plt.show()

