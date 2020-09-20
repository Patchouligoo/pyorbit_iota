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

def getAvgMin(x_list, y_list):

    dis_list = []

    for i in range(len(x_list)):
        dis_list.append(np.sqrt(pow(x_list[i], 2) + pow(y_list[i], 2)))


    dis_list = np.array(dis_list)

    min_ape = dis_list.min()

    return round(dis_list.sum()/len(dis_list), 2), round(min_ape, 2)


def getData(dim, tar):
    f = open("./" + str(dim) +"D_final_" + str(tar) +".dat", 'r')

    for i in range(15):
        f.readline()

    survivor_list = []

    while 1:
        temp = f.readline().split()

        if len(temp) < 1:
            break


        survivor_list.append(int(temp[6]))

    f.close()

    p = open("./" + str(dim) +"D_initial_" + str(tar) + ".dat", 'r')

    x_list = []
    y_list = []

    pre_x = 0
    pre_y = 0
    pre_id = -1


    for i in range(15):
        p.readline()


    while 1:
        temp = p.readline().split()

        if len(temp) < 1:
            break

        id = int(temp[6])

        if id in survivor_list:

            if id != pre_id + 1:
                x_list.append(pre_x)
                y_list.append(pre_y)

            pre_id = id
            pre_x = float(temp[0])
            pre_y = float(temp[2])


    x_list = np.array(x_list)
    y_list = np.array(y_list)

    return x_list/sx, y_list/sy



f = plt.figure()
ax = f.add_subplot(111)

""" ---- 4D ----"""
x_list, y_list = getData(4, "bunch_0p")
avg, min = getAvgMin(x_list, y_list)
plt.scatter(x_list, y_list, color='black')
plt.plot(x_list, y_list, label='py 4D, min: ' + str(min) + " avg: " + str(avg), color='black', linewidth=3)

""" ---- 6D----"""
x_list, y_list = getData(6, "bunch_1p")
avg, min = getAvgMin(x_list, y_list)
plt.scatter(x_list, y_list, color='blue')
plt.plot(x_list, y_list, label='py 6D, min: ' + str(min) + " avg: " + str(avg), color='blue', linewidth=3)




"""----- read MADX 4D-------"""
p = open("4D_survivor_fringe_0p.txt", 'r')
x = []
y = []
while 1:
    temp = p.readline().split()
    if len(temp) < 1:
        break
    x.append(float(temp[0][:-1]))
    y.append(float(temp[1]))
avg, min = getAvgMin(x, y)
plt.scatter(x, y, color='orange')
plt.plot(x, y, label='MADX 4D, min: ' + str(min) + " avg: " + str(avg), color='orange', linewidth='3')


"""----- read MADX 6D -------"""
p = open("6D_survivor_fringe_1p.txt", 'r')
x = []
y = []
while 1:
    temp = p.readline().split()
    if len(temp) < 1:
        break
    x.append(float(temp[0][:-1]))
    y.append(float(temp[1]))
avg, min = getAvgMin(x, y)
plt.scatter(x, y, color='lime')
plt.plot(x, y, label='MADX 6D, min: ' + str(min) + " avg: " + str(avg), color='lime', linewidth='3')

plt.xlabel("x in sig_x")
plt.ylabel("y in sig_y")
plt.xlim((0, 5))
plt.ylim((0, 5))
plt.legend()
ax.set_aspect('equal', adjustable='box')





f2 = plt.figure()
ax = f2.add_subplot(111)


""" ---- 5D ----"""
x_list, y_list = getData(4, "bunch_1p")
avg, min = getAvgMin(x_list, y_list)
plt.scatter(x_list, y_list, color='red')
plt.plot(x_list, y_list, label='py 5D, min: ' + str(min) + " avg: " + str(avg), color='red', linewidth=3)


""" ---- 6D ----"""
x_list, y_list = getData(6, "bunch_1p")
avg, min = getAvgMin(x_list, y_list)
plt.scatter(x_list, y_list, color='blue')
plt.plot(x_list, y_list, label='py 6D, min: ' + str(min) + " avg: " + str(avg), color='blue', linewidth=3)


"""----- read MADX 5D -------"""
p = open("4D_survivor_fringe_1p.txt", 'r')
x = []
y = []
while 1:
    temp = p.readline().split()
    if len(temp) < 1:
        break
    x.append(float(temp[0][:-1]))
    y.append(float(temp[1]))
plt.scatter(x, y, color='cyan')
avg, min = getAvgMin(x, y)
plt.plot(x, y, label='MADX 5D, min: ' + str(min) + " avg: " + str(avg), color='cyan', linewidth='3')

"""----- read MADX 6D-------"""
p = open("6D_survivor_fringe_1p.txt", 'r')
x = []
y = []
while 1:
    temp = p.readline().split()
    if len(temp) < 1:
        break
    x.append(float(temp[0][:-1]))
    y.append(float(temp[1]))

plt.scatter(x, y, color='lime')
avg, min = getAvgMin(x, y)
plt.plot(x, y, label='MADX 6D, min: ' + str(min) + " avg: " + str(avg), color='lime', linewidth='3')


plt.xlabel("x in sig_x")
plt.ylabel("y in sig_y")
plt.xlim((0, 5))
plt.ylim((0, 5))
plt.legend()
ax.set_aspect('equal', adjustable='box')

plt.show()

