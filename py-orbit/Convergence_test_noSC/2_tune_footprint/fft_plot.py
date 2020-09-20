import matplotlib.pyplot as plt
from math import *
import numpy as np
from scipy import fft
plt.style.use("IceCube")

upperlim = 4999
lowerlim = 0
N = upperlim - lowerlim + 1

def getData(index):
    f = open("./track_data/" + str(index) + ".txt", 'r')

    count = 0
    turn = []
    x = []
    y = []


    while 1:
        temp = f.readline().split()
        if len(temp) < 1:
            break

        if count >= lowerlim and count <= upperlim:
            turn.append(float(temp[0]))
            x.append(float(temp[1]))
            y.append(float(temp[2]))

        count += 1

    return np.array(turn), np.array(x), np.array(y)


def singleFFT(index):

    turn, x, y = getData(index)

    if len(turn) < N:
        print("??????")
        return -1, -1

    outx = fft(x)
    outx = outx[1:len(outx)/2 + 1]
    norm_listx = []

    for i in range(len(outx)):
        norm_listx.append(np.absolute(outx[i]))

    """
    f = plt.figure()
    plt.plot(np.linspace(0, 0.5, len(norm_listx)), norm_listx)
    plt.xlabel("tune")
    plt.ylabel("amplitude")
    plt.title("x axis")
    """


    tune_x = float(np.argmax(norm_listx))/N

    outy = fft(y)
    outy = outy[1:len(outy) / 2 + 1]
    norm_listy = []

    for i in range(len(outy)):
        norm_listy.append(np.absolute(outy[i]))

    """
    f1 = plt.figure()
    plt.plot(np.linspace(0, 0.5, len(norm_listy)), norm_listy)
    plt.xlabel("tune")
    plt.ylabel("amplitude")
    plt.title("y axis")
    plt.show()
    """

    tune_y = float(np.argmax(norm_listy))/N

    return tune_x, tune_y



Qx = []
Qy = []
for i in range(10000):
    qx, qy = singleFFT(i)

    if qx != -1:
        Qx.append(qx)
        Qy.append(qy)

    print("doing linear " + str(i) + " " + str(qx) + " " + str(qy))


p = open("dynaptune", 'r')
for i in range(8):
    p.readline()
MADQx = []
MADQy = []
while 1:

    temp = p.readline().split()

    if len(temp) < 1:
        break

    if temp[2] != 'nan' and temp[3] != 'nan':
        MADQx.append(float(temp[2]))
        MADQy.append(float(temp[3]))


f = plt.figure()
plt.scatter(Qx, Qy, label='pyORBIT_linear')
plt.scatter(MADQx, MADQy, label='MADX')

plt.xlabel("Qx")
plt.ylabel("Qy")
plt.legend()
plt.show()


exit()
"""------ plot resonance lines ------"""
p = open("tunespace.txt", 'r')
while 1:
    temp = p.readline().split()
    temp2 = p.readline().split()
    temp3 = p.readline().split()

    if len(temp) < 1:
        break

    beginpoint = [float(temp[0]), float(temp2[0])]

    endpoint = [float(temp[1]), float(temp2[1])]


    plt.plot(beginpoint, endpoint, color='black')


plt.xlabel("Qx")
plt.ylabel("Qy")
plt.legend()
plt.show()
