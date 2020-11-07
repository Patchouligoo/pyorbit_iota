import matplotlib.pyplot as plt
from math import *
import numpy as np
from scipy import fft
from scipy import signal
plt.style.use("IceCube")

upperlim = 1499
lowerlim = 1000
N = upperlim - lowerlim + 1

def getData(index):
    f = open("./Files_shifted/test_particles_intensity_10/" + str(index) + ".txt", 'r')

    count = 0
    x = []
    y = []


    while 1:
        temp = f.readline().split()
        if len(temp) <= 1:
            break

        if count >= lowerlim and count <= upperlim:
            x.append(float(temp[0]))
            y.append(float(temp[1]))

        count += 1

    return np.array(x), np.array(y)


def singleFFT(index):

    x, y = getData(index)

    if len(x) < N:
        #print("??????")
        return -1, -1

    windows = signal.hanning(len(x))
    x = x * windows
    y = y * windows

    outx = fft(x)
    outx = outx[1:len(outx)/2 + 1]
    norm_listx = []

    for i in range(len(outx)):
        norm_listx.append(np.absolute(outx[i]))

    """
    if index == 1502:
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
    if index == 1502:
        f1 = plt.figure()
        plt.plot(np.linspace(0, 0.5, len(norm_listy)), norm_listy)
        plt.xlabel("tune")
        plt.ylabel("amplitude")
        plt.title("y axis")
        plt.show()
    """

    tune_y = float(np.argmax(norm_listy))/N

    #if index == 1502:
    #    print(tune_x, tune_y)
    #    exit()

    return tune_x, tune_y






file_reader = open("ftprint_SC.txt", 'r')
qx_list = []
qy_list = []
while 1:
    temp = file_reader.readline().split()

    if len(temp) < 1:
        break

    qx = (5.31 + float(temp[2]))%1
    qy = (5.32 + float(temp[3]))%1

    if qx > 0.5:
        qx = 1 - qx
    if qy > 0.5:
        qy = 1 - qy

    qx_list.append(qx)
    qy_list.append(qy)



Qx = []
Qy = []

for i in range(3000):
    i = i + 1
    qx, qy = singleFFT(i)

    if qx != -1:
        Qx.append(qx)
        Qy.append(qy)

    print("doing linear " + str(i) + " " + str(qx) + " " + str(qy))



f = plt.figure()
plt.scatter(Qx, Qy, label='FFT result')
#plt.scatter(qx_list, qy_list, label='expected')
plt.xlabel("Qx")
plt.ylabel("Qy")

plt.xlim((0, 0.5))
plt.ylim((0, 0.5))
plt.legend()
plt.show()







