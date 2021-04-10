import matplotlib.pyplot as plt
from math import *
import numpy as np
from scipy import fft
from scipy import signal
plt.style.use("IceCube")


upperlim = 4999
lowerlim = 0
N = upperlim - lowerlim + 1

def getData(index):
    f = open("./9e10_avg_test/test_data/test_" + str(index) + ".txt", 'r')

    count = 0
    z_list = []
    dE_List = []
    r_list = []

    while 1:
        temp = f.readline().split()
        if len(temp) <= 1:
            break

        if count >= lowerlim and count <= upperlim:

            z = float(temp[0]) / 1.6
            dE = float(temp[1]) / (1.85775164151e-05 / 3)

            z_list.append(float(temp[0]))
            dE_List.append(float(temp[1]))

            r_list.append(sqrt(z**2 + dE**2))

        count += 1

    return np.array(z_list) - np.array(z_list).mean(), np.array(dE_List), np.array(r_list) - np.array(r_list).mean()


def singleFFT(r_list, r_or_z):
    if len(r_list) < N:
        if r_or_z == 'r':
            return -1, -1
        else:
            return -1

    windows = signal.hanning(len(r_list))
    r_list = r_list * windows

    outr = fft(r_list)
    outr = outr[1:len(outr) / 2 + 1]
    norm_listr = []

    for i in range(len(outr)):
        norm_listr.append(np.absolute(outr[i]))




    if r_or_z == 'r':


        f = plt.figure()

        plt.plot(np.linspace(0, 0.5, N / 2), norm_listr)

        plt.xlabel("tune")

        plt.show()



        tune1 = float(np.argsort(norm_listr[int((0.005/0.5) * N/2) : int((0.01/0.5) * N/2)])[-1] + int((0.005/0.5) * N/2)) / N
        tune2 = float(np.argsort(norm_listr[int((0.01/0.5) * N/2) : int((0.1/0.5) * N/2)])[-1] + int((0.01/0.5) * N/2)) / N

        #print ((0.005 / 0.5) * N/2) / N
        #print ((0.01 / 0.5) * N/2) / N


        return tune1, tune2


    else:

        """
        f = plt.figure()

        plt.plot(np.linspace(0, 0.5, N / 2), norm_listr)

        plt.xlabel("tune")

        plt.show()
        """

        tune_list = float(np.argsort(norm_listr[1:])[-1]) / N

        return tune_list



FFT_z_list = []
FFT_r_list = []
FFT_r2_list = []

for index in range(101)[1:]:
    z_list, dE_List, r_list = getData(index)

    tuner1, tuner2 = singleFFT(r_list, 'r')

    #print(round(tuner1, 5), round(tuner2, 5))

    tunez = singleFFT(z_list, 'z')

    if tunez != -1:
        FFT_z_list.append(tunez)
        FFT_r_list.append(tuner1)
        FFT_r2_list.append(tuner2)


f1 = plt.figure()
plt.hist(FFT_z_list, bins=50, label='z tune, std='+str(round(np.std(FFT_z_list), 5)) + ", avg="+str(round(np.average(FFT_z_list), 5)))
plt.legend()

f2 = plt.figure()
plt.hist(FFT_r_list, bins=50, label='tune of first r peak, std='+str(round(np.std(FFT_r_list), 5)) + ", avg="+str(round(np.average(FFT_r_list), 5)))
plt.legend()


f3 = plt.figure()
plt.hist(FFT_r2_list, bins=50, label='tune of second r peak, std='+str(round(np.std(FFT_r2_list), 5)) + ", avg="+str(round(np.average(FFT_r2_list), 5)))
plt.legend()
plt.show()

exit()


f = plt.figure()
plt.scatter(z_list/1.6, dE_List/(1.85775164151e-05 / 3))
plt.xlabel("z in sigma_z")
plt.ylabel("dE in 0.62e-5GeV")


f2 = plt.figure()
plt.plot(range(len(r_list)), r_list)
plt.xlabel("turns")
plt.ylabel("r")
plt.show()