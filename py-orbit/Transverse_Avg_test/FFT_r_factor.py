import matplotlib.pyplot as plt
from math import *
import numpy as np
from scipy import fft
from scipy import signal
plt.style.use("IceCube")


upperlim = 4999
lowerlim = 0
N = upperlim - lowerlim + 1

def getData(name, index):
    f = open("./9e10_avg_test" + name + "/test_data/test_" + str(index) + ".txt", 'r')

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


def singleFFT(name, r_list, r_or_z):
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

        """
        f = plt.figure()

        plt.plot(np.linspace(0, 0.5, N / 2), norm_listr)

        plt.xlabel("tune")

        plt.show()
        """

        if name == "_6sigma_trans":
            tune1 = float(np.argsort(norm_listr[int((0.005/0.5) * N/2) : int((0.02/0.5) * N/2)])[-1] + int((0.005/0.5) * N/2)) / N
            tune2 = float(np.argsort(norm_listr[int((0.02/0.5) * N/2) : int((0.1/0.5) * N/2)])[-1] + int((0.02/0.5) * N/2)) / N
        else:
            tune1 = float(np.argsort(norm_listr[int((0.005/0.5) * N/2) : int((0.01/0.5) * N/2)])[-1] + int((0.005/0.5) * N/2)) / N
            tune2 = float(np.argsort(norm_listr[int((0.01/0.5) * N/2) : int((0.1/0.5) * N/2)])[-1] + int((0.01/0.5) * N/2)) / N


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



FFT_z_list_0 = []
FFT_r_list_0 = []
FFT_r2_list_0 = []
err_z_list_0 = []
err_r_list_0 = []
err_r2_list_0 = []

FFT_z_list_3 = []
FFT_r_list_3 = []
FFT_r2_list_3 = []
err_z_list_3 = []
err_r_list_3 = []
err_r2_list_3 = []

FFT_z_list_6 = []
FFT_r_list_6 = []
FFT_r2_list_6 = []
err_z_list_6 = []
err_r_list_6 = []
err_r2_list_6 = []

sig0_loss_ratio = []
sig3_loss_ratio = []
sig6_loss_ratio = []


name = ""  # to 3 sigma_z
for i in range(50):
    print("working on " + name + " " + str(i))
    FFT_temp_z = []
    #FFT_temp_r = []
    #FFT_temp_r2 = []
    for index in range(101)[1:]:
        z_list, dE_List, r_list = getData(name, i*100 + index)

        #tuner1, tuner2 = singleFFT(name, r_list, 'r')

        #print(round(tuner1, 5), round(tuner2, 5))

        tunez = singleFFT(name, z_list, 'z')

        if tunez != -1:
            FFT_temp_z.append(tunez)
            #FFT_temp_r.append(tuner1)
            #FFT_temp_r2.append(tuner2)
    print(i, len(FFT_temp_z), np.array(FFT_temp_z).mean(), np.array(FFT_temp_z).std())

    FFT_z_list_0.append(np.array(FFT_temp_z).mean())
    #FFT_r_list_0.append(np.array(FFT_temp_r).mean())
    #FFT_r2_list_0.append(np.array(FFT_temp_r2).mean())

    sig0_loss_ratio.append(float(len(FFT_temp_z))/100)

    err_z_list_0.append(np.array(FFT_temp_z).std())
    #err_r_list_0.append(np.array(FFT_temp_r).std())
    #err_r2_list_0.append(np.array(FFT_temp_r2).std())


name = "_3sigma_trans"  # to 1 sigma_z
for i in range(50):
    print("working on " + name + " " + str(i))
    FFT_temp_z = []
    #FFT_temp_r = []
    #FFT_temp_r2 = []
    for index in range(201)[1:]:
        z_list, dE_List, r_list = getData("_3sigma_trans", i*200 + index)

        #tuner1, tuner2 = singleFFT(name, r_list, 'r')

        #print(round(tuner1, 5), round(tuner2, 5))

        tunez = singleFFT(name, z_list, 'z')

        if tunez != -1:
            FFT_temp_z.append(tunez)
            #FFT_temp_r.append(tuner1)
            #FFT_temp_r2.append(tuner2)
    print(i, len(FFT_temp_z), np.array(FFT_temp_z).mean(), np.array(FFT_temp_z).std())

    sig3_loss_ratio.append(float(len(FFT_temp_z))/200)

    FFT_z_list_3.append(np.array(FFT_temp_z).mean())
    #FFT_r_list_3.append(np.array(FFT_temp_r).mean())
    #FFT_r2_list_3.append(np.array(FFT_temp_r2).mean())

    err_z_list_3.append(np.array(FFT_temp_z).std())
    #err_r_list_3.append(np.array(FFT_temp_r).std())
    #err_r2_list_3.append(np.array(FFT_temp_r2).std())



name = "_6sigma_trans"  # to 1 sigma_z
for i in range(50):
    print("working on " + name + " " + str(i))
    FFT_temp_z = []
    #FFT_temp_r = []
    #FFT_temp_r2 = []
    for index in range(101)[1:]:
        z_list, dE_List, r_list = getData(name, i*100 + index)

        #tuner1, tuner2 = singleFFT(name, r_list, 'r')

        tunez = singleFFT(name, z_list, 'z')

        if tunez != -1:
            FFT_temp_z.append(tunez)
            #FFT_temp_r.append(tuner1)
            #FFT_temp_r2.append(tuner2)
    print(i, len(FFT_temp_z), np.array(FFT_temp_z).mean(), np.array(FFT_temp_z).std())

    sig6_loss_ratio.append(float(len(FFT_temp_z))/100)

    FFT_z_list_6.append(np.array(FFT_temp_z).mean())
    #FFT_r_list_6.append(np.array(FFT_temp_r).mean())
    #FFT_r2_list_6.append(np.array(FFT_temp_r2).mean())

    err_z_list_6.append(np.array(FFT_temp_z).std())
    #err_r_list_6.append(np.array(FFT_temp_r).std())
    #err_r2_list_6.append(np.array(FFT_temp_r2).std())


f = plt.figure()

plt.errorbar(np.linspace(0, 4.8, 50), FFT_z_list_0, yerr=err_z_list_0, label='$r_{\perp} = 0 \sigma_{\perp}$', linewidth=3.0)
plt.errorbar(np.linspace(0, 4.8, 50), FFT_z_list_3, yerr=err_z_list_3, label='$r_{\perp} = 3 \sigma_{\perp}$', linewidth=3.0)
plt.errorbar(np.linspace(0, 4.8, 50), FFT_z_list_6, yerr=err_z_list_6, label='$r_{\perp} = 6 \sigma_{\perp}$', linewidth=3.0)
plt.xlabel("z[m]")
plt.ylabel("Synchrotron tune")
plt.legend()

"""
f2 = plt.figure()
plt.errorbar(np.linspace(0, 4.8, 50), FFT_r_list_0, yerr=err_r_list_0, label='$r_{\perp} = 0 \sigma_{\perp}$', linewidth=3.0)
plt.errorbar(np.linspace(0, 4.8, 50), FFT_r_list_3, yerr=err_r_list_3, label='$r_{\perp} = 3 \sigma_{\perp}$', linewidth=3.0)
plt.errorbar(np.linspace(0, 4.8, 50), FFT_r_list_6, yerr=err_r_list_6, label='$r_{\perp} = 6 \sigma_{\perp}$', linewidth=3.0)
plt.xlabel("z[m]")
plt.ylabel("Synchrotron tune")
plt.legend()


f3 = plt.figure()
plt.errorbar(np.linspace(0, 4.8, 50), FFT_r_list_0, yerr=err_r2_list_0, label='$r_{\perp} = 0 \sigma_{\perp}$', linewidth=3.0)
plt.errorbar(np.linspace(0, 4.8, 50), FFT_r2_list_3, yerr=err_r2_list_3, label='$r_{\perp} = 3 \sigma_{\perp}$', linewidth=3.0)
plt.errorbar(np.linspace(0, 4.8, 50), FFT_r2_list_6, yerr=err_r2_list_6, label='$r_{\perp} = 6 \sigma_{\perp}$', linewidth=3.0)
plt.xlabel("z[m]")
plt.ylabel("Synchrotron tune")
plt.legend()
"""

f2 = plt.figure()
plt.plot(np.linspace(0, 4.8, 50), sig0_loss_ratio, label='$r_{\perp} = 0 \sigma_{\perp}$', linewidth=3.0)
plt.plot(np.linspace(0, 4.8, 50), sig3_loss_ratio, label='$r_{\perp} = 3 \sigma_{\perp}$', linewidth=3.0)
plt.plot(np.linspace(0, 4.8, 50), sig6_loss_ratio, label='$r_{\perp} = 6 \sigma_{\perp}$', linewidth=3.0)
plt.xlabel("z[m]")
plt.ylabel("remaining fraction")
plt.legend()

plt.show()
