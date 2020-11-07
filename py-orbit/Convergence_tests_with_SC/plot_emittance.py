import numpy as np
import matplotlib.pylab as plt
plt.style.use("IceCube")


alphax = 5.37828125e-09
betax = 0.7905949642
alphay = 5.274041667e-09
betay = 1.163203072
emitlimx = 4.016e-6
emitlimy = 4.016e-6
gx = (1 + alphax**2)/betax
gy = (1 + alphay**2)/betay


def getData(filename):
    turn = []
    emittance_x = []
    emittance_y = []
    beta_x = []
    beta_y = []
    paticle_loss = []


    f = open(filename, 'r')

    while 1:
        temp = f.readline().split()

        if len(temp) < 6:
            break
        turn.append(int(temp[0]))
        emittance_x.append(float(temp[1]))
        emittance_y.append(float(temp[2]))
        beta_x.append(float(temp[3]))
        beta_y.append(float(temp[4]))
        paticle_loss.append(float(temp[5]))


    return turn, emittance_x, emittance_y, beta_x, beta_y, paticle_loss




"""
# for num_particles
filename = './intensity_11/grid_size/tracking_64.txt'
turn, emitx, emity, betax, betay, loss = getData(filename)

filename = './intensity_11/grid_size/tracking_128.txt'
turn1, emitx1, emity1, betax1, betay1, loss1 = getData(filename)

filename = './intensity_11/grid_size/tracking_256.txt'
turn2, emitx2, emity2, betax2, betay2, loss2 = getData(filename)


f = plt.figure()
plt.plot(turn, emitx, label='64')
plt.plot(turn1, emitx1, label='128')
plt.plot(turn2, emitx2, label='256')
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.legend()
plt.xlabel("turns")
plt.ylabel("emittance_x")

f1 = plt.figure()
plt.plot(turn, emity, label='64')
plt.plot(turn1, emity1, label='128')
plt.plot(turn2, emity2, label='256')
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.legend()
plt.xlabel("turns")
plt.ylabel("emittance_y")

f2 = plt.figure()
plt.plot(turn, betax, label='64')
plt.plot(turn1, betax1, label='128')
plt.plot(turn2, betax2, label='256')
plt.legend()
plt.xlabel("turns")
plt.ylabel("betax")

f3 = plt.figure()
plt.plot(turn, betay, label='64')
plt.plot(turn1, betay1, label='128')
plt.plot(turn2, betay2, label='256')
plt.legend()
plt.xlabel("turns")
plt.ylabel("betay")


f4 = plt.figure()
plt.plot(turn, 100 * np.array(loss)/100000, label='64')
plt.plot(turn1, 100 * np.array(loss1)/100000, label='128')
plt.plot(turn2, 100 * np.array(loss2)/100000, label='256')
plt.legend()
plt.xlabel("turns")
plt.ylabel("loss rate")

plt.show()
"""




# for num_particles
filename = './intensity_11/num_particles/tracking_100000.txt'
turn, emitx, emity, betax, betay, loss = getData(filename)

filename = './intensity_11/num_particles/tracking_500000.txt'
turn1, emitx1, emity1, betax1, betay1, loss1 = getData(filename)

filename = './intensity_11/num_particles/tracking_1000000.txt'
turn2, emitx2, emity2, betax2, betay2, loss2 = getData(filename)

f = plt.figure()
plt.plot(turn, emitx, label='100000')
plt.plot(turn1, emitx1, label='500000')
plt.plot(turn2, emitx2, label='1000000')
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.legend()
plt.xlabel("turns")
plt.ylabel("emittance_x")

f1 = plt.figure()
plt.plot(turn, emity, label='100000')
plt.plot(turn1, emity1, label='500000')
plt.plot(turn2, emity2, label='1000000')
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.legend()
plt.xlabel("turns")
plt.ylabel("emittance_y")

f2 = plt.figure()
plt.plot(turn, betax, label='100000')
plt.plot(turn1, betax1, label='500000')
plt.plot(turn2, betax2, label='1000000')
plt.legend()
plt.xlabel("turns")
plt.ylabel("betax")

f3 = plt.figure()
plt.plot(turn, betay, label='100000')
plt.plot(turn1, betay1, label='500000')
plt.plot(turn2, betay2, label='1000000')
plt.legend()
plt.xlabel("turns")
plt.ylabel("betay")


f4 = plt.figure()
plt.plot(turn, 100 * np.array(loss)/20000, label='100000')
plt.plot(turn1, 100 * np.array(loss1)/100000, label='500000')
plt.plot(turn2, 100 * np.array(loss2)/200000, label='1000000')
plt.legend()
plt.xlabel("turns")
plt.ylabel("loss rate")

plt.show()





"""
# for num_slice
filename = './intensity_11/n_parts/tracking_slice_05.txt'
turn, emitx, emity, betax, betay, loss = getData(filename)

filename = './intensity_11/n_parts/tracking_slice_03.txt'
turn1, emitx1, emity1, betax1, betay1, loss1 = getData(filename)

filename = './intensity_11/n_parts/tracking_slice_02.txt'
turn2, emitx2, emity2, betax2, betay2, loss2 = getData(filename)

filename = './intensity_11/n_parts/tracking_slice_01.txt'
turn3, emitx3, emity3, betax3, betay3, loss3 = getData(filename)

f = plt.figure()
plt.plot(turn, emitx, label='0.5')
plt.plot(turn1, emitx1, label='0.3')
plt.plot(turn2, emitx2, label='0.2')
plt.plot(turn3, emitx3, label='0.1')
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.legend()
plt.xlabel("turns")
plt.ylabel("emittance_x")

f1 = plt.figure()
plt.plot(turn, emity, label='0.5')
plt.plot(turn1, emity1, label='0.3')
plt.plot(turn2, emity2, label='0.2')
plt.plot(turn3, emity3, label='0.1')
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.legend()
plt.xlabel("turns")
plt.ylabel("emittance_y")

f2 = plt.figure()
plt.plot(turn, betax, label='0.5')
plt.plot(turn1, betax1, label='0.3')
plt.plot(turn2, betax2, label='0.2')
plt.plot(turn3, betax3, label='0.1')
plt.legend()
plt.xlabel("turns")
plt.ylabel("betax")

f3 = plt.figure()
plt.plot(turn, betay, label='0.5')
plt.plot(turn1, betay1, label='0.3')
plt.plot(turn2, betay2, label='0.2')
plt.plot(turn3, betay3, label='0.1')
plt.legend()
plt.xlabel("turns")
plt.ylabel("betay")


f4 = plt.figure()
plt.plot(turn, 100 * np.array(loss)/100000, label='0.5')
plt.plot(turn1, 100 * np.array(loss1)/100000, label='0.3')
plt.plot(turn2, 100 * np.array(loss2)/100000, label='0.2')
plt.plot(turn3, 100 * np.array(loss3)/100000, label='0.1')
plt.legend()
plt.xlabel("turns")
plt.ylabel("loss rate in %")

plt.show()
"""
