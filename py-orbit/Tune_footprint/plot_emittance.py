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


filename = './Files/tracking_intensity_10.txt'
turn, emitx, emity, betax, betay, loss = getData(filename)


f = plt.figure()
plt.plot(turn, emitx, label='slice2')
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.legend()
plt.xlabel("turns")
plt.ylabel("emittance_x")

f1 = plt.figure()
plt.plot(turn, emity, label='slice2')
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.legend()
plt.xlabel("turns")
plt.ylabel("emittance_y")

f2 = plt.figure()
plt.plot(turn, betax, label='slice2')
plt.legend()
plt.xlabel("turns")
plt.ylabel("betax")

f3 = plt.figure()
plt.plot(turn, betay, label='slice2')
plt.legend()
plt.xlabel("turns")
plt.ylabel("betay")


f4 = plt.figure()
plt.plot(turn, np.array(loss)/100000, label='slice2')
plt.legend()
plt.xlabel("turns")
plt.ylabel("loss rate")

plt.show()
