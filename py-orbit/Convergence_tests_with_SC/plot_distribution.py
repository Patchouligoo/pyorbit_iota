import matplotlib.pylab as plt
from math import *
import numpy as np
plt.style.use("IceCube")
import matplotlib.colors as colors



alphax = 5.37828125e-09
betax = 0.7905949642
alphay = 5.274041667e-09
betay = 1.163203072
emitlimx = 4.016e-6
emitlimy = 4.016e-6
gx = (1 + alphax**2)/betax
gy = (1 + alphay**2)/betay

sx = np.sqrt(betax * emitlimx)
spx = np.sqrt(gx * emitlimx)
sy = np.sqrt(betay * emitlimy)
spy = np.sqrt(gy * emitlimy)

sz = 1.7
sdE = 0.5e-5

f = open("./intensity_11/n_parts/initial_bunch_slice_03.dat", 'r')
for i in range(16):
    f.readline()
x = []
px = []
y = []
py = []
z = []
dE = []
while 1:
    temp = f.readline().split()

    if len(temp) < 1:
        break

    x.append(float(temp[0])/sx)
    px.append(float(temp[1])/spx)
    y.append(float(temp[2])/sy)
    py.append(float(temp[3])/spy)

    z.append(float(temp[4]))
    dE.append(float(temp[5])/sdE)



f0 = plt.figure()
plt.hist2d(z, dE, bins=[50, 50])
plt.xlabel("z in meter")
plt.ylabel("dE/E in sigmadE/E")
c = plt.colorbar()
#c.set_label('number of particles', rotation=270)


f01 = plt.figure()
plt.hist(x, bins=200)
plt.xlabel("x in sig_x")

f02 = plt.figure()
plt.hist(y, bins=200)
plt.xlabel("y in sig_y")

f1 = plt.figure()
plt.hist2d(x, px, bins=[100, 100], norm=colors.LogNorm())
c = plt.colorbar()
#c.set_label('# particles', rotation=270)
plt.xlabel("x in sig_x")
plt.ylabel("px in sig_px")


f2 = plt.figure()
plt.hist2d(y, py, bins=[100, 100], norm=colors.LogNorm())
c = plt.colorbar()
#c.set_label('# particles', rotation=270)
plt.xlabel("y in sig_y")
plt.ylabel("py in sig_py")

f3 = plt.figure()
plt.hist2d(x, y, bins=[100, 100], norm=colors.LogNorm())
c = plt.colorbar()
#c.set_label('# particles', rotation=270)
plt.xlabel("x in sig_x")
plt.ylabel("y in sig_y")
plt.show()

