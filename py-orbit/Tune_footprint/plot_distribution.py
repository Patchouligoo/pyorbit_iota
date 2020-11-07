
import matplotlib.pylab as plt
from math import *
import numpy as np
plt.style.use("IceCube")
import matplotlib.colors as colors

"""

from math import *
import os
from bunch import Bunch
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import TEAPOT_MATRIX_Lattice
from bunch import BunchTuneAnalysis
from orbit.utils.orbit_mpi_utils import bunch_orbit_to_pyorbit, bunch_pyorbit_to_orbit
from orbit.matrix_lattice import MATRIX_Lattice
from orbit.space_charge.sc2p5d import scAccNodes, scLatticeModifications
from spacecharge import SpaceChargeCalc2p5D, Boundary2D
from orbit.aperture import TeapotApertureNode, CircleApertureNode, EllipseApertureNode, RectangleApertureNode
from orbit.aperture import addTeapotApertureNode
from bunch import BunchTwissAnalysis
from orbit.diagnostics import TeapotStatLatsNode, TeapotMomentsNode, TeapotTuneAnalysisNode
from orbit.diagnostics import addTeapotDiagnosticsNode
from orbit.rf_cavities import RFNode, RFLatticeModifications

import orbit_mpi
from orbit_mpi import mpi_datatype
from orbit_mpi import mpi_op
import orbit

b = Bunch()

b.readBunch("./Files/initial_bunch_intensity_10.dat")

z_max = 0
z_min = 0

for i in range(b.getSize()):
    if b.z(i) > z_max:
        z_max = b.z(i)

    if b.z(i) < z_min:
        z_min = b.z(i)

interval = z_max - z_min


up_bound = float(z_max-z_min)/6
low_bound = -1 * float(z_max-z_min)/6

print(up_bound, low_bound)

k = 0
while k < b.getSize():
    print(k)
    if b.z(k) < up_bound and b.z(k) > low_bound:
        k += 1
    else:
        b.deleteParticle(k)

ana = BunchTwissAnalysis()
ana.analyzeBunch(b)


print("size of MPs: " + str(b.getSizeGlobal()))

print("emittance: x = " + str(ana.getEmittanceNormalized(0)) + ", y = " + str(ana.getEmittanceNormalized(1)))

#b.dumpBunch("WTF.dat")

exit()
"""

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

f = open("./Files/final_bunch_intensity_10.dat", 'r')
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


f00 = plt.figure()
N, bins, patches = plt.hist(z, bins=61)
"""
for i in range(0,20):
    patches[i].set_facecolor('b')
for i in range(20,40):
    patches[i].set_facecolor('r')
for i in range(40, len(patches)):
    patches[i].set_facecolor('b')
"""
plt.xlabel("z in meter")


f0 = plt.figure()
plt.hist2d(z, dE, bins=[50, 50])
plt.xlabel("z in meter")
plt.ylabel("dE/E in sigmadE/E")
c = plt.colorbar()
c.set_label('number of particles', rotation=270)


f01 = plt.figure()
plt.hist(x, bins=200)
plt.xlabel("x in sig_x")

f02 = plt.figure()
plt.hist(y, bins=200)
plt.xlabel("y in sig_y")

f1 = plt.figure()
plt.hist2d(x, px, bins=[100, 100], norm=colors.LogNorm())
c = plt.colorbar()
c.set_label('# particles', rotation=270)
plt.xlabel("x in sig_x")
plt.ylabel("px in sig_px")


f2 = plt.figure()
plt.hist2d(y, py, bins=[100, 100], norm=colors.LogNorm())
c = plt.colorbar()
c.set_label('# particles', rotation=270)
plt.xlabel("y in sig_y")
plt.ylabel("py in sig_py")

f3 = plt.figure()
plt.hist2d(x, y, bins=[100, 100], norm=colors.LogNorm())
c = plt.colorbar()
c.set_label('# particles', rotation=270)
plt.xlabel("x in sig_x")
plt.ylabel("y in sig_y")
plt.show()



z = np.array(z)

z_max = z.max()
z_min = z.min()

interval_up = (z_max - z_min)/6
interval_down = -1 * (z_max - z_min)/6

print(interval_down, interval_up)

count = 0
for pos in z:
    if pos > interval_down and pos < interval_up:
        count += 1

print(count)
print(float(count) * 20000)

