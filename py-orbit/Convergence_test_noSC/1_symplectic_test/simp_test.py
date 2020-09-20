import numpy as np
from math import *
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


alphax = 5.37828125e-09
betax = 0.7906
alphay = 5.274041667e-09
betay = 1.163203072
emitlimx = 4e-6
emitlimy = 4e-6
gx = (1 + alphax**2)/betax
gy = (1 + alphay**2)/betay


def addAperture2(lattice):

    addTeapotApertureNode(lattice, 0, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 0.5, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 1.5, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 2.1, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 3.05, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 3.6, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 4.45, CircleApertureNode(0.025))

    addTeapotApertureNode(lattice, 5.610940777, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 5.710940777, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 6.610940777, CircleApertureNode(0.025))

    addTeapotApertureNode(lattice, 7.5, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 9, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 9.6, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 10.7, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 11, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 12, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 12.5, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 13.5, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 13.8, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 14.7, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 15.5, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 16.3, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 17, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 18, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 19.1, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 20.15, CircleApertureNode(0.025))

    addTeapotApertureNode(lattice, 21.2, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 22.6, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 23, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 24, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 25, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 25.5, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 26, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 27.5, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 28, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 29, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 30.7, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 31, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 31.9, CircleApertureNode(0.025))

    addTeapotApertureNode(lattice, 33.25, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 33.4530644, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 34.2, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 34.6, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 34.95, CircleApertureNode(0.025))

    addTeapotApertureNode(lattice, 36.4, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 37.8, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 38.3, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 39.5, CircleApertureNode(0.025))


"""
This method read a madx file and set charge
"""
def readLattice(file_name, lattice_name):
    lattice = TEAPOT_Lattice("lattice")

    print "Generate Lattice using MADX (new) parser from .SEQ file"

    lattice.readMADX(file_name,lattice_name)

    lattice.setUseRealCharge(1) # change to + 1 for proton in li ring

    #lattice.setUseRealCharge(1)

    return lattice


# electron bunch for iota

lattice = readLattice("iota_oct4.seq","iota")
#addAperture2(lattice)


for node in lattice.getNodes():
    if node.getType() == "drift teapot":
        node.setnParts(8)
    if node.getType() == "quad teapot":
        node.setnParts(8)
    if node.getType() == "bend teapot":
        node.setnParts(8)
    if node.getType() == "multipole teapot":
        node.setnParts(2)


"""---------------------------- prepare bunch ----------------------------"""
def track_test(factor1, factor2, index):
    b = Bunch()
    b.macroSize(100000)
    energy = 0.0025
    b.mass(0.93827)
    syncPart=b.getSyncParticle()
    syncPart.kinEnergy(energy)

    sx = np.sqrt(betax * emitlimx)
    spx = np.sqrt(gx * emitlimx)
    sy = np.sqrt(betay * emitlimy)
    spy = np.sqrt(gy * emitlimy)


    paramsDict = {}
    lostbunch = Bunch()
    paramsDict["lostbunch"]=lostbunch
    paramsDict["bunch"]= b
    lostbunch.addPartAttr("LostParticleAttributes")

    if index == 0:

        x0 = factor1 * sx - sx * factor2
        px0 = 0 #- spx * 0.01
        y0 = factor1 * sy #- sy * 0.01
        py0 = 0 #- spy * 0.01
        b.addParticle(x0, px0, y0, py0, 0, 0)

        x2 = factor1 * sx + sx * factor2
        px2 = 0 #+ spx * 0.01
        y2 = factor1 * sy #+ sy * 0.01
        py2 = 0 #+ spy * 0.01
        b.addParticle(x2, px2, y2, py2, 0, 0)

        lattice.trackBunch(b, paramsDict)

        wx0 = b.x(0)
        wpx0 = b.px(0)
        wy0 = b.y(0)
        wpy0 = b.py(0)

        wx1 = b.x(1)
        wpx1 = b.px(1)
        wy1 = b.y(1)
        wpy1 = b.py(1)

        return x0, x2, [wx0, wpx0, wy0, wpy0, wx1, wpx1, wy1, wpy1]

    if index == 1:
        x0 = factor1 * sx #- sx * 0.01
        px0 = 0 - spx * factor2
        y0 = factor1 * sy  # - sy * 0.01
        py0 = 0  # - spy * 0.01
        b.addParticle(x0, px0, y0, py0, 0, 0)

        x2 = factor1 * sx #+ sx * 0.01
        px2 = 0 + spx * factor2
        y2 = factor1 * sy  # + sy * 0.01
        py2 = 0  # + spy * 0.01
        b.addParticle(x2, px2, y2, py2, 0, 0)

        lattice.trackBunch(b, paramsDict)

        wx0 = b.x(0)
        wpx0 = b.px(0)
        wy0 = b.y(0)
        wpy0 = b.py(0)

        wx1 = b.x(1)
        wpx1 = b.px(1)
        wy1 = b.y(1)
        wpy1 = b.py(1)

        return px0, px2, [wx0, wpx0, wy0, wpy0, wx1, wpx1, wy1, wpy1]

    if index == 2:
        x0 = factor1 * sx #- sx * 0.01
        px0 = 0  # - spx * 0.01
        y0 = factor1 * sy - sy * factor2
        py0 = 0  # - spy * 0.01
        b.addParticle(x0, px0, y0, py0, 0, 0)

        x2 = factor1 * sx #+ sx * 0.01
        px2 = 0  # + spx * 0.01
        y2 = factor1 * sy + sy * factor2
        py2 = 0  # + spy * 0.01
        b.addParticle(x2, px2, y2, py2, 0, 0)

        lattice.trackBunch(b, paramsDict)

        wx0 = b.x(0)
        wpx0 = b.px(0)
        wy0 = b.y(0)
        wpy0 = b.py(0)

        wx1 = b.x(1)
        wpx1 = b.px(1)
        wy1 = b.y(1)
        wpy1 = b.py(1)

        return y0, y2, [wx0, wpx0, wy0, wpy0, wx1, wpx1, wy1, wpy1]

    if index == 3:
        x0 = factor1 * sx #- sx * 0.01
        px0 = 0  # - spx * 0.01
        y0 = factor1 * sy  # - sy * 0.01
        py0 = 0 - spy * factor2
        b.addParticle(x0, px0, y0, py0, 0, 0)

        x2 = factor1 * sx #+ sx * 0.01
        px2 = 0  # + spx * 0.01
        y2 = factor1 * sy  # + sy * 0.01
        py2 = 0 + spy * factor2
        b.addParticle(x2, px2, y2, py2, 0, 0)

        lattice.trackBunch(b, paramsDict)

        wx0 = b.x(0)
        wpx0 = b.px(0)
        wy0 = b.y(0)
        wpy0 = b.py(0)

        wx1 = b.x(1)
        wpx1 = b.px(1)
        wy1 = b.y(1)
        wpy1 = b.py(1)

        return py0, py2, [wx0, wpx0, wy0, wpy0, wx1, wpx1, wy1, wpy1]


def getderiv(x1, x2, y1, y2):
    return (y2 - y1)/(x2 - x1)

def getJ(factor1, factor2):
    x0, x2, res1 = track_test(factor1, factor2, 0)

    px0, px2, res2 = track_test(factor1, factor2, 1)

    y0, y2, res3 = track_test(factor1, factor2, 2)

    py0, py2, res4 = track_test(factor1, factor2, 3)

    j11 = getderiv(x0, x2, res1[0], res1[4])

    j12 = getderiv(px0, px2, res2[0], res2[4])

    j13 = getderiv(y0, y2, res3[0], res3[4])

    j14 = getderiv(py0, py2, res4[0], res4[4])

    j21 = getderiv(x0, x2, res1[1], res1[5])

    j22 = getderiv(px0, px2, res2[1], res2[5])

    j23 = getderiv(y0, y2, res3[1], res3[5])

    j24 = getderiv(py0, py2, res4[1], res4[5])

    j31 = getderiv(x0, x2, res1[2], res1[6])

    j32 = getderiv(px0, px2, res2[2], res2[6])

    j33 = getderiv(y0, y2, res3[2], res3[6])

    j34 = getderiv(py0, py2, res4[2], res4[6])

    j41 = getderiv(x0, x2, res1[3], res1[7])

    j42 = getderiv(px0, px2, res2[3], res2[7])

    j43 = getderiv(y0, y2, res3[3], res3[7])

    j44 = getderiv(py0, py2, res4[3], res4[7])

    Jacobin = np.matrix([[j11, j12, j13, j14], [j21, j22, j23, j24], [j31, j32, j33, j34], [j41, j42, j43, j44]])

    return Jacobin

"""
sx = np.sqrt(betax * emitlimx)
spx = np.sqrt(gx * emitlimx)
sy = np.sqrt(betay * emitlimy)
spy = np.sqrt(gy * emitlimy)
print(sx, spx, sy, spy)
exit()
"""

factor1 = 5
factor2 = 1e-3
S = np.matrix([[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]])
J = getJ(factor1, factor2)
Jt = np.transpose(J)
print("Jacobin is")
print(J)
print("")
print("det is " + str(np.linalg.det(J)))
print("")
print("J_tran S J is")
print(np.dot(Jt, np.dot(S, J)))



exit()