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
from orbit.diagnostics import TeapotStatLatsNode, TeapotMomentsNode, TeapotTuneAnalysisNode
from orbit.diagnostics import addTeapotDiagnosticsNode
from orbit.rf_cavities import RFNode, RFLatticeModifications

import orbit_mpi
from orbit_mpi import mpi_datatype
from orbit_mpi import mpi_op
import orbit

comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
rank = orbit_mpi.MPI_Comm_rank(comm)
size = orbit_mpi.MPI_Comm_size(comm)


def addAperture2(lattice):
    addTeapotApertureNode(lattice, 0, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 0.5, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 1.5, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 2.1, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 2.16, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 3.05, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 3.6, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 3.84, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 4.41, CircleApertureNode(0.025))

    addTeapotApertureNode(lattice, 5.610940777, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 5.710940777, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 6.610940777, CircleApertureNode(0.025))

    addTeapotApertureNode(lattice, 7.2, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 7.5, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 7.775, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 9, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 9.6, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 10.12, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 10.7, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 11, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 12, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 12.5, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 13, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 13.5, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 13.8, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 14.7, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 14.95, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 15.5, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 16.3, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 17, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 17.4, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 18, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 19.1, CircleApertureNode(0.025))

    addTeapotApertureNode(lattice, 20.1, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 20.15, CircleApertureNode(0.025))

    addTeapotApertureNode(lattice, 21.2, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 22.57, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 22.6, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 23, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 24, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 25, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 25.5, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 26, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 26.9, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 27.45, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 28, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 29, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 30.12, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 30.7, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 31, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 31.9, CircleApertureNode(0.025))

    addTeapotApertureNode(lattice, 32.5, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 32.8, CircleApertureNode(0.025))

    addTeapotApertureNode(lattice, 33.25, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 33.4530644, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 34.2, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 34.6, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 34.95, CircleApertureNode(0.025))

    addTeapotApertureNode(lattice, 35.5, CircleApertureNode(0.025))

    addTeapotApertureNode(lattice, 36.12, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 36.4, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 37.78, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 38.3, CircleApertureNode(0.025))
    addTeapotApertureNode(lattice, 38.42, CircleApertureNode(0.025))
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
addAperture2(lattice)


"""------------------------------- add RF Cavity -----------------------------"""

"""
ZtoPhi = 2.0 * pi / lattice.getLength()
dESync = 0.0
RFHNum = 4
RFVoltage = 400*pow(10, -9)
RFPhase = 0.0
length = 0.05
name = "rfnode"
rf_node = RFNode.Harmonic_RFNode(ZtoPhi, dESync, RFHNum, RFVoltage, RFPhase, length, name)
position = 5.3
RFLatticeModifications.addRFNode(lattice, position, rf_node)
"""

sizeX = 64   #number of grid points in horizontal direction
sizeY = 64  #number of grid points in vertical direction
sizeZ = 3  #number of longitudinal slices in the 2.5D space charge solver
calc2p5d = SpaceChargeCalc2p5D(sizeX,sizeY,sizeZ)

for node in lattice.getNodes():
    if node.getType() == "drift teapot":
        if node.getLength() > 0.01 and node.getLength() < 0.1:
            node.setnParts(2)
        elif node.getLength() >= 0.1 and node.getLength() < 0.5:
            node.setnParts(3)
        elif node.getLength() > 0.5:
            node.setnParts(8)
    if node.getType() == "quad teapot":
        node.setnParts(8)
    if node.getType() == "bend teapot":
        node.setnParts(8)
    if node.getType() == "multipole teapot":
        node.setnParts(2)

sc_path_length_min = 0.00000001
#scLatticeModifications.setSC2p5DAccNodes(lattice, sc_path_length_min, calc2p5d)


""" ------------------------- prepare the bunch ---------------------------"""
num_charge_max = 20000


b = Bunch()
b.macroSize(num_charge_max/400)
energy = 0.0025
b.mass(0.93827)
syncPart=b.getSyncParticle()
syncPart.kinEnergy(energy)
b.addPartAttr("ParticleIdNumber")

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

f_list = []

if 1:
    particle_num = 10000
    n = 0
    for i in range(100):
        for j in range(100):
            print("initializing: " + str(n))

            x = float(i + 1)/100 * 5 * sx

            y = float(j + 1)/100 * 5 * sy

            b.addParticle(x, 0, y, 0, 0, 0)

            b.partAttrValue("ParticleIdNumber", n, 0, n)

            f = open("./track_data/" + str(n) + ".txt", 'w')

            f_list.append(f)

            n += 1



"""--------------------------------- end ---------------------------------"""

paramsDict = {}
lostbunch = Bunch()
paramsDict["lostbunch"]=lostbunch
paramsDict["bunch"]= b



ana = BunchTwissAnalysis()


for i in range(1000):
    if rank == 0:
        print("==========================")
        print("turn: " + str(i))
        print("lost: " + str(lostbunch.getSize()))
        print("remain: " + str(b.getSize()))
        print("charge: " + str(b.macroSize()))
        print("intensity: " + str(b.macroSize() * b.getSize() * size))
        print("==========================")
    lattice.trackBunch(b, paramsDict)


for i in range(5000):

    ana.analyzeBunch(b)

    #if rank == 0:
    #    f.write(str(i) + " " + str(ana.getEmittanceNormalized(0)) + " " + str(ana.getEmittanceNormalized(1)) + " " + str(ana.getBeta(0)) + " " + str(ana.getBeta(1)) + " " + str(lostbunch.getSize()) + "\n")

    lattice.trackBunch(b, paramsDict)

    if rank == 0:
        print("==========================")
        print("turn: " + str(i))
        print("lost: " + str(lostbunch.getSize()))
        print("remain: " + str(b.getSize()))
        print("charge: " + str(b.macroSize()))
        print("intensity: " + str(b.macroSize() * b.getSize() * size))
        print(str(ana.getEmittanceNormalized(0)) + " " + str(ana.getEmittanceNormalized(1)))
        print("==========================")

    for j in range(b.getSize()):
        id = int(b.partAttrValue("ParticleIdNumber", j, 0))

        f = f_list[id]

        f.write(str(i) + " " + str(b.x(j)) + " " + str(b.y(j)) + "\n")


