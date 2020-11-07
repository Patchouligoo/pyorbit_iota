import numpy as np
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

comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
rank = orbit_mpi.MPI_Comm_rank(comm)
size = orbit_mpi.MPI_Comm_size(comm)


def compare(z_list, dE_list, z):
    index = 0
    for i in range(len(z_list) - 1):
        if z < z_list[i]:
            index = i
            break

    E_min = dE_list[index - 1]
    E_max = dE_list[index]

    slope = (E_max - E_min)/(z_list[index] - z_list[index - 1])

    dE_max = slope * (z - z_list[index - 1]) + E_min

    return dE_max

def check_bucket(z, E):
    if abs(z) > 6:
        return False
    E_max = 1.85e-5 * cos(z/4.7 * pi/2)

    if abs(E) >  E_max:
        return False
    return True


def addAperture2(lattice):

    ape_size = 0.06

    addTeapotApertureNode(lattice, 0, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 0.5, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 1.5, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 2.1, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 2.16, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 3.05, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 3.6, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 3.84, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 4.41, CircleApertureNode(ape_size))

    addTeapotApertureNode(lattice, 5.610940777, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 5.710940777, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 6.610940777, CircleApertureNode(ape_size))

    addTeapotApertureNode(lattice, 7.2, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 7.5, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 7.775, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 9, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 9.6, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 10.12, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 10.7, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 11, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 12, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 12.5, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 13, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 13.5, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 13.8, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 14.7, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 14.95, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 15.5, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 16.3, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 17, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 17.4, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 18, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 19.1, CircleApertureNode(ape_size))

    addTeapotApertureNode(lattice, 20.1, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 20.15, CircleApertureNode(ape_size))

    addTeapotApertureNode(lattice, 21.2, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 22.57, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 22.6, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 23, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 24, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 25, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 25.5, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 26, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 26.9, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 27.45, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 28, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 29, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 30.12, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 30.7, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 31, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 31.9, CircleApertureNode(ape_size))

    addTeapotApertureNode(lattice, 32.5, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 32.8, CircleApertureNode(ape_size))

    addTeapotApertureNode(lattice, 33.25, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 33.4530644, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 34.2, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 34.6, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 34.95, CircleApertureNode(ape_size))

    addTeapotApertureNode(lattice, 35.5, CircleApertureNode(ape_size))

    addTeapotApertureNode(lattice, 36.12, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 36.4, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 37.78, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 38.3, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 38.42, CircleApertureNode(ape_size))
    addTeapotApertureNode(lattice, 39.5, CircleApertureNode(ape_size))


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

lattice = readLattice("iota_linear_shifted.seq","iota")
addAperture2(lattice)


""" ----------------- add rf cavity ------------------ """
ZtoPhi = 2.0 * pi / lattice.getLength()
dESync = 0.0
RFHNum = 4
RFVoltage = 400*pow(10, -9)
RFPhase = 0.0
length = 0.05
name = "rfnode"
rf_node = RFNode.Harmonic_RFNode(ZtoPhi, dESync, RFHNum, RFVoltage, RFPhase, length, name)
position = 25.79485606
RFLatticeModifications.addRFNode(lattice, position, rf_node)


""" ----------------- add SC effect -------------------"""
sizeX = 128   #number of grid points in horizontal direction
sizeY = 128  #number of grid points in vertical direction
sizeZ = 3  #number of longitudinal slices in the 2.5D space charge solver
calc2p5d = SpaceChargeCalc2p5D(sizeX,sizeY,sizeZ)


interval = 0.2

for node in lattice.getNodes():
    if node.getType() == "drift teapot":
        node.setnParts(int(float(node.getLength())/interval) + 1)
    if node.getType() == "quad teapot":
        node.setnParts(int(float(node.getLength())/interval) + 2)
    if node.getType() == "bend teapot":
        node.setnParts(int(float(node.getLength())/interval) + 2)
    if node.getType() == "multipole teapot":
        node.setnParts(int(float(node.getLength())/interval) + 2)
    if node.getType() == "solenoid teapot":
        node.setnParts(int(float(node.getLength())/interval) + 2)



sc_path_length_min = 0.00000001
scLatticeModifications.setSC2p5DAccNodes(lattice, sc_path_length_min, calc2p5d)



count = 0
dis = 0
pos_list = []
for node in lattice.getNodes():

    for sub in node.getAllChildren():
        #print(sub.getType() + " " + str(dis))

        if sub.getType() == "SC2p5D":
            count += 1
            dis = dis + sub.getLengthOfSC()
            pos_list.append(dis)
print(len(pos_list))
file_writer = open("./SC_KICKS.txt", 'w')

for i in pos_list:
    file_writer.write(str(i) + "\n")

file_writer.close()



num_charge_max = 20000


b = Bunch()
b.macroSize(num_charge_max/200)
energy = 0.0025
b.mass(0.93827)
syncPart=b.getSyncParticle()
syncPart.kinEnergy(energy)
b.addPartAttr("ParticleIdNumber")

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


"""--------------------------------- end ---------------------------------"""

"""
initialize particles
"""
DX = -3.32348914 * b.getSyncParticle().beta()
DPX = 9.759764092e-08 * (b.getSyncParticle().beta() ** 2)


p = open("./z_tuple.txt", 'r')
z_list = []
dE_list = []
while 1:
    temp = p.readline().split()
    if len(temp) < 1:
        break
    z_list.append(float(temp[0]))
    dE_list.append(float(temp[1]))



if 1:
    particle_num = 100000
    for n in range(particle_num):
        print("initializing: " + str(n))
        while 1:

            x = np.random.normal(0, sx)
            px = np.random.normal(0, spx)

            y = np.random.normal(0, sy)
            py = np.random.normal(0, spy)

            z = np.random.random() * (3.5 * 2) - 3.5
            dE = (np.random.random() * (1.84e-5 * 2) - 1.84e-5)

            test_E = compare(z_list, dE_list, z)

            dP_P = (dE / 0.94077) * np.power(b.getSyncParticle().beta(), -2)
            x = x + dP_P * DX
            px = px + dP_P * DPX

            if abs(dE) < test_E: # and sqrt((x/sx)**2 + (y/sy)**2) < 1:
                break

        b.addParticle(x, px, y, py, z, dE)

        b.partAttrValue("ParticleIdNumber", n, 0, n + particle_num * rank)


#b = Bunch()
#b.readBunch('./intensity_11/n_parts/initial_bunch.dat')


paramsDict = {}
lostbunch = Bunch()
paramsDict["lostbunch"]=lostbunch
paramsDict["bunch"]= b


#todo:
pre_list = './Files_shifted'
test_name = 'intensity_10'



b.dumpBunch(pre_list + "/initial_bunch_" + str(test_name) + ".dat")


ana = BunchTwissAnalysis()


f = open(pre_list + "/tracking_" + str(test_name) + ".txt", 'w')


for i in range(500):

    ana.analyzeBunch(b)

    if rank == 0:
        f.write(str(i) + " " + str(ana.getEmittanceNormalized(0)) + " " + str(ana.getEmittanceNormalized(1)) + " " + str(ana.getBeta(0)) + " " + str(ana.getBeta(1)) + " " + str(lostbunch.getSize()) + "\n")

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

    if i%100 == 0:
        lostbunch.dumpBunch(pre_list + "/lost_data/turn_" + str(i) + ".dat")


    if 1:
        k = 0
        while k < b.getSize():
            if check_bucket(b.z(k), b.dE(k)) is False:
                b.deleteParticle(k)
            else:
                k += 1

    if i < 200:
        b.macroSize(num_charge_max/200 * (i + 1))







ana.analyzeBunch(b)
curr_emitx = ana.getEmittance(0)
curr_betax = ana.getBeta(0)
curr_sx = sqrt(curr_betax * curr_emitx)

if rank == 0:
    if not os.path.exists(pre_list + "/test_particles_" + str(test_name)):
        os.mkdir(pre_list + "/test_particles_" + str(test_name))


    p_list = []
    count = 1
    for ang in np.linspace(0, pi/2, 30):
        for r in np.linspace(0, 10*curr_sx, 100):

            x = cos(ang) * r
            px = 0
            if x == 0:
                x = curr_sx * 0.01

            y = sin(ang) * r
            py = 0
            if y == 0:
                y = curr_sx * 0.01

            b.addParticle(x, px, y, py, 0, 0)

            p = open(pre_list + "/test_particles_" + str(test_name) + "/" + str(count) + ".txt", 'w')

            p_list.append(p)

            b.partAttrValue("ParticleIdNumber", b.getSize()-1, 0, -1 * count)
            count += 1



for i in range(2000):

    ana.analyzeBunch(b)

    if rank == 0:
        f.write(str(i + 500) + " " + str(ana.getEmittanceNormalized(0)) + " " + str(ana.getEmittanceNormalized(1)) + " " + str(ana.getBeta(0)) + " " + str(ana.getBeta(1)) + " " + str(lostbunch.getSize()) + "\n")

        for j in range(b.getSize()):
            index = int(b.partAttrValue("ParticleIdNumber", j, 0))

            if index < 0:
                p = p_list[abs(index) - 1]
                p.write(str(b.x(j)) + " " + str(b.y(j)) + "\n")


    lattice.trackBunch(b, paramsDict)

    if rank == 0:
        print("==========================")
        print("turn: " + str(i + 500))
        print("lost: " + str(lostbunch.getSize()))
        print("remain: " + str(b.getSize()))
        print("charge: " + str(b.macroSize()))
        print("intensity: " + str(b.macroSize() * b.getSize() * size))
        print(str(ana.getEmittanceNormalized(0)) + " " + str(ana.getEmittanceNormalized(1)))
        print("==========================")

    if i%100 == 0:
        lostbunch.dumpBunch(pre_list + "/lost_data/turn_" + str(i + 500) + ".dat")


    if 1:
        k = 0
        while k < b.getSize():
            if check_bucket(b.z(k), b.dE(k)) is False:
                b.deleteParticle(k)
            else:
                k += 1


b.dumpBunch(pre_list + "/final_bunch_" + str(test_name) + ".dat")
