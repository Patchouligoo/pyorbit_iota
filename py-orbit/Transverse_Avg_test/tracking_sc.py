import numpy as np
from math import *
import random
import os
from bunch import Bunch
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import TEAPOT_MATRIX_Lattice
from bunch import BunchTuneAnalysis
from orbit.utils.orbit_mpi_utils import bunch_orbit_to_pyorbit, bunch_pyorbit_to_orbit
from orbit.matrix_lattice import MATRIX_Lattice
from orbit.space_charge.sc2p5d import scAccNodes, scLatticeModifications
from spacecharge import SpaceChargeCalc2p5D, Boundary2D, LSpaceChargeCalc
from orbit.aperture import TeapotApertureNode, CircleApertureNode, EllipseApertureNode, RectangleApertureNode
from orbit.aperture import addTeapotApertureNode, addCircleApertureSet
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


"""
used when initializing the z-dE distribution. This method will give the abs(dE)_max from the contour
we choose, as represented by the z_list and dE_list
"""
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

"""
This method will check if the given z, dE position is outside the separatrix modeled by a cos function
"""
def check_bucket(z, E):
    if abs(z) > 5:
        return False
    E_max = 1.9e-05 * cos(z/5 * pi/2)

    #E_max = 2.64e-05 * cos(z / 5 * pi / 2)

    if abs(E) >  E_max:
        return False
    return True


"""
This method read a madx file and set charge
"""
def readLattice(file_name, lattice_name):
    lattice = TEAPOT_Lattice("lattice")

    print "Generate Lattice using MADX (new) parser from .SEQ file"

    lattice.readMADX(file_name,lattice_name)

    lattice.setUseRealCharge(1) # change to + 1 for proton

    return lattice


# electron bunch for iota

lattice = readLattice("./iota_linear.seq","iota")
# add circular aperture 25mm to all elements
addCircleApertureSet(0.1, lattice, 0, 40)


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
sizeZ = 13  #number of longitudinal slices in the 2.5D space charge solver
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


"""--------------------initialize the bunch-----------------"""

Intensity = 9*pow(10, 10)

# num of MPs
particle_num_total = 500000

# store the output file as: pre_list/test_name.txt
pre_list = './9e10_avg_test'
test_name = 'tracking_transverse_tune'


num_charge_max = Intensity / particle_num_total


b = Bunch()
b.macroSize(num_charge_max/200)  # charge/200 for 200 turns of slow initialization
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

sx = sqrt(betax * emitlimx)
spx = sqrt(gx * emitlimx)
sy = sqrt(betay * emitlimy)
spy = sqrt(gy * emitlimy)


DX = -3.32348914 * b.getSyncParticle().beta()
DPX = 9.759764092e-08 * (b.getSyncParticle().beta() ** 2)


# read the contour defined in z_tuple file
# todo: 3
p = open("./z_tuple_45.txt", 'r')
z_list = []
dE_list = []
while 1:
    temp = p.readline().split()
    if len(temp) < 1:
        break
    z_list.append(float(temp[0]))
    dE_list.append(float(temp[1]))



"""---------------------- add particles -------------------------"""
if 1:
    particle_num = particle_num_total / size

    for n in range(particle_num):
        print("initializing: " + str(n))
        while 1:

            x = random.gauss(0, sx)
            px = random.gauss(0, spx)

            y = random.gauss(0, sy)
            py = random.gauss(0, spy)

            z = random.random() * (4.5 * 2) - 4.5

            dE = (random.random() * (1.84e-5 * 2) - 1.84e-5)

            test_E = compare(z_list, dE_list, z)

            dP_P = (dE / 0.94077) * pow(b.getSyncParticle().beta(), -2)
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



if rank == 0:
    if not os.path.exists(pre_list):
        os.mkdir(pre_list)

# record initial bunch
b.dumpBunch(pre_list + "/initial_bunch_" + str(test_name) + ".dat")


ana = BunchTwissAnalysis()

# record particle loss due to debunch
if rank == 0:
    f = open(pre_list + "/tracking_" + str(test_name) + ".txt", 'w')
    f2 = open(pre_list + "/z_loss_" + str(test_name) + ".txt", 'w')


# create folder to store lost particles using CPU0
if rank == 0:
    if not os.path.exists(pre_list + "/test_data"):
        os.mkdir(pre_list + "/test_data")




"""--------------- tracking for n turns----------------"""
for i in range(501):

    ana.analyzeBunch(b)

    if rank == 0:
        f.write(str(i) + " " + str(ana.getEmittanceNormalized(0)) + " " + str(ana.getEmittanceNormalized(1)) + " " + str(ana.getBeta(0)) + " " + str(ana.getBeta(1)) + " " + str(lostbunch.getSize()) + "\n")

    lattice.trackBunch(b, paramsDict)

    # record using CPU0
    if rank == 0:
        print("==========================")
        print("turn: " + str(i))
        print("lost: " + str(lostbunch.getSize()))
        print("remain: " + str(b.getSize()))
        print("charge: " + str(b.macroSize()))
        print("intensity: " + str(b.macroSize() * b.getSize() * size))
        print(str(ana.getEmittanceNormalized(0)) + " " + str(ana.getEmittanceNormalized(1)))
        print("==========================")



    # check for debunch particles, then remove and record them
    if 1:
        k = 0
        while k < b.getSize():
            if check_bucket(b.z(k), b.dE(k)) is False:
                if rank == 0:
                    f2.write(str(i) + " " + str(b.partAttrValue("ParticleIdNumber", k, 0)) + " " + str(b.z(k)) + " " + str(b.dE(k)) + "\n")
                b.deleteParticle(k)
            else:
                k += 1

        if rank == 0:
            f2.flush()


    # slow initialization for the first 200 turns
    if i < 200:
        b.macroSize(num_charge_max/200 * (i + 1))



b.dumpBunch(pre_list + "/ready_bunch_" + str(test_name) + ".dat")




"""----------------------add test particles----------------------"""

ana.analyzeBunch(b)
curr_emitx = ana.getEmittance(0)
curr_betax = ana.getBeta(0)
curr_alphax = ana.getAlpha(0)

curr_emity = ana.getEmittance(1)
curr_betay = ana.getBeta(1)
curr_alphay = ana.getAlpha(1)


if rank == 0:


    count = 1


    for n_factor in np.linspace(0, 3, 51)[1:]:

        for theta in np.linspace(0, 2*pi, 101)[:-1]:

            Jx = n_factor**2 * curr_emitx / 2

            Jy = n_factor**2 * curr_emity / 2

            x_temp = sqrt(2 * curr_betax * Jx) * cos(theta)
            px_temp = -1 * sqrt(2 * Jx / curr_betax) * (sin(theta) + curr_alphax * cos(theta))

            y_temp = sqrt(2 * curr_betay * Jy) * cos(theta)
            py_temp = -1 * sqrt(2 * Jy / curr_betay) * (sin(theta) + curr_alphay * cos(theta))


            b.addParticle(x_temp, px_temp, y_temp, py_temp, 0, 0)

            p = open(pre_list + "/test_data/test_" + str(count) + ".txt", 'w')

            p.close()

            b.partAttrValue("ParticleIdNumber", b.getSize() - 1, 0, -1 * count)

            count += 1




for i in range(200):

    ana.analyzeBunch(b)

    if rank == 0:
        f.write(str(i + 501) + " " + str(ana.getEmittanceNormalized(0)) + " " + str(ana.getEmittanceNormalized(1)) + " " + str(ana.getBeta(0)) + " " + str(ana.getBeta(1)) + " " + str(lostbunch.getSize()) + "\n")

        count_num_left = 0

        for j in range(b.getSize()):
            index = int(b.partAttrValue("ParticleIdNumber", j, 0))

            if index < 0:

                count_num_left += 1

                p = open(pre_list + "/test_data/test_" + str(abs(index)) + ".txt", 'a')

                p.write(str(b.x(j)) + " " + str(b.px(j)) + " " + str(b.y(j)) + " " + str(b.py(j)) + " " + str(b.z(j)) + " " + str(b.dE(j)) + "\n")

                p.close()

        print("number of test particles left: " + str(count_num_left))


    lattice.trackBunch(b, paramsDict)

    if rank == 0:
        print("==========================")
        print("turn: " + str(i + 501))
        print("lost: " + str(lostbunch.getSize()))
        print("remain: " + str(b.getSize()))
        print("charge: " + str(b.macroSize()))
        print("intensity: " + str(b.macroSize() * b.getSize() * size))
        print(str(ana.getEmittanceNormalized(0)) + " " + str(ana.getEmittanceNormalized(1)))
        print("==========================")



    if 1:
        k = 0
        while k < b.getSize():
            if check_bucket(b.z(k), b.dE(k)) is False:
                if rank == 0:
                    f2.write(str(i + 501) + " " + str(b.partAttrValue("ParticleIdNumber", k, 0)) + " " + str(b.z(k)) + " " + str(b.dE(k)) + "\n")
                b.deleteParticle(k)
            else:
                k += 1

        if rank == 0:
            f2.flush()


b.dumpBunch(pre_list + "/final_bunch_" + str(test_name) + ".dat")

lostbunch.dumpBunch(pre_list + "/lost_bunch_" + str(test_name) + ".dat")




