"""
Module. Includes functions that will modify the accelerator lattice by inserting the one teapot node accelerator node.
"""
# import the auxiliary classes
from orbit.utils import orbitFinalize

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer, AccNodeBunchTracker

# import Teapot Aperture node
from aperture import Aperture
from TeapotApertureNode import TeapotApertureNode, CircleApertureNode, EllipseApertureNode, RectangleApertureNode
from ApertureLatticeModifications import addTeapotApertureNode

# import teapot drift class
from orbit.teapot import DriftTEAPOT

#This create a set of circular apertures. a is the radius of the apertures, s is the starting position, e is the ending position, and  c is the x offset and d is the y offset of the apertures.
def addCircleApertureSet(a, lattice, s = 0, e = 0, c = 0, d = 0):

	if e > lattice.getLength():
		e = lattice.getLength()
		print 'Warning, end position exceeding lattice length. Resetting to lattice length.'
	
	position = s
	position_list = []
	for node in lattice.getNodes():
			if(isinstance(node,DriftTEAPOT) and node.getLength() >= 0.01):
					positiont = lattice.getNodePositionsDict()[node][0]
					if positiont > position and positiont < e:
							position = positiont
							position_list.append(position + 0.001)
					positiont = lattice.getNodePositionsDict()[node][-1]
					if positiont > position and positiont < e:
							position = positiont
							position_list.append(position - 0.001)
                                        #print(position)
					#Aperturenode = CircleApertureNode(a, c, d)
					#addTeapotApertureNode(lattice, position, Aperturenode)
	for pos in position_list:
                #print(pos)
                Aperturenode = CircleApertureNode(a, c, d)
		addTeapotApertureNode(lattice, pos, Aperturenode)


#This create a set of eliptic apertures. a is the x radius and b is the y radius of the apertures, s is the starting position, e is the ending position, and c is the x offset and d is the y offset of the apertures.
def addEllipseApertureSet(a, b, lattice, s = 0, e = 0, c = 0, d = 0):

	if e > lattice.getLength():
		e = lattice.getLength()
		print 'Warning, end position exceeding lattice length. Resetting to lattice length.'
        position = s
	position_list = []
	for node in lattice.getNodes():
			if(isinstance(node,DriftTEAPOT) and node.getLength() >= 0.01):
					positiont = lattice.getNodePositionsDict()[node][0]
					if positiont > position and positiont < e:
							position = positiont
							position_list.append(position + 0.001)
					positiont = lattice.getNodePositionsDict()[node][-1]
					if positiont > position and positiont < e:
							position = positiont
							position_list.append(position - 0.001)
                                        #print(position)
					#Aperturenode = CircleApertureNode(a, c, d)
					#addTeapotApertureNode(lattice, position, Aperturenode)
	for pos in position_list:
                Aperturenode = EllipseApertureNode(a, b, c, d)
                addTeapotApertureNode(lattice, position, Aperturenode)



#This create a set of rectangular apertures. a is the x half width and b is the y half hight of the apertures, s is the starting position, e is the ending position, and c is the x offset and d is the y offset of the apertures.
def addRectangleApertureSet(a, b, lattice, s = 0, e = 0, c = 0, d = 0):

	if e > lattice.getLength():
		e = lattice.getLength()
		print 'Warning, end position exceeding lattice length. Resetting to lattice length.'
	
	position = s
	position_list = []
	for node in lattice.getNodes():
			if(isinstance(node,DriftTEAPOT) and node.getLength() >= 0.01):
					positiont = lattice.getNodePositionsDict()[node][0]
					if positiont > position and positiont < e:
							position = positiont
							position_list.append(position + 0.001)
					positiont = lattice.getNodePositionsDict()[node][-1]
					if positiont > position and positiont < e:
							position = positiont
							position_list.append(position - 0.001)
                                        #print(position)
					#Aperturenode = CircleApertureNode(a, c, d)
					#addTeapotApertureNode(lattice, position, Aperturenode)
	for pos in position_list:
                Aperturenode = RectangleApertureNode(a, b, c, d)
		addTeapotApertureNode(lattice, position, Aperturenode)



