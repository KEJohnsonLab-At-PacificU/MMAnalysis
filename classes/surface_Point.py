##########################################################################
# This is a class definition for a surface point. This contains the data
# for a point in a matrix describing a surface, for example the surface of
# water at the interface in a water/organics simulation.
# Author: Chris McCormick
# Date: Summer 2016
##########################################################################

from global_Vars import *
from numpy import sqrt, square

class surface_Point(object):
    def __init__(self):
        self.atom_Type = None
        self.z_Intersect = None
        
    def calc_Intersect(self, atom, x_Index, y_Index):
        intersect = sqrt(square(atom.VDW) - square(atom.position[X] * 10 - x_Index) - square(atom.position[Y] * 10 - y_Index)) + (atom.position[Z] * 10)
        return intersect
    
    # a rough attempt at organic bottom surface. I may comment this out if it doesn't work in time
    def calc_Intersect_Inverse(self, atom, x_Index, y_Index):
        intersect = sqrt(square(atom.VDW) - square(atom.position[X] * 10 - x_Index) - square(atom.position[Y] * 10 - y_Index)) + (atom.position[Z] * -10)
        return intersect