##########################################################################
# This is a class definition for an atom. Note that the units are actually
# whatever is output by openMM, but the defaults are indicated in comments
# here
# Author: Chris McCormick
# Date: Summer 2016
##########################################################################
from global_Vars import *

class atom_Data(object):
    #def __init__(self, atom_ID, atom_Name, atom_Mass, position, velocity, force):
    def __init__(self, atom_Name, atom_Mass, position, velocity, force):
        #self.atom_ID = atom_ID			# we don't use this now, but it might be nice later
        self.atom_Name = atom_Name      # in letters (wah-wah)
        self.atom_Mass = atom_Mass		# in daltons
        self.position = position		# in nanometers
        self.velocity = velocity		# in nanometers/picosecond
        self.force = force				# in kilojoules/mole/nanometer
        
    ''' This is our aborted attempt to calculate interfacial pressure. We found that these pressure values
    grew to large for the program to handle. The source for these equations does not explain them very clearly,
    so we are not sure what we are missing.
    def calc_Pressure(self):
        pressure_X = (self.atom_Mass * self.velocity[X] * self.velocity[X]) + (self.position[X] * self.force[X])
        pressure_Y = (self.atom_Mass * self.velocity[Y] * self.velocity[Y]) + (self.position[Y] * self.force[Y])
        pressure_Z = (self.atom_Mass * self.velocity[Z] * self.velocity[Z]) + (self.position[Z] * self.force[Z])
        self.pressure = (pressure_X, pressure_Y, pressure_Z)
    '''