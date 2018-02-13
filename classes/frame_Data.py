##########################################################################
# This is a class definition for a simulation - this contains data about a
# specific frame, or 'step' being analysed.
# Author: Chris McCormick
# Date: Summer 2016
##########################################################################

import numpy as np
from classes.z_Slice import Z
from classes.z_Slice import z_Slice
from global_Vars import *

class frame_Data(object):

    def __init__(self, frame_ID, x_len, y_len, z_len):
        self.frame_ID = frame_ID
        self.x_len = x_len # nanometers
        self.y_len = y_len # nanometers
        self.z_len = z_len # nanometers
        '''
        commented out are attributes that are now only assigned when needed within other function calls.
        again, this is because I wanted to see the impact on memory usage, and it seems to have helped a bit.
        I keep the comments here just as a reminded of what attributes can exist
        '''
        # self.interface_Z_Pos = None       # nanometers
        # self.num_Negative_Slices = None
        # self.num_Positive_Slices = None
        self.molecules = []
        self.slices_Above_Interface = []
        self.slices_Below_Interface = []
        # we hope to calculate and store the interface pressure
        
    # Figures out what the top (128) oxygen atoms are in the simulation, and returns the average of their z-position
    # which is used as an approximation of the interface location, and the z-origin in analysis.
    def get_Interface_Z(self):
        all_Oxy_Z_Values = []
        num_Interface_Atoms = 128
        for molecule in self.molecules:
            if molecule.molecule_Name == 'HOH':
                for atom in molecule.atoms:
                    if atom.atom_Name == 'O':
                        all_Oxy_Z_Values.append(atom.position[Z])
        all_Oxy_Z_Values = np.sort(all_Oxy_Z_Values, kind = 'mergesort')
        self.interface_Z_Pos = np.average(all_Oxy_Z_Values[-num_Interface_Atoms:])
        return self.interface_Z_Pos
    
    # Figures out the proportion of space above and below the interface, then creates the appropriate number
    # of slices for each.
    def build_Slices(self, my_Type_List, my_Atom_Groups):
        type_List = my_Type_List
        atom_Groups = my_Atom_Groups
        
        num_Negative_Slices = int(NUM_SLICES * (self.interface_Z_Pos / self.z_len))
        num_Positive_Slices = NUM_SLICES - num_Negative_Slices
        for i in range(num_Negative_Slices):
            self.slices_Below_Interface.append(z_Slice(((i + 1) * -1), self.interface_Z_Pos - ((i + 1) * SLICE_SIZE), self.x_len, self.y_len, type_List, atom_Groups))
        for i in range(num_Positive_Slices):
            self.slices_Above_Interface.append(z_Slice(i, self.interface_Z_Pos + (i * SLICE_SIZE), self.x_len, self.y_len, type_List, atom_Groups))

    '''
    This is our aborted attempt to calculate interfacial pressure. We found that these pressure values
    grew to large for the program to handle. The source for these equations does not explain them very clearly,
    so we are not sure what we are missing.

    def calculate_Pressure (self):
        x_Pressure_Total = 0;
        y_Pressure_Total = 0;
        z_Pressure_Total = 0;
        for molecule in self.molecules:
            for atom in molecule.atoms:
                x_Pressure_Total += x_Pressure_Total + atom.pressure[X]
                y_Pressure_Total += y_Pressure_Total + atom.pressure[Y]
                z_Pressure_Total += z_Pressure_Total + atom.pressure[Z]
        self.pressure_Total = (((x_Pressure_Total + y_Pressure_Total) / 2 ) -  z_Pressure_Total) / (self.x_len * self.y_len)
    '''