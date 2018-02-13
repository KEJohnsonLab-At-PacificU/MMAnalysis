##########################################################################
# This is a class definition for atom_Type. An atom type is
# a unique combination of molecule name and atom name, i.e. all H1 atoms belonging to any 
# water(HOH) molecule 
# Author: Chris McCormick
# Date: Summer 2016
##########################################################################

import numpy as np
from cmath import sqrt

class atom_Type(object):

    def __init__(self):
        '''
        commented out are attributes that are now only assigned when needed within other function calls.
        again, this is because I wanted to see the impact on memory usage, and it seems to have helped a bit.
        I keep the comments here just as a reminded of what attributes each atom Type can have
        '''
        #self.type_Name = None
        #self.atom_Name = None
        #self.molecule_Name = None
        self.tot_Mass_In_Slice = 0 # daltons
        #self.density_In_Slice = None
        #self.avg_Slice_Density_In_Sim = None   # this is the average of density_In_Slice across all frames of a simulation
        #self.stdErr_Density = None
        self.all_Densities_In_Slice = []    # this contains an array of the densities of this atom type in this slice, across all frames 

    def calculate_Density(self, my_Volume):
        volume = my_Volume
        self.density_In_Slice = self.tot_Mass_In_Slice / volume   # Da / nm^3
        
    def calculate_Avg_Slice_Density(self):
        self.avg_Slice_Density_In_Sim = np.average(self.all_Densities_In_Slice)
        
    # calculates the standard error based on the standard deviation.
    def calculate_stdErr_Density(self):
        self.stdErr_Density = np.std(self.all_Densities_In_Slice) / abs(sqrt(len(self.all_Densities_In_Slice)))