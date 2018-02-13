##########################################################################
# This is a class definition for a z_Slice
# Author: Chris McCormick
# Date: Summer 2016
##########################################################################
from copy import deepcopy
from classes.atom_Type import atom_Type
from global_Vars import *

class z_Slice(object):
    slice_Depth    = SLICE_SIZE # defined in global_Vars.py
    
    def __init__(self, slice_ID, bottom, x_Len, y_Len, type_List, atom_Groups):
        self.slice_ID = slice_ID
        self.bottom_Bound   = bottom # bottom bound of slice in nm
        self.volume         = self.slice_Depth * x_Len * y_Len # nm Cubed
        self.atom_Types = []
        self.atom_Groups = deepcopy(atom_Groups)    # each slice gets it's own list of atom groups 
        
    # goes through all of the atoms in a frame, and builds a list of atom types that exist in that slice.
    # As it does this, it also keeps a running total of the mass of that atom type in that slice.
    def load_Slice_Masses(self, my_Frame):
        frame = my_Frame
        
        for molecule in frame.molecules:
            for atom in molecule.atoms:
                if atom.position[Z] >= self.bottom_Bound and atom.position[Z] < self.bottom_Bound + self.slice_Depth:
                    add_Me = True
                    if not self.atom_Types:
                        new_A_Type = atom_Type()
                        new_A_Type.molecule_Name = molecule.molecule_Name
                        new_A_Type.type_Name = molecule.molecule_Name
                        new_A_Type.atom_Name = atom.atom_Name
                        new_A_Type.type_Name = new_A_Type.type_Name + '-' + atom.atom_Name
                        self.atom_Types.append(new_A_Type)
                        self.atom_Types[-1].tot_Mass_In_Slice = self.atom_Types[-1].tot_Mass_In_Slice + atom.atom_Mass
                        add_Me = False
                    else:
                        for i in range(len(self.atom_Types)):
                            if molecule.molecule_Name == self.atom_Types[i].molecule_Name and atom.atom_Name == self.atom_Types[i].atom_Name:
                                self.atom_Types[i].tot_Mass_In_Slice = self.atom_Types[i].tot_Mass_In_Slice + atom.atom_Mass
                                add_Me = False
                        if add_Me == True:
                            new_A_Type = atom_Type()
                            new_A_Type.molecule_Name = molecule.molecule_Name
                            new_A_Type.type_Name = molecule.molecule_Name
                            new_A_Type.atom_Name = atom.atom_Name
                            new_A_Type.type_Name = new_A_Type.type_Name + '-' + atom.atom_Name
                            self.atom_Types.append(new_A_Type)
                            self.atom_Types[-1].tot_Mass_In_Slice = self.atom_Types[-1].tot_Mass_In_Slice + atom.atom_Mass

    # calls the calculate density function for each atom type that exists in a given slice
    def calculate_Densities(self):
        for a_Type in self.atom_Types:
            a_Type.calculate_Density(self.volume)

    # returns the density in this slice of the atom type passed in. If the atom type does not exist in this
    # slice, return 0.
    def get_Density_Or_Zero(self, check_Type):
        check_Type_Name = check_Type
        if not self.atom_Types:
            return 0
        for each_Type in self.atom_Types:
            if check_Type_Name == each_Type.type_Name:
                return each_Type.density_In_Slice
        return 0