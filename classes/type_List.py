##########################################################################
# This is a class definition for a list of atom_Types
# Author: Chris McCormick
# Date: Summer 2016
##########################################################################
from classes.atom_Type import atom_Type

class type_List(object):
    def __init__(self):
        self.atom_Types = []
            
    # method returns true if the atom passed in exists in the list of atom_Types
    def is_In_Type_List (self, my_Compare_Atom):
        compare_Atom = my_Compare_Atom

        if not self.atom_Types:
            return False
        for a_Type in self.atom_Types:
            if a_Type.type_Name == compare_Atom.type_Name:
                return True
        return False
    
    # constructs the list of all atom types that exist in a simulation based on a single frame 
    # (since this is a closed simulation, we only have to look at one frame). An atom type is
    # a unique combination of molecule name and atom name, i.e. all H1 atoms belonging to any 
    # water(HOH) molecule 
    def build_Type_List(self, my_Frame):
        frame = my_Frame

        for molec in frame.molecules:
            for atom in molec.atoms:
                new_A_Type = atom_Type()
                new_A_Type.molecule_Name = molec.molecule_Name
                new_A_Type.type_Name = molec.molecule_Name
                new_A_Type.atom_Name = atom.atom_Name
                new_A_Type.type_Name = new_A_Type.type_Name + '-' + atom.atom_Name
                if not self.is_In_Type_List(new_A_Type):
                    self.atom_Types.append(new_A_Type)