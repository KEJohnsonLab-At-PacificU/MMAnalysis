##########################################################################
# This is a class definition for analysis of groups of atom_Types
# This also holds the functions for building the group list at the bottom
# I've gone back and forth about putting the definitions in a separate file
# Author: Chris McCormick
# Date: Summer 2016
##########################################################################

from cmath import sqrt
from numpy import std

class atom_Type_Group(object):

    def __init__(self, group_Name):
        self.group_Name = group_Name
        self.types_In_Group = []    # this is a list of all atom types in a particular group of atom types
        self.group_Density = 0      # in daltons / nanometer^2
        self.group_StErr = 0
        
    # This function goes through all of the atoms in a given slice, for each atom that belongs to this group of
    # atoms, it adds it's density in the slice to the running sum of densities for all atoms in this group.
    def sum_Density(self, my_Z_Slice):
        z_Slice = my_Z_Slice

        for a_Type in z_Slice.atom_Types:
            for group_Types in self.types_In_Group:
                if a_Type.type_Name == group_Types:
                    self.group_Density = self.group_Density + a_Type.avg_Slice_Density_In_Sim
    
    # calculates the combined standard error of all atom types belonging to the group.
    def calc_Group_StErr(self, my_Z_Slice):
        z_Slice = my_Z_Slice

        std_Err_Divisor = 0
        for a_Type in z_Slice.atom_Types:   # for every atom type in the slice
            for group_Types in self.types_In_Group: # for every atom type in the atom_Type_Group
                if a_Type.type_Name == group_Types: # if the atom type is in the type group
                    type_STD = std(a_Type.all_Densities_In_Slice)   # calculate the standard error of it's density values
                    self.group_StErr = self.group_StErr + (type_STD * type_STD) # and add to the combined standard error
                    std_Err_Divisor += len(a_Type.all_Densities_In_Slice)       # for the group
        if not std_Err_Divisor == 0:    # make sure we don't divide by zero. If it is, error value remains 0 because there is nothing to find the standard error of
            self.group_StErr = abs(sqrt(self.group_StErr)) / abs(sqrt(std_Err_Divisor))


###########################################################################
#               Definition of possible atom_Type_Groups
# This is not technically part of the class, but I decided it was the best
# place to define the possible atom type groups. atom_Type_Groups are
# molecules, or parts of molecules we are interested in, such as 
# 'head group' or 'tail'
###########################################################################

# Here at the top, we define each group. Below, we have a function which calls these to build the group list
# depending on which are present in the simulation.

def define_Water():
    water_Group = atom_Type_Group("water")
    water_Group.types_In_Group.append("HOH-H1")
    water_Group.types_In_Group.append("HOH-H2")
    water_Group.types_In_Group.append("HOH-O")
    return water_Group

def define_CTC():
    ctc_Group = atom_Type_Group("ctc")
    ctc_Group.types_In_Group.append("CTC-C")
    ctc_Group.types_In_Group.append("CTC-Cl1")
    ctc_Group.types_In_Group.append("CTC-Cl2")
    ctc_Group.types_In_Group.append("CTC-Cl3")
    ctc_Group.types_In_Group.append("CTC-Cl4")
    return ctc_Group

def define_HXN():
    hexane_Group = atom_Type_Group("hexane")
    hexane_Group.types_In_Group.append("HXN-C1")
    hexane_Group.types_In_Group.append("HXN-H01")
    hexane_Group.types_In_Group.append("HXN-H02")
    hexane_Group.types_In_Group.append("HXN-H03")
    hexane_Group.types_In_Group.append("HXN-C2")
    hexane_Group.types_In_Group.append("HXN-H04")
    hexane_Group.types_In_Group.append("HXN-H05")
    hexane_Group.types_In_Group.append("HXN-C3")
    hexane_Group.types_In_Group.append("HXN-H06")
    hexane_Group.types_In_Group.append("HXN-H07")
    hexane_Group.types_In_Group.append("HXN-C4")
    hexane_Group.types_In_Group.append("HXN-H08")
    hexane_Group.types_In_Group.append("HXN-H09")
    hexane_Group.types_In_Group.append("HXN-C5")
    hexane_Group.types_In_Group.append("HXN-H10")
    hexane_Group.types_In_Group.append("HXN-H11")
    hexane_Group.types_In_Group.append("HXN-C6")
    hexane_Group.types_In_Group.append("HXN-H12")
    hexane_Group.types_In_Group.append("HXN-H13")
    hexane_Group.types_In_Group.append("HXN-H14")
    return hexane_Group
    
def define_DOD():
    dod_Group = atom_Type_Group("dod")
    dod_Group.types_In_Group.append("dod-C1")
    dod_Group.types_In_Group.append("dod-H1")
    dod_Group.types_In_Group.append("dod-H2")
    dod_Group.types_In_Group.append("dod-H3")
    dod_Group.types_In_Group.append("dod-C2")
    dod_Group.types_In_Group.append("dod-H4")
    dod_Group.types_In_Group.append("dod-H5")
    dod_Group.types_In_Group.append("dod-C3")
    dod_Group.types_In_Group.append("dod-H6")
    dod_Group.types_In_Group.append("dod-H7")
    dod_Group.types_In_Group.append("dod-C4")
    dod_Group.types_In_Group.append("dod-H8")
    dod_Group.types_In_Group.append("dod-H9")
    dod_Group.types_In_Group.append("dod-C5")
    dod_Group.types_In_Group.append("dod-H10")
    dod_Group.types_In_Group.append("dod-H11")
    dod_Group.types_In_Group.append("dod-C6")
    dod_Group.types_In_Group.append("dod-H12")
    dod_Group.types_In_Group.append("dod-H13")
    dod_Group.types_In_Group.append("dod-C7")
    dod_Group.types_In_Group.append("dod-H14")
    dod_Group.types_In_Group.append("dod-H15")
    dod_Group.types_In_Group.append("dod-C8")
    dod_Group.types_In_Group.append("dod-H16")
    dod_Group.types_In_Group.append("dod-H17")
    dod_Group.types_In_Group.append("dod-C9")
    dod_Group.types_In_Group.append("dod-H18")
    dod_Group.types_In_Group.append("dod-H19")
    dod_Group.types_In_Group.append("dod-C10")
    dod_Group.types_In_Group.append("dod-H20")
    dod_Group.types_In_Group.append("dod-H21")
    dod_Group.types_In_Group.append("dod-C11")
    dod_Group.types_In_Group.append("dod-H22")
    dod_Group.types_In_Group.append("dod-H23")
    dod_Group.types_In_Group.append("dod-CO")
    dod_Group.types_In_Group.append("dod-OC")
    dod_Group.types_In_Group.append("dod-OH")
    dod_Group.types_In_Group.append("dod-HO")
    return dod_Group
    
def define_DOD_Head():
    dodHead_Group = atom_Type_Group("dodHead")
    dodHead_Group.types_In_Group.append("dod-CO")
    dodHead_Group.types_In_Group.append("dod-OC")
    dodHead_Group.types_In_Group.append("dod-OH")
    dodHead_Group.types_In_Group.append("dod-HO")
    return dodHead_Group
    
def define_DOD_Tail():
    dodTail_Group = atom_Type_Group("dodTail")
    dodTail_Group.types_In_Group.append("dod-C1")
    dodTail_Group.types_In_Group.append("dod-H1")
    dodTail_Group.types_In_Group.append("dod-H2")
    dodTail_Group.types_In_Group.append("dod-H3")
    dodTail_Group.types_In_Group.append("dod-C2")
    dodTail_Group.types_In_Group.append("dod-H4")
    dodTail_Group.types_In_Group.append("dod-H5")
    dodTail_Group.types_In_Group.append("dod-C3")
    dodTail_Group.types_In_Group.append("dod-H6")
    dodTail_Group.types_In_Group.append("dod-H7")
    dodTail_Group.types_In_Group.append("dod-C4")
    dodTail_Group.types_In_Group.append("dod-H8")
    dodTail_Group.types_In_Group.append("dod-H9")
    dodTail_Group.types_In_Group.append("dod-C5")
    dodTail_Group.types_In_Group.append("dod-H10")
    dodTail_Group.types_In_Group.append("dod-H11")
    dodTail_Group.types_In_Group.append("dod-C6")
    dodTail_Group.types_In_Group.append("dod-H12")
    dodTail_Group.types_In_Group.append("dod-H13")
    dodTail_Group.types_In_Group.append("dod-C7")
    dodTail_Group.types_In_Group.append("dod-H14")
    dodTail_Group.types_In_Group.append("dod-H15")
    dodTail_Group.types_In_Group.append("dod-C8")
    dodTail_Group.types_In_Group.append("dod-H16")
    dodTail_Group.types_In_Group.append("dod-H17")
    dodTail_Group.types_In_Group.append("dod-C9")
    dodTail_Group.types_In_Group.append("dod-H18")
    dodTail_Group.types_In_Group.append("dod-H19")
    dodTail_Group.types_In_Group.append("dod-C10")
    dodTail_Group.types_In_Group.append("dod-H20")
    dodTail_Group.types_In_Group.append("dod-H21")
    dodTail_Group.types_In_Group.append("dod-C11")
    dodTail_Group.types_In_Group.append("dod-H22")
    dodTail_Group.types_In_Group.append("dod-H23")
    return dodTail_Group

def define_DDA():
    dda_Group = atom_Type_Group("dda")
    dda_Group.types_In_Group.append("dda-C1")
    dda_Group.types_In_Group.append("dda-H1")
    dda_Group.types_In_Group.append("dda-H2")
    dda_Group.types_In_Group.append("dda-H3")
    dda_Group.types_In_Group.append("dda-C2")
    dda_Group.types_In_Group.append("dda-H4")
    dda_Group.types_In_Group.append("dda-H5")
    dda_Group.types_In_Group.append("dda-C3")
    dda_Group.types_In_Group.append("dda-H6")
    dda_Group.types_In_Group.append("dda-H7")
    dda_Group.types_In_Group.append("dda-C4")
    dda_Group.types_In_Group.append("dda-H8")
    dda_Group.types_In_Group.append("dda-H9")
    dda_Group.types_In_Group.append("dda-C5")
    dda_Group.types_In_Group.append("dda-H10")
    dda_Group.types_In_Group.append("dda-H11")
    dda_Group.types_In_Group.append("dda-C6")
    dda_Group.types_In_Group.append("dda-H12")
    dda_Group.types_In_Group.append("dda-H13")
    dda_Group.types_In_Group.append("dda-C7")
    dda_Group.types_In_Group.append("dda-H14")
    dda_Group.types_In_Group.append("dda-H15")
    dda_Group.types_In_Group.append("dda-C8")
    dda_Group.types_In_Group.append("dda-H16")
    dda_Group.types_In_Group.append("dda-H17")
    dda_Group.types_In_Group.append("dda-C9")
    dda_Group.types_In_Group.append("dda-H18")
    dda_Group.types_In_Group.append("dda-H19")
    dda_Group.types_In_Group.append("dda-C10")
    dda_Group.types_In_Group.append("dda-H20")
    dda_Group.types_In_Group.append("dda-H21")
    dda_Group.types_In_Group.append("dda-C11")
    dda_Group.types_In_Group.append("dda-H22")
    dda_Group.types_In_Group.append("dda-H23")
    dda_Group.types_In_Group.append("dda-CO")
    dda_Group.types_In_Group.append("dda-OC")
    dda_Group.types_In_Group.append("dda-OH")
    return dda_Group

def define_DDA_Head():
    ddaHead_Group = atom_Type_Group("ddaHead")
    ddaHead_Group.types_In_Group.append("dda-CO")
    ddaHead_Group.types_In_Group.append("dda-OC")
    ddaHead_Group.types_In_Group.append("dda-OH")
    return ddaHead_Group

def define_DDA_Tail():
    ddaTail_Group = atom_Type_Group("ddaTail")
    ddaTail_Group.types_In_Group.append("dda-C1")
    ddaTail_Group.types_In_Group.append("dda-H1")
    ddaTail_Group.types_In_Group.append("dda-H2")
    ddaTail_Group.types_In_Group.append("dda-H3")
    ddaTail_Group.types_In_Group.append("dda-C2")
    ddaTail_Group.types_In_Group.append("dda-H4")
    ddaTail_Group.types_In_Group.append("dda-H5")
    ddaTail_Group.types_In_Group.append("dda-C3")
    ddaTail_Group.types_In_Group.append("dda-H6")
    ddaTail_Group.types_In_Group.append("dda-H7")
    ddaTail_Group.types_In_Group.append("dda-C4")
    ddaTail_Group.types_In_Group.append("dda-H8")
    ddaTail_Group.types_In_Group.append("dda-H9")
    ddaTail_Group.types_In_Group.append("dda-C5")
    ddaTail_Group.types_In_Group.append("dda-H10")
    ddaTail_Group.types_In_Group.append("dda-H11")
    ddaTail_Group.types_In_Group.append("dda-C6")
    ddaTail_Group.types_In_Group.append("dda-H12")
    ddaTail_Group.types_In_Group.append("dda-H13")
    ddaTail_Group.types_In_Group.append("dda-C7")
    ddaTail_Group.types_In_Group.append("dda-H14")
    ddaTail_Group.types_In_Group.append("dda-H15")
    ddaTail_Group.types_In_Group.append("dda-C8")
    ddaTail_Group.types_In_Group.append("dda-H16")
    ddaTail_Group.types_In_Group.append("dda-H17")
    ddaTail_Group.types_In_Group.append("dda-C9")
    ddaTail_Group.types_In_Group.append("dda-H18")
    ddaTail_Group.types_In_Group.append("dda-H19")
    ddaTail_Group.types_In_Group.append("dda-C10")
    ddaTail_Group.types_In_Group.append("dda-H20")
    ddaTail_Group.types_In_Group.append("dda-H21")
    ddaTail_Group.types_In_Group.append("dda-C11")
    ddaTail_Group.types_In_Group.append("dda-H22")
    ddaTail_Group.types_In_Group.append("dda-H23")
    return ddaTail_Group

# a function builds the group list depending on which things are present in the simulation.
def build_Group_List(atom_Groups, molecules):
    addWater = False
    addCTC = False
    addHXN = False
    addDOD = False
    addDDA = False
    for each_Molecule in molecules:
        if each_Molecule.molecule_Name == "HOH":
            addWater = True
        if each_Molecule.molecule_Name == "CTC":
            addCTC = True
        if each_Molecule.molecule_Name == "HXN":
            addHXN = True
        if each_Molecule.molecule_Name == "dod":
            addDOD = True
        if each_Molecule.molecule_Name == "dda":
            addDDA = True
    
    if addWater:
        water_Group = define_Water()
        atom_Groups.append(water_Group)
    if addCTC:
        ctc_Group = define_CTC()
        atom_Groups.append(ctc_Group)
    if addHXN:
        hexane_Group = define_HXN()
        atom_Groups.append(hexane_Group)
    if addDOD:
        dod_Group = define_DOD()
        dodHead_Group = define_DOD_Head()
        dodTail_Group = define_DOD_Tail()
        atom_Groups.append(dod_Group)
        atom_Groups.append(dodHead_Group)
        atom_Groups.append(dodTail_Group)
    if addDDA:
        dda_Group = define_DDA()
        ddaHead_Group = define_DDA_Head()
        ddaTail_Group = define_DDA_Tail()
        atom_Groups.append(dda_Group)
        atom_Groups.append(ddaHead_Group)
        atom_Groups.append(ddaTail_Group)