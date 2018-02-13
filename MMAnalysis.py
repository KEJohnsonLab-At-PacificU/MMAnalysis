#!/usr/local/bin/PyPy

##########################################################################
# This is a driver for analysis of data produced through simulations using
# openMM
# Author: Chris McCormick
# Date: Summer 2016
##########################################################################

#we're actually doing tab delimited text, but csv is useful
import sys
import csv
import driver_Functions as driver
from global_Vars import *
from classes.file_Reader import file_Reader
from classes.file_Writer import file_Writer
from classes.atom_Type_Group import build_Group_List
from classes.frame_Group import frame_Group
from datetime import datetime
from classes.surface_Point import surface_Point
from numpy import sqrt, square

time_Start = 0                              # Holds the start time of some process
time_End = 0                                # Holds the end time of some process
tot_Time_Start = datetime.now()             # Holds the start time of some process
tot_Time_End = datetime.now()               # Holds the end time of some process
time_Difference = time_End - time_Start

print("##########################################################################")
print("#                             MMAnalysis.py                              #")
print("##########################################################################")

atom_Groups = []

input_File_Path = "input.txt"       # default file if no argument is given
file_Reader = file_Reader()
file_Writer = file_Writer()

if len(sys.argv) > 2:
    print ("This program accepts one argument, which is the file path to begin the analysis. If no argument is provided, the default 'input.txt' is used")
    sys.exit()
if len(sys.argv) == 2:
    input_File_Path = str(sys.argv[1])

print("Loading Simulation Data")
time_Start = datetime.now()
sim_Data = file_Reader.load_Sim_Data(input_File_Path)
time_End = datetime.now()
time_Difference = time_End - time_Start
print ('\tRan for ' + str(time_Difference))

print("Calculating Interfaces")
time_Start = datetime.now()
for frame in sim_Data.frames:
    frame.get_Interface_Z()
time_End = datetime.now()
time_Difference = time_End - time_Start
print ('\tRan for ' + str(time_Difference))

file_Writer.write_Interfaces(sim_Data, input_File_Path)

print("Building Group List")
time_Start = datetime.now()
build_Group_List(atom_Groups, sim_Data.frames[0].molecules)
time_End = datetime.now()
time_Difference = time_End - time_Start
print ('\tRan for ' + str(time_Difference))

print("Building Type List")
time_Start = datetime.now()
sim_Data.atom_Type_List.build_Type_List(sim_Data.frames[0])
time_End = datetime.now()
time_Difference = time_End - time_Start
print ('\tRan for ' + str(time_Difference))
    
print("Building Slices")
time_Start = datetime.now()
for each_Frame in sim_Data.frames:
    each_Frame.build_Slices(sim_Data.atom_Type_List, atom_Groups)
time_End = datetime.now()
time_Difference = time_End - time_Start
print ('\tRan for ' + str(time_Difference))

print("Calculating Densities")
time_Start = datetime.now()
for each_Frame in sim_Data.frames:
    for z_Slice in each_Frame.slices_Above_Interface:
        z_Slice.load_Slice_Masses(each_Frame)
        z_Slice.calculate_Densities()
    for z_Slice in each_Frame.slices_Below_Interface:
        z_Slice.load_Slice_Masses(each_Frame)
        z_Slice.calculate_Densities()
time_End = datetime.now()
time_Difference = time_End - time_Start
print ('\tRan for ' + str(time_Difference))

print("Making Frame Group")
time_Start = datetime.now()
all_Frames = frame_Group()
time_End = datetime.now()
time_Difference = time_End - time_Start
print ('\tRan for ' + str(time_Difference))

print("Summing Frame Group Data")
time_Start = datetime.now()
all_Frames.compute_Frame_Group_Data(sim_Data, atom_Groups)
time_End = datetime.now()
time_Difference = time_End - time_Start
print ('\tRan for ' + str(time_Difference))

print("Computing AVG Frame Group Data")
time_Start = datetime.now()
for each_Slice in all_Frames.slices_Above_Interface:
    for each_Type in each_Slice.atom_Types:
        each_Type.calculate_Avg_Slice_Density()
        each_Type.calculate_stdErr_Density()
for each_Slice in all_Frames.slices_Below_Interface:
    for each_Type in each_Slice.atom_Types:
        each_Type.calculate_Avg_Slice_Density()
        each_Type.calculate_stdErr_Density()

time_End = datetime.now()
time_Difference = time_End - time_Start
print ('\tRan for ' + str(time_Difference))

print("summing Densities for groups")
time_Start = datetime.now()
for each_Slice in all_Frames.slices_Above_Interface:
    for group in each_Slice.atom_Groups:
        for each_Group_Definition in atom_Groups: 
            if group.group_Name == each_Group_Definition.group_Name:
                group.sum_Density(each_Slice)
                group.calc_Group_StErr(each_Slice)
for each_Slice in all_Frames.slices_Below_Interface:
    for group in each_Slice.atom_Groups:
        for each_Group_Definition in atom_Groups: 
            if group.group_Name == each_Group_Definition.group_Name:
                group.sum_Density(each_Slice)
                group.calc_Group_StErr(each_Slice)
time_End = datetime.now()
time_Difference = time_End - time_Start
print ('\tRan for ' + str(time_Difference))

file_Writer.write_Avg_Frame_Densities(all_Frames, input_File_Path)

'''
This is our aborted attempt to calculate interfacial pressure.
# Pressure calc #
print ("let's try calculating interface properties?")
for each_Molecule in sim_Data.frames[0].molecules:
    for each_Atom in each_Molecule.atoms:
        each_Atom.calc_Pressure()
sim_Data.frames[0].calculate_Pressure()
print ("Pressure for frame 0 calculated as :" + str(sim_Data.frames[0].pressure_Total))
'''
print("describing the surface...")
time_Start = datetime.now()
# here we define a water surface. This represents an xy grid, where point is 0.1nm apart. 0.1nm was chosen because 
# the van der waals of hydrogen is just over that, so we won't miss any atoms when calculating the VDW intercepts.

# This whole mess should be it's own function, but I'm trying to wrap this up before the semester starts

# for each frame
for each_Frame in sim_Data.frames:
    water_Surface = [[surface_Point() for j in range(int((sim_Data.frames[0].x_len * 10)))] for i in range(int(sim_Data.frames[0].y_len * 10))]
    # Look at each atom
    for each_Molecule in each_Frame.molecules:
        # being in a specific atom group (now water or organic)
        if each_Molecule.molecule_Name == "HOH":
            for each_Atom in each_Molecule.atoms:
                # decide which points in surface this atom is candidate for
                if each_Atom.position[Z] > SURF_BOTTOM_TOLERANCE or each_Atom.position[Z] < SURF_TOP_TOLERANCE:
                    # We check the 16 nearest points with root(dX^2 + dY^2) < VDW
                    hits = []
                    if "H" in each_Atom.atom_Name:
                        each_Atom.VDW = H_VDW
                    if each_Atom.atom_Name == "O":
                        each_Atom.VDW = O_VDW
                    probe_X = int((each_Atom.position[X] * 10) - 2)
                    probe_Y = int((each_Atom.position[Y] * 10) - 2)
                    for i in range(4):
                        for j in range(4):
                            if sqrt(square(((each_Atom.position[X] * 10) - (probe_X + i))) + square(((each_Atom.position[Y] * 10) - (probe_Y + j)))) < each_Atom.VDW:
                                hits.append([probe_X + i, probe_Y + j])
                                intersect_X_Probe = probe_X + i
                                intersect_Y_Probe = probe_Y + j
                                # Here we are checking whether candidate points are across the periodic boundary
                                if probe_X + i < 0 or probe_Y + j < 0 or probe_X + i >= 40 or probe_Y + j >= 40:
                                    #print ("I wrapped!")
                                    if probe_X + i < 0:
                                        intersect_X_Probe = 40 + (probe_X + i)
                                    if probe_Y + j < 0:
                                        intersect_Y_Probe = 40 + (probe_Y + j)
                                        
                                    if probe_X + i >= 40:
                                        intersect_X_Probe = (probe_X + i) - 40
                                    if probe_Y + j >= 40:
                                        intersect_Y_Probe = (probe_Y + j) - 40
                                        
                                # for each point, compare intersect height with that point in matrix of surface description
                                # if max (for water), or min (for org) save atom type and intersect Z value in matrix
                                if water_Surface[intersect_X_Probe][intersect_Y_Probe].z_Intersect < water_Surface[intersect_X_Probe][intersect_Y_Probe].calc_Intersect(each_Atom, probe_X + i, probe_Y + j):
                                    water_Surface[intersect_X_Probe][intersect_Y_Probe].z_Intersect = water_Surface[intersect_X_Probe][intersect_Y_Probe].calc_Intersect(each_Atom, probe_X + i, probe_Y + j)
                                    water_Surface[intersect_X_Probe][intersect_Y_Probe].atom_Type = each_Atom.atom_Name
    # write out the heights and types for the frame
    file_Writer.write_Surface_Heights(water_Surface, input_File_Path, "HOH")
    file_Writer.write_Surface_Types(water_Surface, input_File_Path, "HOH")
                        
# I have an idea for doing the bottom surface of the organic with minimal modifications. We can just change the molecule name we're checking for above
# and then, for finding the bottom surface, just invert the z dimension of the simulation and return the absolute value result.

# We will also need the Van Der Waals for the atoms in the organic, which I don't have access to at the moment.
# for each frame
for each_Frame in sim_Data.frames:
    org_Surface = [[surface_Point() for j in range(int((sim_Data.frames[0].x_len * 10)))] for i in range(int(sim_Data.frames[0].y_len * 10))]
    # Look at each atom
    for each_Molecule in each_Frame.molecules:
        # being in a specific atom group (now water or organic)
        if each_Molecule.molecule_Name == "HXN" or each_Molecule.molecule_Name == "CTC":
            for each_Atom in each_Molecule.atoms:
                # decide which points in surface this atom is candidate for
                if each_Atom.position[Z] > SURF_BOTTOM_TOLERANCE or each_Atom.position[Z] < SURF_TOP_TOLERANCE:
                    # We check the 16 nearest points with root(dX^2 + dY^2) < VDW
                    hits = []
                    if "H" in each_Atom.atom_Name:
                        each_Atom.VDW = H_VDW
                    if "C" in each_Atom.atom_Name:
                        each_Atom.VDW = C_VDW
                    if "Cl" in each_Atom.atom_Name:
                        each_Atom.VDW = CL_VDW
                    probe_X = int((each_Atom.position[X] * 10) - 2)
                    probe_Y = int((each_Atom.position[Y] * 10) - 2)
                    for i in range(4):
                        for j in range(4):
                            if sqrt(square(((each_Atom.position[X] * 10) - (probe_X + i))) + square(((each_Atom.position[Y] * 10) - (probe_Y + j)))) < each_Atom.VDW:
                                hits.append([probe_X + i, probe_Y + j])
                                intersect_X_Probe = probe_X + i
                                intersect_Y_Probe = probe_Y + j
                                # Here we are checking whether candidate points are across the periodic boundary
                                if probe_X + i < 0 or probe_Y + j < 0 or probe_X + i >= 40 or probe_Y + j >= 40:
                                    #print ("I wrapped!")
                                    if probe_X + i < 0:
                                        intersect_X_Probe = 40 + (probe_X + i)
                                    if probe_Y + j < 0:
                                        intersect_Y_Probe = 40 + (probe_Y + j)
                                        
                                    if probe_X + i >= 40:
                                        intersect_X_Probe = (probe_X + i) - 40
                                    if probe_Y + j >= 40:
                                        intersect_Y_Probe = (probe_Y + j) - 40
                                        
                                # for each point, compare intersect height with that point in matrix of surface description
                                # if max (for water), or min (for org) save atom type and intersect Z value in matrix
                                if org_Surface[intersect_X_Probe][intersect_Y_Probe].z_Intersect < org_Surface[intersect_X_Probe][intersect_Y_Probe].calc_Intersect_Inverse(each_Atom, probe_X + i, probe_Y + j):
                                    org_Surface[intersect_X_Probe][intersect_Y_Probe].z_Intersect = org_Surface[intersect_X_Probe][intersect_Y_Probe].calc_Intersect_Inverse(each_Atom, probe_X + i, probe_Y + j)
                                    org_Surface[intersect_X_Probe][intersect_Y_Probe].atom_Type = each_Atom.atom_Name
    # write out the heights and types for the frame
    file_Writer.write_Surface_Heights(org_Surface, input_File_Path, "ORG")
    file_Writer.write_Surface_Types(org_Surface, input_File_Path, "ORG")


time_End = datetime.now()
time_Difference = time_End - time_Start
print ('\tRan for ' + str(time_Difference))


tot_Time_End = datetime.now()
time_Difference = tot_Time_End - tot_Time_Start
print ('total time: ' + str(time_Difference))
print("##########################################################################")
print("#                          Analysis Completed                            #")
print("##########################################################################")
