##########################################################################
# This is a class definition for analysis of groups of frames
# Currently, the only group of frames we are interested in is the group of
# all frames, for averaging accross the durration of a simulation
# Author: Chris McCormick
# Date: Summer 2016
##########################################################################
from classes.z_Slice import *
from copy import deepcopy

class frame_Group(object):
    def __init__(self):
        self.slices_Above_Interface = []    # These slices store data that is aggregated accross the frame group
        self.slices_Below_Interface = []    # These slices store data that is aggregated accross the frame group
        
    # This is the ugliest function in the program. It goes through each frame of simulation, and adds it's data
    # to the aggregate slices for the frame group.
    def compute_Frame_Group_Data(self, my_Simulation, my_Atom_Groups):
        simulation = my_Simulation
        atom_Groups = my_Atom_Groups

        for each_Frame in simulation.frames:
            for each_Slice in each_Frame.slices_Above_Interface:
                if not self.slices_Above_Interface:
                    self.slices_Above_Interface.append(z_Slice(each_Slice.slice_ID, None, each_Frame.x_len, each_Frame.y_len, simulation.atom_Type_List, atom_Groups))
                    self.slices_Above_Interface[-1].atom_Types = deepcopy(simulation.atom_Type_List.atom_Types)
                    for i in range(len(simulation.atom_Type_List.atom_Types)):
                        self.slices_Above_Interface[-1].atom_Types[i].all_Densities_In_Slice.append(each_Slice.get_Density_Or_Zero(simulation.atom_Type_List.atom_Types[i].type_Name))
                else:
                    if each_Slice.slice_ID >= len(self.slices_Above_Interface):
                        self.slices_Above_Interface.append(z_Slice(each_Slice.slice_ID, None, each_Frame.x_len, each_Frame.y_len, simulation.atom_Type_List, atom_Groups))
                        self.slices_Above_Interface[-1].atom_Types = deepcopy(simulation.atom_Type_List.atom_Types)
                        for i in range(len(simulation.atom_Type_List.atom_Types)):
                            self.slices_Above_Interface[-1].atom_Types[i].all_Densities_In_Slice.append(each_Slice.get_Density_Or_Zero(simulation.atom_Type_List.atom_Types[i].type_Name))
                    else:
                        for i in range(len(simulation.atom_Type_List.atom_Types)):
                            self.slices_Above_Interface[each_Slice.slice_ID].atom_Types[i].all_Densities_In_Slice.append(each_Slice.get_Density_Or_Zero(simulation.atom_Type_List.atom_Types[i].type_Name))
            
            for each_Slice in each_Frame.slices_Below_Interface:
                if not self.slices_Below_Interface:
                    self.slices_Below_Interface.append(z_Slice(each_Slice.slice_ID, None, each_Frame.x_len, each_Frame.y_len, simulation.atom_Type_List, atom_Groups))
                    self.slices_Below_Interface[-1].atom_Types = deepcopy(simulation.atom_Type_List.atom_Types)
                    for i in range(len(simulation.atom_Type_List.atom_Types)):
                        self.slices_Below_Interface[-1].atom_Types[i].all_Densities_In_Slice.append(each_Slice.get_Density_Or_Zero(simulation.atom_Type_List.atom_Types[i].type_Name))
                else:
                    if ((each_Slice.slice_ID * -1) - 1) >= len(self.slices_Below_Interface):
                        self.slices_Below_Interface.append(z_Slice(each_Slice.slice_ID, None, each_Frame.x_len, each_Frame.y_len, simulation.atom_Type_List, atom_Groups))
                        self.slices_Below_Interface[-1].atom_Types = deepcopy(simulation.atom_Type_List.atom_Types)
                        for i in range(len(simulation.atom_Type_List.atom_Types)):
                            self.slices_Below_Interface[-1].atom_Types[i].all_Densities_In_Slice.append(each_Slice.get_Density_Or_Zero(simulation.atom_Type_List.atom_Types[i].type_Name))
                    else:
                        for i in range(len(simulation.atom_Type_List.atom_Types)):
                            self.slices_Below_Interface[(each_Slice.slice_ID * -1) - 1].atom_Types[i].all_Densities_In_Slice.append(each_Slice.get_Density_Or_Zero(simulation.atom_Type_List.atom_Types[i].type_Name))