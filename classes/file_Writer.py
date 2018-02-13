##########################################################################
# This is a class containing functions for parsing a data file into usable
# data classes
# openMM
# Author: Chris McCormick
# Date: Summer 2016
##########################################################################

import csv
from numpy import abs

class file_Writer(object):
    def write_Interfaces(self, my_Sim_Data, my_Input_File_Path):
        sim_Data = my_Sim_Data
        input_File_Path = my_Input_File_Path
        with open(input_File_Path + '.interfaces.txt', 'wb') as interfaces_csv:
            interface_Writer = csv.writer(interfaces_csv, delimiter = '\t')
            print_Array = []
            print_Array.append(["Frame", "Interface"])
            for frame in sim_Data.frames:
                print_Array.append([frame.frame_ID, frame.interface_Z_Pos])
            interface_Writer.writerows(print_Array)
    
    def write_Avg_Frame_Densities(self, my_All_Frames, my_Input_File_Path):
        all_Frames = my_All_Frames
        input_File_Path = my_Input_File_Path
        with open(input_File_Path + '.Density_Analysis.txt', 'wb') as AVG_Frame_Densities_csv:
            data_Writer = csv.writer(AVG_Frame_Densities_csv, delimiter = '\t')
            print_Array = []
            print_Row = []
            print_Row.append("SliceID")
            print_Row.append("Z-Position")
            for each_Type in all_Frames.slices_Above_Interface[0].atom_Types:
                print_Row.append(each_Type.type_Name + ' Density')
                print_Row.append(each_Type.type_Name + ' Error')
            for each_Group in all_Frames.slices_Above_Interface[0].atom_Groups:
                print_Row.append(each_Group.group_Name + ' Density')
                print_Row.append(each_Group.group_Name + ' Error')
            print_Array.append(print_Row)
            for each_Slice in all_Frames.slices_Above_Interface:
                print_Array.append(self.build_Density_Row(each_Slice))
            for each_Slice in all_Frames.slices_Below_Interface:
                print_Array.append(self.build_Density_Row(each_Slice))
            data_Writer.writerows(print_Array)
            
    def build_Density_Row(self, the_Slice):
        my_Slice = the_Slice
        Da_PER_nm3_CONVERSION = 1.66053e-3
        print_Row = []
        print_Row.append(my_Slice.slice_ID)
        print_Row.append(my_Slice.slice_ID * my_Slice.slice_Depth * 10) # output in angstroms
        for each_Type in my_Slice.atom_Types:
            print_Row.append(each_Type.avg_Slice_Density_In_Sim * Da_PER_nm3_CONVERSION)
            print_Row.append(each_Type.stdErr_Density * Da_PER_nm3_CONVERSION)
        for each_Group in my_Slice.atom_Groups:
            print_Row.append(each_Group.group_Density * Da_PER_nm3_CONVERSION)
            print_Row.append(each_Group.group_StErr * Da_PER_nm3_CONVERSION)
        return print_Row
    
    def write_Surface_Heights(self, surface, input_File_Path, surface_Name):
        # the write method for these is append - make sure you don't have an existing file with the same name (you shouldn't generally because 
        # file names are generated with time stamps)
        with open(input_File_Path + '.' + surface_Name + '.surface_Height.txt', 'a') as surface_Height_csv:
            surface_Height_Writer = csv.writer(surface_Height_csv, delimiter = '\t')
            print_Surface = [[0 for i in range(40)] for j in range(40)]
            for i in range(40):
                for j in range(40):
                    print_Surface[i][j] = surface[i][j].z_Intersect
                    if print_Surface[i][j] < 0 and print_Surface[i][j] is not None:
                        print_Surface[i][j] = print_Surface[i][j] * -1
            surface_Height_Writer.writerows(print_Surface)
            # this inserts a separator blank row between each frame of description.
            surface_Height_Writer.writerows(" ")
        
    def write_Surface_Types(self, surface, input_File_Path, surface_Name):
        # the write method for these is append - make sure you don't have an existing file with the same name (you shouldn't generally because 
        # file names are generated with time stamps)
        with open(input_File_Path + '.' + surface_Name + '.surface_Type.txt', 'a') as surface_Type_csv:
            surface_Type_Writer = csv.writer(surface_Type_csv, delimiter = '\t')
            print_Surface = [[0 for i in range(40)] for j in range(40)]
            for i in range(40):
                for j in range(40):
                    print_Surface[i][j] = surface[i][j].atom_Type
            surface_Type_Writer.writerows(print_Surface)
            # this inserts a separator blank row between each frame of description.
            surface_Type_Writer.writerows(" ")