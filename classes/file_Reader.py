##########################################################################
# This is a class containing functions for parsing a data file into usable
# data classes
# Author: Chris McCormick
# Date: Summer 2016
##########################################################################

import csv
from classes.simulation_Data import simulation_Data
from classes.frame_Data import frame_Data
from classes.molecule_Data import molecule_Data
from classes.atom_Data import atom_Data

class file_Reader(object):
    
    """
    Input:        path to the file to read
    Function:     saves data from the input file into data structures for use in analysis
                  Calls load_Frames_Data
    Return:       A simulation_Data object, with properties (including frames>molecules>atoms) populated by data from input files
    """
    def load_Sim_Data(self, input_File_Path):
        # Constant definitions, these associate indexes for input files with usable labels.
        # Simulation Data Constants
        TOTAL_STEPS     = 0
        STEP_DURATION   = 1
        OUTPUT_INTERVAL = 2
        LOGFILE_PATH    = 3
        X_DIM           = 4
        Y_DIM           = 5
        Z_DIM           = 6
        
        # this input file has the simulation variables, including the file name of the datafile.
        with open(input_File_Path) as input_File:
            #skip the header line
            input_File.readline()
            csv_Input = csv.reader(input_File, delimiter = '\t')
            for row in csv_Input:
                sim_Data = simulation_Data(#int(row[TOTAL_STEPS]), 
                                           #row[STEP_DURATION], 
                                           #int(row[OUTPUT_INTERVAL]), 
                                           row[LOGFILE_PATH])
                
                self.x_len, x_units = row[X_DIM].split()
                self.x_len = float(self.x_len)
                self.y_len, y_units = row[Y_DIM].split()
                self.y_len = float(self.y_len)
                self.z_len, z_units = row[Z_DIM].split()
                self.z_len = float(self.z_len)
        with open(sim_Data.log_Location) as data_File:
            #skip the header line
            data_File.readline()
            csv_Data = csv.reader(data_File, delimiter = '\t')
            # populate all data in each frame>molecule>atom
            sim_Data.frames = self.load_Frames_Data (csv_Data)

        return sim_Data

    """
    Input:    csv.reader object which has the data file to read, and delimiter already defined
    Return:   a list of frame_Data objects with properties (including molecules>atoms) populated by data from input file
    """
    def load_Frames_Data(self, csv_Data):
        # Frames Data Constants
        STEP            = 0
        MOLECULE_ID     = 1
        MOECULE_NAME    = 2
        ATOM_ID         = 3
        ATOM_NAME       = 4
        ATOM_MASS       = 5
        A_POS_X         = 6
        A_POS_Y         = 7
        A_POS_Z         = 8
        A_VEL_X         = 9
        A_VEL_Y         = 10
        A_VEL_Z         = 11
        A_FOR_X         = 12
        A_FOR_Y         = 13
        A_FOR_Z         = 14
        
        prev_Frame = 0
        prev_Molecule = -1
        frames_Data = []
        for row in csv_Data:
            # load the frame data into frame_Data objects
            if prev_Frame != row[STEP]:
                # here we do int(float()) because the step is expressed in scientific notation for values over a million
                frames_Data.append(frame_Data(int(float(row[STEP])), self.x_len, self.y_len, self.z_len))
                prev_Frame = row[STEP]
            # load the molecule data into molecule_Data objects
            if prev_Molecule != row[MOLECULE_ID]:
                frames_Data[-1].molecules.append(molecule_Data(#int(row[MOLECULE_ID]), 
                                                               row[MOECULE_NAME]))
                prev_Molecule = row[MOLECULE_ID]
            # load the atom data into atom_Data objects
            # notice that Position, Velocity, and Force are 3 value tuples
            frames_Data[-1].molecules[-1].atoms.append(atom_Data(#int(row[ATOM_ID]), 
                                                                 row[ATOM_NAME], 
                                                                 float(row[ATOM_MASS]), 
                                                                 (float(row[A_POS_X]), float(row[A_POS_Y]), float(row[A_POS_Z])), 
                                                                 (float(row[A_VEL_X]), float(row[A_VEL_Y]), float(row[A_VEL_Z])), 
                                                                 (float(row[A_FOR_X]), float(row[A_FOR_Y]), float(row[A_FOR_Z]))
                                                                 )
                                                       )
        return frames_Data