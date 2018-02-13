##########################################################################
# This is a class definition for a simulation - this contains data about a
# specific simulation being analysed.
# Author: Chris McCormick
# Date: Summer 2016
##########################################################################

from classes.type_List import type_List

class simulation_Data(object):
    #def __init__(self, num_Steps, frame_Duration, report_Interval, log_Location):
    def __init__(self, log_Location):
        #self.num_Steps = num_Steps					# not used
        #self.frame_Duration = frame_Duration		# not used
        #self.report_Interval = report_Interval		# not used
        self.log_Location = log_Location
        self.frames = []
        self.atom_Type_List = type_List()
