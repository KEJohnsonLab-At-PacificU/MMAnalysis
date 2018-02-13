##########################################################################
# This is a class definition for a molecule - this contains data about a
# specific molecule
# Author: Chris McCormick
# Date: Summer 2016
##########################################################################

class molecule_Data(object):
    #def __init__(self, molecule_ID, molecule_Name):
    def __init__(self, molecule_Name):
        #self.molecule_ID = molecule_ID		# we don't currently use molecule_ID, but it seems like a good idea to have
        self.molecule_Name = molecule_Name
        self.atoms = []