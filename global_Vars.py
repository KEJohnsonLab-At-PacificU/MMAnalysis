##########################################################################
# This file stores values used in various parts of MMAnalysis
# Author: Chris McCormick
# Date: Summer 2016
##########################################################################

# Larger slice size yields smoother graphs
# For this we use 50 slices at 0.2 nm each
NUM_SLICES = 500
SLICE_SIZE = 0.02  # nanometers

# commonly used to reference X Y and Z indexes.
X = 0
Y = 1
Z = 2

#VDWS
# These are from parameter files provided by Kjersti
H_VDW = 1.062
O_VDW = 1.362

# These I looked up on the internet here: http://www.lenntech.com/periodic-chart-elements/vanderwaals.htm
C_VDW = 1.85
CL_VDW = 1.81

SURF_BOTTOM_TOLERANCE = 20
SURF_TOP_TOLERANCE = 100