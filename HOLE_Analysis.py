# 	Program: HOLE_Analysis.py
#	Author:	 Bassam Haddad
import numpy as np
from glob import glob
from sys import argv
import pickle as pkl

script, globstring = argv

# Create list of the appropriate HOLE files

hole_file_list = glob(globstring)
hole_file_list.sort()

# Load HOLE output files
up_lim  = 64
low_lim = -64
Pore_Radii = []
Pore_Axis  = np.arange(-64,65)

# Extract relevant data from HOLE output files
for hole_file in hole_file_list:
    with open(hole_file) as FileIN:
        temp_radii = []
        for line in FileIN:
            val = line.split()
            if (float(val[0]) <= up_lim and float(val[0]) >= low_lim):
                temp_radii.append(float(val[1]))
        if len(temp_radii) > 0:
            Pore_Radii.append(temp_radii)

# Save Extracted data for future processing
with open(str(globstring + '_data.pkl'), 'wb') as out:
    pkl.dump(Pore_Radii, out)
    pkl.dump(Pore_Axis, out)
