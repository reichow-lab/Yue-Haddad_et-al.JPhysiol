#
#    Program    : Propogator.py
#    Author    : Bassam Haddad
#
#    This program replaces both 'initiator.py' & 'continuator.py' from the previous iteration of DetBa_2.sh.
#    Propogator.py has two methods built into it:
#
#        - initialize    : Initializes pop_mat based off the needs of the user: len-by-len for rates-pmf, or len-by-2 for histograms.
#
#            Inputs    : bin_min, bin_max, bin_size, Array_dimensions
#            Outputs    : Initialized pop_mat, log_file containing the input information, and the number of bins
#
#
#        - populate     : Populates the previously created pop_mat
#
#            Inputs    : Initialized pop_mat, numfiles, file-prefix
#            Outputs    : Populated pop_mat
#

import numpy as np
from tqdm import tqdm
import itertools as it

# Initializes pop_mat for subsequent calculations. MxN pop_mat, such that M is number of bins (num_bins) long, and either 2 or num_bins wide. The N=2 pop_mat is for simple histograms, and the N=num_bins pop_mat is for transition matrices.
# Additionally, initialize() creates a dictionary where the keys are the int(zcoord) and the values are the bin indeces.

def initialize(file_list, bin_size, outname, array_dim, d_col, bin_lim='auto'):
    bin_min = 0
    bin_max = 0

    # Find the min-max of the system ad-hoc
    if bin_lim == 'auto':
        for file in file_list:
            Data = open(file,'r')
            for line in Data:
                val = line.split()
                if val[0] == 'IonID:':
                    pass
                elif float(val[d_col]) <= bin_min:
                    bin_min = int(float(val[d_col]))
                elif float(val[d_col]) >= bin_max:
                    bin_max = int(float(val[d_col]))
                else:
                    pass
        bin_dim = min(abs(bin_min),abs(bin_max))
    else:
        bin_dim = int(bin_lim)

    # Pair it down 'til it's symmetric and evenly divisible by bin_size
    while (bin_dim*2 + 1) % bin_size != 0:
        bin_dim -= 1
    bin_min = bin_dim * -1
    bin_max = bin_dim
    print(f'bin_min = {bin_min}, bin_max = {bin_max}')
    pop_mat_length  =    int(abs(bin_min) + abs(bin_max) + 1)
    num_bins        =    pop_mat_length/bin_size
    num_bins        =    int(num_bins)

    # Create dictionary with int(zcoord) and bin_index
    ZtoBin  = {}
    bin     = 0
    counter = 1
    for z in range(bin_min,bin_max + 1):
        ZtoBin[z] = bin
        if counter % bin_size == 0:
            bin += 1
            counter += 1
        else:
            counter += 1

    with open(str(outname + '.log'), 'w') as log:
        log.write('bin_min: ' + str(bin_min) + '\n' + 'bin_max: ' + str(bin_max) + '\n' + 'bin_size: ' + str(bin_size) + '\n' + 'num_bins: ' + str(num_bins))
    if array_dim == 0:
        # output for histograms
        pop_mat = np.zeros((num_bins,2))
        return pop_mat,bin_min,bin_max,num_bins,ZtoBin,bin_dim
    else:
        # output for transition matrix
        pop_mat = np.zeros((num_bins,num_bins))
        return pop_mat,bin_min,bin_max,num_bins,ZtoBin,bin_dim

# This method populates a transition pop_mat from 1D data (e.g. ion trajectories along z-coordinate)

def populate(file_list, pop_mat, bin_min, bin_max, bin_s, num_bins, array_dim, d_col, lag_step, ZtoBin):

    def hist_pop(bin_now):
        pop_mat[bin_now, 1] = pop_mat[bin_now,1] + 1
        return 0

#####################################
#                                    #
#            Main Program            #
#                                    #
#####################################

    # Open file, and read all lines in
    for file in tqdm(file_list):
        with open(file, "r") as f:
            all_lines = f.read().splitlines()
        # generate list of indexes denoting new ions
        start_list = []
        for i in range(len(all_lines)):
            if all_lines[i].split()[0] == "IonID:":
                start_list.append(i)
            else:
                pass
        if len(start_list) == 0:
            start_list.append(0)
        start_list.append(-1)
        # Loop through each ion's index and process their data
        for i in range(len(start_list)-1):
            bin_j = 'new'
            for line in all_lines[start_list[i]:start_list[i+1]:lag_step]:
                if line.split()[0] == "IonID:":
                    pass
                elif abs(float(line.split()[d_col])) > bin_max:
                    pass
                elif bin_j == 'new':
                    bin_j = ZtoBin[int(float(line.split()[d_col]))]
                else:
                    bin_i = bin_j
                    bin_j = ZtoBin[int(float(line.split()[d_col]))]
                    if array_dim == 1:    # Rates calculation: choice == 'R' or 'M'
                        pop_mat[bin_j,bin_i] += 1
                    elif array_dim == 0:    # Histogram calculation: choice == 'H'
                        hist_pop(bin_j)
    return pop_mat
