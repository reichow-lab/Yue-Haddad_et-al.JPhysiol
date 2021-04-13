#
#    Program: GapJ_Analysis.py
#
#    The purpose of this program is to replace the current version GapJ_Analysis_2.sh, a bash script that strung together a series of python scripts in an
#    ad-hoc fashion that was clunky, and potentially prone to error. This version will be entirely contained within this python program that uses the
#    aforementioned analysis scripts as modules.
#
#    initializer.py and continuator.py will be replaced by a singel program that does both jobs, since I can write a python conditional that can distinguish
#    between the first and not-first ion trajectories. This means that I will have to go into thos scripts (more likely create a new one) that creates functions that
#    can be called by this (GapJ_Analysis.py) program.

import numpy as np
import matplotlib.pylab as plt
import Ion_Tracker
from Propagator import initialize,populate
from Calculator import sympop,rate2gibbs,hist_write,mfpt,check_SS
from Current_Calculator import Current,Text2PMF,VoltPMF
from Diffusion_Calc import normalize,Diff_Calc
from Edge_Erase import edge_erase,tri_diag
from PMF_Prep import Prep
from sys import argv
from glob import glob

script, globstring = argv

#################################
#                               #
#    Initial Variables          #
#                               #
#################################

file_list       =    glob(globstring)
bin_size        =    float(input("What is the desired bin size? "))
outname         =    str(input("What would you like to name this project? "))
choice          =    str(input("What would you like to do? MSM PMF (M), Histogram (H), Ion Tracker (T), I-V approximator (I) "))
choice          =    choice.upper()
out_pop_mat     =    outname + "_pop.mat"
out_rate_mat    =    outname + "_rate.mat"
out_IV          =    outname + "_I-V.data"
out_final       =    outname + "_Pss_final.txt"
lag_base        =    int(input("What is the step-size per frame? (answer in picoseconds) "))  # the base lag_time, or frame is 2ps...this can be softcoded later.

#################################
#                               #
#    Main Program               #
#                               #
#################################

END    =    False

while END == False:
    if choice == 'M':
        #het         =   bool(int(input("Homotypic (0), or Heterotypic (1)? ")))
        het         =   True
        d_col       =   int(input("Which column from your data_file will you use? "))
        lag_step    =   int(int(input(f"Choose a lag time. (multiple of {lag_base}ps) "))/lag_base)
        lag_time    =   lag_step * lag_base
        bin_lim     = input('What is the Bin limit? ') # "auto" is acceptable
        array_dim   =    1
        init_matrix,bin_min,bin_max,num_bins,ZtoBin,bin_dim = initialize(file_list, bin_size, outname, array_dim, d_col, bin_lim)
        first_center = bin_min + (bin_size/2)
        pop_matrix  =    populate(file_list, init_matrix, bin_min, bin_max, bin_size, num_bins, array_dim, d_col, lag_step, ZtoBin)
        pop_matrix.dump(out_pop_mat)
        rate_matrix = normalize(pop_matrix)
        rate_matrix.dump(out_rate_mat)
        gibbs       =   rate2gibbs(num_bins, first_center, rate_matrix, bin_size, str(outname + '_rate'))
        Prep(gibbs, str(outname + '_rate_final.txt'), bin_dim, het)
        source      = int(input("Which bin is the source? "))
        sink        = int(input("which bin is the sink? "))
        gibbs,K_AB,MFPT,MSM,Pss =   mfpt(pop_matrix,num_bins,outname,source,sink,bin_min,bin_max,bin_size,ZtoBin,lag_time)
        Prep(gibbs, out_final, bin_dim, het)
        check_SS(MSM,Pss,num_bins,lag_time,outname,bin_size)
    elif choice == 'H':
        d_col       = int(input("Which column from your data_file will you use? "))
        array_dim   = 0
        lag_step    = 1
        bin_lim     = input('What is the Bin limit? ')
        init_matrix,bin_min,bin_max,num_bins,ZtoBin,bin_dim = initialize(file_list, bin_size, outname, array_dim, d_col, bin_lim)
        pop_matrix  =    populate(file_list, init_matrix, bin_min, bin_max, bin_size, num_bins, array_dim, d_col, lag_step, ZtoBin)
        write_mat    =    hist_write(bin_min, pop_matrix, outname, bin_size, num_bins)
    elif choice == 'T':
        d_col   =    int(input("Which column from your data_file will you use? "))
        ION        =    Ion_Tracker.ION()
        ION.tracker(file_list, d_col, outname, lag_base)
        Ion_Tracker.process(outname, lag_base)
    elif choice == 'A':
        ION        =    Ion_Tracker_DEV.ION()
        ION.tracker(num_files, prefix, outname)
    elif choice == 'I':
        PMF_txt = input("which PMF file would you like to use? ")
        PMF_0,num_bins  = Text2PMF(PMF_txt)
        pop_mat_EE      = edge_erase(np.load(out_pop_mat,allow_pickle=True),bin_size)
        trans_mat       = normalize(pop_mat_EE)
        Diff_pore       = Diff_Calc(trans_mat,bin_size)                         # Calculates in units (A^2/second)
        #print    "Diff_pore = %s" % Diff_pore
        Voltages        =    [-150,-100,-75,-50,-25,0,25,50,75,100,150]
        out            =    open(out_IV, 'w')
        out.write('dV (mV)' + ',' + 'Current (pA)' + ',' + 'Forward MFPT' + ',' + 'Reverse MFPT' + '\n')
        for V in Voltages:
            dV        =    float(V)
            pmf_v        =    VoltPMF(PMF_0,dV,num_bins)
            I,Tau_f,Tau_r    =    Current(pmf_v,num_bins,bin_size,Diff_pore)
            out.write(str(dV) + ',' + str(I) + ',' + str(Tau_f) + ',' + str(Tau_r) + '\n')
        out.close()
    choice        =    str(input("What would you like to do? MSM PMF (M), Histogram (H), Ion Tracker (T), I-V approximator (I), Exit (E) "))
    choice        =    choice.upper()
    if choice    ==    'E':
        END    =    True
    elif choice    in    {"M","H","T","I"}:
        END    =    False
        bin_size        =    float(input("What is the desired bin size? "))
    else:
        END    =    True
