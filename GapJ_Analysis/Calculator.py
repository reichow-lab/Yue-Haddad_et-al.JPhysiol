#
#    Program    : Calculator.py
#    Author    : Bassam Haddad
#
#    This module replaces both 'pop2matrix.py', and 'rate2gibbs.py' from the previous iteration of DetBa.sh.
#    Calculator.py has three methods build into it:
#
#        - pop2rate    : Creates a rate-matrix, calculated from the population matrix; which was calculated in 'Propogator.py'
#
#            Inputs    : Initial bin (init), Final bin (fina), population matrix (pop_matrix), bin size (bin_size)
#            Outputs    : Rate Matrix (markov)
#
#
#        - rate2gibbs    : Uses values from the matrix matrix, and calculates the gibbs free energy for the transitions between bins.
#
#            Inputs     : Initial bin (init), Final bin (fina), rates matrix (matrix), bin size (bin_size), output file (outname)
#            Outputs    : 'Free-Energy vs. Position'

import numpy as np
from tqdm import tqdm
from Diffusion_Calc import normalize
# This isn't actually called in the program, but is here as tool for testing random transition matrices.
def sympop(bin_min, bin_size, pop_matrix, ZtoBin):
    bin_list = np.arange(bin_min,abs(bin_min),bin_size)
    print(bin_list)
    for i in bin_list:
        i = ZtoBin[int(i)]
        for j in bin_list:
            j = ZtoBin[int(j)]
            pop_matrix[i,j] = np.mean([pop_matrix[i,j],pop_matrix[-i,-j]])
            pop_matrix[-i,-j] = pop_matrix[i,j]
    return pop_matrix

def rate2gibbs(num_bins, center, tran_matrix, bin_size, outname):
    sum     = 0
    gibbs     = np.zeros((num_bins,2))
    i         = 0
    def det_bal(forw, rev):
        rel        = forw/rev
        gibb     = -0.00198588*310*np.log(rel)
        return gibb
    for j in range(1,num_bins):    # This loop is the first to populate the matrix[i,0] column. I need to incrememnt center by bin_size instead of 1.
        forw    = tran_matrix[j,i]
        rev        = tran_matrix[i,j]
        gib        = det_bal(forw, rev)
        gibbs[i,0] = center
        gibbs[i,1] = gib
        center   = center + bin_size
        i += 1
    for b in range(num_bins):
        sum = sum + gibbs[b,1]
        gibbs[b,1] = sum
    out_file = str(outname + '_penult.txt')
    with open(out_file, 'w') as out:
        for c in range(num_bins):
            bin = gibbs[c,0]
            gib = gibbs[c,1]
            out.write(str(bin) + ' ' + str(gib) + ' ' + '\n')
        return out_file                # This allows PMF_Prep.py to take the name of the created text-file and use it as an input.

def mfpt(count_mat, num_bins, outname, source, sink, bin_min, bin_max, bin_size, ZtoBin, lag_time):
    # source and sink are bins (A,B) that that will direct where counts in the transition matrix will transfer from and to.
    # Counts from T_sink->j, to T_sink->source
    # Modify the transition (count) matrix by moving all counts from sink to source
    #if source != sink:
        # Convert from z (Å) to bin (i)
        #source  = ZtoBin[source]
        #sink    = ZtoBin[sink]
        # John Russo's Method
        #for j in range(num_bins):
        #    count_mat[source,j] += count_mat[sink,j]
        #    count_mat[sink,j]   = 0
        #for j in range(num_bins):
            #count_mat[source,sink] += count_mat[]
    # recalculate the transition (probability/rate) matrix
    tran_mat = normalize(count_mat)
    tran_mat.dump(str(outname + '_MFPT.mat'))
    # Calculate the stationary state probabilities (Pss)
    P   = np.zeros([num_bins,1])
    P[-1]   = 1
    ones    = np.ones([1,num_bins])
    for i in range(num_bins):
        tran_mat[i][i] -= 1
    temp_mat    = np.delete(tran_mat,-1,axis=0)
    aug_mat     = np.append(temp_mat,ones,axis=0)
    Pss         = np.linalg.solve(aug_mat,P)
    with open(str(outname + '_penult.txt'), 'w') as out:
        i = bin_min + (bin_size/2)
        for p in Pss:
            out.write(str(i)+'\t'+str(-0.00198588*310*np.log(p[0]))+'\n')
            i += bin_size
    # This step only occurs if system is in "non-equilibrium"
    if source != sink:
        # Calculate the rate from source to sink (K_AB)
        K_AB = 0
        # Create the i-list (bins in state A) & j-list (bins in state B)
        # (recall, T_ji corrosponds to the conditional probability that and ion
        # transitions from bin i to j)
        i_list,j_list = [],[]
        if sink < source:
            pref = 1
        else:
            pref = -1
        for i in range(0,num_bins,1):
            if i < (pref*sink):
                i_list.append(i)
            else:
                j_list.append(i)
        # Perform double sum over i's & j's
        for i,p_i in zip(i_list,Pss):
            for j in j_list:
                K_AB    = K_AB + p_i[0]*tran_mat[j,i]
        K_AB    = (1/lag_time)*K_AB
        # Calculate the mean first passage time
        MFPT    = 1/K_AB
        print(f"\n The MFPT is {MFPT} ps.")
        return str(outname + '_penult.txt'),K_AB,MFPT,tran_mat,Pss
    return str(outname + '_penult.txt'),0,0,tran_mat,Pss

def check_SS(MSM,Pss,num_bins,lag_time,outname,bin_size):
    """
    When a system is in a steady state (SS) then the flux of probability from
    state i -> i+1 is equal to the flux into bin-i must be equal to the flux out
    of bin-i. This method takes in a transition matrix (i.e. MSM) and calculates
    the distribution of flux to the nearest neighbor. The final output will be a
    coefficient between 0 and 1 where 1 is a perfect steady state.
    """

    with open(outname + '_cSS.txt', 'w') as outss:
        outss.write("Pss\tI-flux/O-flux\tJ_i-j\n")
        State_list = list(range(num_bins))
        hold = 0
        for state in State_list:
            # Calculate the flux in and out of each state (i.e. bin)
            iflux  = 0
            oflux = 0
            # Create the i-list (bins in state A) & j-list (bins in state B) (recall, T_ji corrosponds to the conditional probability that and ion transitions from bin i to j)
            j_list = list(range(num_bins))
            j_list.pop(state)
            #flux into 'state' and out of
            for j in j_list:
                iflux += (Pss[j,0])*MSM[state,j]
                oflux += (Pss[state,0])*MSM[j,state]

            #flux through the i->j surface
            J_ij = 0
            edge_bins1, edge_bins2 = list(range(int(20/bin_size))), list(range(int(180/bin_size),num_bins))
            for i in range(state+1):
                for j in range(state+1,num_bins):
                    if (i in edge_bins1 and j in edge_bins2) or (i in edge_bins2 and j in edge_bins1):
                        pass
                    else:
                        J_ij += (Pss[j,0] * MSM[i,j]) - (Pss[i,0] * MSM[j,i])

            outss.write(f"{Pss[state][0]}\t{iflux/oflux}\t{J_ij}\n")


def hist_write(init, pop_matrix, outname, bin_size, num_bins):
    bin_init = int(init)
    counter = 0
    row_tot = 0
    num_bins = int(num_bins)
    outfile        = outname + '_hist.dat'
    out = open(outfile, 'w')
    for i in tqdm(range(0,num_bins,1)):
        val1    = bin_init
        val2    = pop_matrix[i,1]
        val_1    = str(val1)
        val_2    = str(val2)
        out.write(val_1 + ' ' + val_2 + '\n')
        bin_init = bin_init + bin_size
    out.close()
