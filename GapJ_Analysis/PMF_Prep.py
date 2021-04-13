import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from math import sqrt
Error    = {}            # Dictionary keeping track of which "cut_num" provides the smallest error
#################################################################
def trim(PMF_in, cut_num, limit, final=False):    # cut_num is the number of values from PMF_for[pore-axis] to cut.
    PMF_for        = [[],[]]               # Forward PMF, [[pore-axis],[PMF_for],[MT]]
    PMF_for[0]     = list(PMF_in[0])
    PMF_for[1]     = list(PMF_in[1])
    PMF_rev        = [[],[]]               # Reverse PMF, [[pore-axis],[PMF_rev],[MT]]
    for i in range(cut_num):                    # This loop cuts the [pore-axis] accordingly, if cut_num == 0 then nothing happens
        if cut_num < 0:
            del PMF_for[0][-1]
        else:
            del PMF_for[0][0]
    while len(PMF_for[0]) < len(PMF_for[1]):    # This ensures that the two sub-lists are the same length prior to trimming to the limit
        if cut_num < 0:
            del PMF_for[1][0]
        else:
            del PMF_for[1][-1]
    while PMF_for[0][0] < -limit:            # These two loops trim the PMF to the pre-defined limits
        del PMF_for[0][0]
        del PMF_for[1][0]
    while (PMF_for[0][-1] > limit) or (PMF_for[0][-1] == 0):
        del PMF_for[0][-1]
        del PMF_for[1][-1]

    PMF_rev[0]    =    list(PMF_for[0][:])        # Creates the PMF reverse PMF...same pore axis
    PMF_rev[1]    =    list(PMF_for[1][::-1])    # but reversed PMF values...the list slice [::-1] is reversing the second column
    #return PMF_for, PMF_rev, cut_num        # Putting this on hold, I am thinking of calling the error function from here directly instead of returning these PMF values
    PMF_avg        =    error(PMF_for, PMF_rev, cut_num)    # I don't actually need the new PMF at this point..once I find the optimal cut number, I'll call it later.
    if final == True:
        return PMF_avg, PMF_for, PMF_rev
#################################################################
def error(PMF_for, PMF_rev, cut_num):
    PMF_avg        = [[],[],[]]    # Average PMF, [[pore-axis],[PMF_avg],[SEM]]. I am reassigning this everytime I call it to clear out the columns
    PMF_avg[0]     = list(PMF_for[0][:])
    hold           = np.zeros([1,2])
    for PMFf,PMFr in zip(PMF_for[1],PMF_rev[1]):
        PMF_avg[1].append(np.mean([PMFf,PMFr]))
        PMF_avg[2].append((np.std([PMFf,PMFr]))/sqrt(2))
    Error[cut_num] = sum(PMF_avg[2])
    return PMF_avg
#################################################################
def final(PMF_avg):
    if len(PMF_avg) == 3:
        PMF_fin = [[],[],[],[]]         # Final PMF  , [[pore-axis],[PMF_avg_adj],[Avg+SEM],[Avg-SEM]]
        PMF_fin[0]     = PMF_avg[0][:]
        for i in range(0,len(PMF_avg[0]),1):
            PMF_fin[1].append(PMF_avg[1][i] - PMF_avg[1][0])
            PMF_fin[2].append(PMF_fin[1][i] + PMF_avg[2][i])
            PMF_fin[3].append(PMF_fin[1][i] - PMF_avg[2][i])
        return PMF_fin
    elif len(PMF_avg) == 2:
        PMF_fin = [[],[]]
        PMF_fin[0] = PMF_avg[0][:]
        for i in range(len(PMF_avg[0])):
            PMF_fin[1].append(PMF_avg[1][i] - (((PMF_avg[1][0])+(PMF_avg[1][-1]))/2))
        return PMF_fin
#################################################################
def interp(PMF_in):
    PMF_IN  =   [[],[]]
    PMF_fix =   [[],[]]
    with open(PMF_in, 'r') as data:
        for line in data:
            val     = line.split()
            PMF_IN[0].append(float(val[0]))
            PMF_IN[1].append(float(val[1]))
    del PMF_IN[0][-1]  # not sure why I had these here in the first place
    del PMF_IN[1][-1]  # but the code won't run without it.
    xin     =   np.array(PMF_IN[0])
    yin     =   np.array(PMF_IN[1])
    f       =   interpolate.CubicSpline(xin,yin)  # Cubic-spline interpolation
    xout    =   np.arange(PMF_IN[0][0],PMF_IN[0][-1] + 1,1) # Creating the correct dimension for the data output.
    yout    =   f(xout)
    PMF_fix[0]  =   xout.tolist()
    PMF_fix[1]  =   yout.tolist()
    return  PMF_fix
#################################################################
def HetCenter(PMF_in):
    PeakList,hold = [],1000
    for i in range(3,len(PMF_in[0])-3):
        if (PMF_in[1][i-3] < PMF_in[1][i]) and (PMF_in[1][i+3] < PMF_in[1][i]):
            PeakList.append(i)
        else:
            pass
    # Find the peak nearest x = 0
    print(PeakList)
    for i in PeakList:
        d = abs(PMF_in[0][i])
        if d <= hold:
            hold = d
            holdi = i
        else:
            pass
    return -round(PMF_in[0][holdi]+0.5)

#################################################################
#                                                               #
#                         Main Program                          #
#                                                               #
#################################################################
def Prep(PMF_in, outname, bin_dim, het):
    
    CUT_NUMS = list(range(-10,11))    # Eventually I want to have it more dynamically search for cut nums, but just performing the calculation for all
    limit = min(85,bin_dim)
    PMF_fix  =  interp(PMF_in)
    PMF_trim =  list(PMF_fix)
    if het == False:
        for x in CUT_NUMS:
            trim(PMF_trim, x, limit)
        BESTCUT = min(Error, key=Error.get)
    elif het == True:
        BESTCUT = HetCenter(PMF_trim)
    Average_PMF,PMF_for,PMF_rev    = trim(PMF_fix, BESTCUT, limit, True)
    Final_PMF    = final(Average_PMF)
    Final_For    = final(PMF_for)
    Final_Rev    = final(PMF_rev)
    with open(outname, "w") as ofile:
        ofile.write("Pore Axis\tAvg PMF\tAvg + SEM\tAvg - SEM\n")
        for i in range(len(Final_PMF[0])):
            ofile.write(str(Final_PMF[0][i])+"\t"+str(Final_PMF[1][i])+"\t"+str(Final_PMF[2][i])+"\t"+str(Final_PMF[3][i])+"\n")
    with open(str(outname + "_asym.txt"), "w") as ofile:
        ofile.write("Pore Axis\tForward PMF\tReverse PMF\n")
        for i in range(len(Final_For[0])):
            ofile.write(str(Final_For[0][i])+"\t"+str(Final_For[1][i])+"\t"+str(Final_Rev[1][i])+"\n")
