import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from glob import glob
from sys import argv
import pandas as pd
from scipy.interpolate import interp1d
from scipy.signal import find_peaks
#script, system, start, outname, palette, WS, obs = argv
def Interp(xin,yin,LT):
    f = interp1d(xin,yin,kind="previous")
    xnew = np.arange(0,LT,1)
    ynew = f(xnew)
    return xnew,ynew
def count_dist_peaks(series, bins=50):
    count, division = np.histogram(series, bins=bins)
    peaks, props = find_peaks(count, prominence=50)
    return peaks, count
def TrackerPlot(system,start,outname,palette,WS,obs,LT,d_col,ObString):
    LT = LT + WS
    if bool(obs) == True:
        ObsFileList = glob(ObString)
        ObsFileList.sort()
    WinS = int(WS)
    FileList = glob(system+"*Tracking.log")
    FileList.sort()
    # Final: [time] [cum. permeations] [label] [hue] [Average Current]
    Final = [[],[],[],[],[],[]]
    labels, hues = [],[]
    elemC,convF = 1.60217662e-19,(1e-9/1e-12)
    for i in range(len(FileList)):
        hues.append(i)
        labels.append(input("Label? "))
    # Import Data and interpolate. Also import the dts
    # fpt or First Passage Time (properly the first transition time)
    fptList = [[],[]]
    for FILE,label,hue in zip(FileList,labels,hues):
        with open(FILE, 'r') as file:
            Semi = [[],[],[],[],[]]
            for line in file:
                val = line.split()
                if val[0] == "Time":
                    pass
                elif val[0] == "Total" or val[0] == "Last":
                    pass
                elif float(val[0]) < start:
                    pass
                else:
                    Semi[0].append(float(val[0]))
                    Semi[1].append(int(val[2]))
                    Semi[2].append(str(label))
                    Semi[3].append(int(hue))
                    if float(val[0]) <= 0:
                        pA = 0
                    elif float(val[0]) > 0:
                        pA = (int(val[2])*elemC)/(float(val[0])*1e-9)/(1e-12)
                    Semi[4].append(pA)
                    fptList[0].append(float(val[1]))
                    fptList[1].append(0)
            xnew,ynew = Interp(Semi[0],Semi[1],LT)
            # Process the interpolated data
            for i in range(len(xnew)):
                Final[0].append(float(xnew[i]))
                Final[1].append(int(ynew[i]))
                Final[2].append(str(label))
                Final[3].append(int(hue))
                if i == 0:
                    Final[4].append(0)
                # Calculate instantaneous Slopes
                else:
                    Curr = float((ynew[i]-ynew[i-1])/(xnew[i]-xnew[i-1]))
                    Final[4].append(Curr)
                # Calculate Cumulative Current
                if xnew[i] <= 0:
                    pA = 0
                elif xnew[i] > 0:
                    pA = (int(ynew[i])*elemC)/(float(xnew[i])*1e-9)/(1e-12)
                Final[5].append(pA)
    # Separate to calculate running-averages
    HoldSep = []
    for i in range(len(FileList)):
        HoldSep.append([])
    n = -1
    for i in range(len(Final[0])):
        if float(Final[0][i]) <= 0:
            n += 1
        HoldSep[n].append(Final[4][i])
    # RunAvg: Time  Current  Hue  Label
    RunAvg = [[],[],[],[]]
    h = 0
    for n in HoldSep:
        for i in range(len(n)-WinS):
            RunAvg[0].append(Final[0][i])
            hold = []
            for j in range(WinS):
                hold.append(float(n[i+j]))
            RunAvg[1].append(np.mean(hold)*160)
            RunAvg[2].append(h)
            RunAvg[3].append(labels[h])
        h += 1
    with open(outname+'_wa.txt', 'w') as out:
        for i in range(len(RunAvg[0])):
            out.write(f"{RunAvg[0][i]}\t{RunAvg[1][i]}\t{RunAvg[2][i]}\n")
    if bool(obs) == True:
        # Perform the same data formatting for the new observable such that it has the same window averaging as the current
        Obs, Semi, WinObs = [[]], [[]], [[]]
        # Open file, and read all lines in
        for file in ObsFileList:
            with open(file, "r") as f:
                all_lines = f.read().splitlines()
            i = 0
            for line in all_lines:
                if line.split()[0] == "Chain:":
                    Obs.append([])
                    Semi.append([])
                    WinObs.append([])
                    i += 1
                elif i == 1:
                    Semi[0].append(float(line.split()[0])/10)
                else:
                    pass
            i = 0
            for line in all_lines:
                if line.split()[0] == "Chain:":
                    i += 1
                else:
                    Semi[i].append(float(line.split()[d_col]))
            # Process the interpolated data
            for i in range(1,len(Semi)):
                xnew,ynew = Interp(Semi[0],Semi[i],LT)
                print(xnew)
                for ii in range(len(xnew)):
                    if i == 1:
                        Obs[0].append(float(xnew[ii]))
                        Obs[i].append(int(ynew[ii]))
                    else:
                        Obs[i].append(int(ynew[ii]))
        # Data is already in wide-format so no need for separating as we do above.
        for c in range(1,len(Obs)):
            for i in range(len(Obs[0])-WinS):
                if c == 1:
                    WinObs[0].append(Obs[0][i])
                    hold = []
                    for j in range(WinS):
                        hold.append(float(Obs[c][i+j]))
                    WinObs[c].append(np.mean(hold))
                else:
                    hold = []
                    for j in range(WinS):
                        hold.append(float(Obs[c][i+j]))
                    WinObs[c].append(np.mean(hold))

    ################################################################################
    plot_data1 = pd.DataFrame({"Time (ns)": Final[0], "Ion Permeations": Final[1]})
    plot_data2 = pd.DataFrame({"Time (ns)": RunAvg[0], "current (pA)": RunAvg[1]})
    plot_data3 = pd.DataFrame({"Time (ns)": Final[0], "<current> (pA)": Final[5]})
    plot_data4 = pd.DataFrame({"Transition Times (ns)": fptList[0]})
    if bool(obs) == True:
        #print(len(RunAvg[1]),len(WinObs[1]))
        plot_data5 = pd.DataFrame({"Current (pA)": RunAvg[1]})
        for i in range(1,len(WinObs)):
            plot_data5[str(i)] = pd.Series(WinObs[i])
        print(plot_data5)
        plt.title("Current vs. Observable")
        plt.xlabel("Current")
        plt.ylabel("Observable")
        for i in range(1,len(WinObs)):
            sns.scatterplot(data=plot_data5, x="Current (pA)", y=WinObs[i], linewidth=0)
        plt.savefig(outname+"_ObsVsCurr.png")
        plt.clf()
    fig, ax = plt.subplots()
    #plt.xlim(0,260)
    #plt.ylim(0,200)
    plt.title("Current")
    plt.xlabel("Time (ns)")
    plt.ylabel('Windowed Avg. Current (pA)')
    sns.lineplot(data=plot_data2, x="Time (ns)", y="current (pA)", hue=RunAvg[3], palette=sns.color_palette(palette, n_colors=len(FileList)))
    plt.savefig(outname+"_current.png", dpi=400)
    plt.clf()
    sns.displot(data=RunAvg[1], kind="kde")
    plt.xlim(0,160)
    plt.ylim(0,0.025)
    plt.savefig(outname+"_current-hist.png", dpi=400)
    plt.clf()
    #plt.xlim(0,260)
    #plt.ylim(0,500)
    #sns.set_palette(palette, 4)
    fig, ax = plt.subplots()
    plt.title("Flux")
    plt.xlabel("Time (ns)")
    plt.ylabel('Ion Permeations')
    sns.lineplot(data=plot_data1, x="Time (ns)", y="Ion Permeations", hue=Final[2], palette=sns.color_palette(palette, n_colors=len(FileList)))
    plt.savefig(outname+"_flux.png", dpi=400)
    plt.clf()
    fig, ax = plt.subplots()
    #sns.set_palette(palette, 4)
    plt.title("Cumulative Current")
    plt.xlabel("Time (ns)")
    plt.ylabel("Running Avg. Current")
    sns.lineplot(data=plot_data3, x="Time (ns)", y="<current> (pA)", hue=Final[2], palette=sns.color_palette(palette, n_colors=len(FileList)))
    plt.savefig(outname+"_CumCurrent.png", dpi=400)
    plt.clf()
    fig, ax = plt.subplots()
    ax2 = ax.twinx()
    np.set_printoptions(precision=3)
    plt.title("Ionic Transition Times")
    plt.xlabel("Transition Times (ns)")
    ax.set_ylabel("Probability Mass Function")
    ax2.set_ylabel("Cumulative Distribution Function")
    ax = sns.histplot(data=plot_data4,stat='probability',palette=sns.color_palette(palette, n_colors=1))
    ax2 = sns.ecdfplot(data=plot_data4,stat='proportion',palette=sns.color_palette("Greys_r", n_colors=1))
    plt.savefig(outname+"_FPT.png", dpi=400)
    plt.clf()

    return RunAvg
