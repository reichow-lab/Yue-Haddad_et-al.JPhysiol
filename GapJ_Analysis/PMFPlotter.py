import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from glob import glob
from sys import argv
import pandas as pd

#script, system, outname, palette = argv
#het = bool(int(input("Homotypic (0) or Heterotypic (1)? ")))
def PMFPlotter(system,outname,palette):
    het = True # We've decided to not symmetrize any of our data.
    RateList = glob(system+"*rate_final.txt_asym.txt")
    if het == False:
        PssList = glob(system+"*Pss_final.txt")
    elif het == True:
        print("debug")
        PssList = glob(system+"*Pss_final.txt_asym.txt")
    print(RateList,PssList)
    # Final: [Pore] [PMF_Pss] [PMF_rate] [PMF_rate-Pss] [label]
    # Final:   0        1         2            3           4
    Final = [[],[],[],[],[]]
    labels = ["Avg. PMF"]

    for FILE in RateList:
        with open(FILE, 'r') as file:
            for line in file:
                val = line.split()
                if val[0] == "Pore":
                    pass
                else:
                    Final[0].append(float(val[0]))
                    Final[2].append(float(val[1]))
                    Final[4].append(labels[0])
    for FILE in PssList:
        with open(FILE, 'r') as file:
            for line in file:
                val = line.split()
                if val[0] == "Pore":
                    pass
                else:
                    Final[1].append(float(val[1]))
    for i in range(len(Final[0])):
        Final[3].append(Final[2][i]-Final[1][i])

    # Calculate the effective voltage drop accross the system by taking the diff-
    # erence between the DiffPont(zmin) and DiffPont(zmax) and converting to mV.
    MinList,MaxList,MinAvg,MaxAvg = [[],[]],[[],[]],[[],[]],[[],[]]
    n_counter, n = 0, 0
    for i in range(len(Final[0])):
        if Final[0][i] <= -75:
            MinList[0].append(Final[3][i])
            MinList[1].append(n)
        elif Final[0][i] >= 75:
            MaxList[0].append(Final[3][i])
            MaxList[1].append(n)
        n_counter += 1
        if n_counter == (len(Final[0])/len(RateList)):
            n += 1
            n_counter = 0
    # In order to properly calculate the variance from these data, I needed to separate
    # the averages into their respective PMFs.
    holdMin,holdMax = [],[]
    for n in RateList:
        holdMin.append([])
        holdMax.append([])
    for x in range(len(MaxList[0])):
        holdMin[int(MinList[1][x])].append(MinList[0][x])
        holdMax[int(MaxList[1][x])].append(MaxList[0][x])
    for n in range(len(holdMin)):
        MinAvg[0].append(np.mean(holdMin[n]))
        MinAvg[1].append(np.var(holdMin[n]))
        MaxAvg[0].append(np.mean(holdMax[n]))
        MaxAvg[1].append(np.var(holdMax[n]))
    print(MinAvg[1],MaxAvg[1])
    # 0.04336 (V*mol)/Kcal
    voltageAvg = (np.absolute(np.mean(MaxAvg[0])) + np.absolute(np.mean(MinAvg[0])))*0.04336*1000
    voltageStD = np.sqrt((np.var(MaxAvg[0]) + np.var(MinAvg[0])))*(0.04336)*1000

    with open(outname + "_volt.log", 'w') as out:
        out.write("Pore-Axis\tPMF\tDriving-Potential\tDifference-Potential\tVoltage (Avg/StdDev)\n")
        for i in range(len(Final[0])):
            if i == 0:
                out.write(f"{Final[0][i]}\t{Final[1][i]}\t{Final[2][i]}\t{Final[3][i]}\t{voltageAvg}\n")
            elif i == 1:
                out.write(f"{Final[0][i]}\t{Final[1][i]}\t{Final[2][i]}\t{Final[3][i]}\t{voltageStD}\n")
            else:
                out.write(f"{Final[0][i]}\t{Final[1][i]}\t{Final[2][i]}\t{Final[3][i]}\n")
    print(len(Final[0]),len(Final[1]),len(Final[2]),len(Final[3]),len(Final[4]))
    PMFPot = pd.DataFrame({"Pore-Axis (A)":Final[0], "Energy (Kcal/mol)": Final[1]})
    DrivePot = pd.DataFrame({"Pore-Axis (A)":Final[0], "Energy (Kcal/mol)": Final[2]})
    DiffPot = pd.DataFrame({"Pore-Axis (A)":Final[0], "Energy (Kcal/mol)": Final[3]})
    plt.xlim(-85,85)
    plt.ylim(-0.5,5)
    sns.set_palette(palette)
    plt.title('PMF')
    plt.xlabel("Pore Axis (A)")
    plt.ylabel('Energy (Kcal/mol)')
    sns.lineplot(data=PMFPot, x="Pore-Axis (A)", y="Energy (Kcal/mol)", hue=Final[4], style=Final[4])
    plt.savefig(outname+"_PssPMF.png", dpi=400)
    plt.clf()
    plt.xlim(-85,85)
    plt.ylim(-2,5)
    plt.title('Driving Potential')
    plt.xlabel("Pore Axis (A)")
    plt.ylabel('Energy (Kcal/mol)')
    sns.lineplot(data=DrivePot, x="Pore-Axis (A)", y="Energy (Kcal/mol)", hue=Final[4], style=Final[4])
    plt.savefig(outname+"_RatesPMF.png", dpi=400)
    plt.clf()
    plt.xlim(-85,85)
    plt.ylim(-2.5,2.5)
    plt.title('Difference Potential')
    plt.xlabel("Pore Axis (A)")
    plt.ylabel('Energy (Kcal/mol)')
    sns.lineplot(data=DiffPot, x="Pore-Axis (A)", y="Energy (Kcal/mol)", hue=Final[4], style=Final[4])
    plt.savefig(outname+"_DiffPMF.png", dpi=400)
    plt.clf()
