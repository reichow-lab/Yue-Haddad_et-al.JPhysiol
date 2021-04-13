import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from glob import glob
from sys import argv
import pandas as pd
script, system, b = argv
FileList = glob(system+"*Tracking.log")
FileList.sort()
N = int(b)
# Bins: bin_range   count
Bins = [[]]
for i in range(1,N+1):
    Bins[0].append(i)
F = 1
for FILE in FileList:
    Bins.append([])
    with open(FILE, 'r') as f:
        all_lines = f.read().splitlines()
    # find the final Time, and cal
    FinalTime = float(all_lines[-3].split()[0])
    hold = 0
    for i in range(1,N+1):
        count = 0
        for line in all_lines:
            if line.split()[0] == "Time" or line.split()[0] == "Total" or line.split()[0] == "Last" or line.split()[0] == "Time(ns)":
                pass
            elif float(line.split()[0]) > (i-1)*(FinalTime/N) and float(line.split()[0]) <= (i)*(FinalTime/N):
                count += (int(line.split()[2]) - hold)
                hold = int(line.split()[2])
        Bins[F].append((count/(FinalTime/N))*160)
    F += 1
with open(system+"_Block.txt", 'w') as out:
    for i in range(N):
        out.write(f"{Bins[0][i]}\t")
        for j in range(1,len(FileList)+1):
            if (j) == len(FileList):
                out.write(f"{Bins[j][i]}\n")
            else:
                out.write(f"{Bins[j][i]}\t")
