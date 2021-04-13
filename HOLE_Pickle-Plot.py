#	Program: HOLE_Pickle-Plot.py
#	Author:	 Bassam Haddad

import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt
import sys
import argparse

# Parse inputs
parser = argparse.ArgumentParser()
parser.add_argument("-dat", dest = "datfile", action = "store")
parser.add_argument("-out", dest = "outname", action = "store", default = "OUTFILE")
parser.add_argument("-c", "--choice", dest = "choice", action = "store")
parser.add_argument("-exp", dest = "exp", action = "store")
parser.add_argument("-color", dest = "color", action = "store")
parser.add_argument("-sym")
args = parser.parse_args()

if args.datfile == None:
	sys.exit("Missing input file.")
if args.outname == None:
	print("Output name is unspecified. Using default: OUTFILE.png/txt")
if args.choice != 'apbs' and args.choice != 'hole':
	print("make a choice! hole or apbs")
if args.color == None:
	print("Pick a color")

# Set color scheme

Color_Options = {'red': ['firebrick','lightcoral','lightcoral'], 'blue': ['royalblue','skyblue','skyblue'], 'green': ['olivedrab','palegreen','palegreen'], 'purple': ['rebeccapurple','mediumpurple','mediumpurple']}

Color_Scheme  = Color_Options[args.color]

with open(args.datfile, 'rb') as datin:

	CenterPots = pkl.load(datin)
	Pore_Axes  = pkl.load(datin)

# Calculate mean and error
if args.sym:
	CenterPotsSym = []
	for i in range(0,len(CenterPots),1):
		hold = []
		hold.append(CenterPots[i])
		hold.append(hold[0][::-1])
		CenterPotsSym.append(np.mean(hold, axis=0))
	Final = []
	Final.append(np.mean(CenterPotsSym, axis=0))
	Final.append(np.mean(CenterPotsSym, axis=0) + np.std(CenterPotsSym, axis=0))
	Final.append(np.mean(CenterPotsSym, axis=0) - np.std(CenterPotsSym, axis=0))

else:
	print('testing not sym')
	Final = []
	Final.append(np.mean(CenterPots, axis=0))
	Final.append(np.mean(CenterPots, axis=0) + np.std(CenterPots, axis=0))
	Final.append(np.mean(CenterPots, axis=0) - np.std(CenterPots, axis=0))

# Import Experimental data

if args.exp != None:
	with open(args.exp, 'rb') as datain:
		experiment = pkl.load(datain)
# Write out .txt file for Excel

	with open((args.outname + '.txt'), 'w') as excelout:
		excelout.write('Pore-Axis\tExp\tMD\t+SEM\t-SEM\n')
		for i in range(0,len(Pore_Axes),1):
			excelout.write(str(Pore_Axes[i])+'\t'+str(experiment[0][i])+'\t'+str(Final[0][i])+'\t'+str(Final[1][i])+'\t'+str(Final[2][i])+'\n')
		excelout.close()

elif args.exp == None:
	with open((args.outname + '.txt'), 'w') as excelout:
		excelout.write('Pore-Axis\tMD\t+SEM\t-SEM\n')
		for i in range(0,len(Pore_Axes),1):
			excelout.write(str(Pore_Axes[i])+'\t'+str(Final[0][i])+'\t'+str(Final[1][i])+'\t'+str(Final[2][i])+'\n')
		excelout.close()
# Plot Data

if (args.choice == "apbs") and (args.exp != None):
	for i in range(0,len(Final),1):
		plt.plot(Pore_Axes,Final[i],color=Color_Scheme[i])
	plt.plot(Pore_Axes,experiment[0],color='black')
	plt.title(input("Title? "))
	plt.xlabel('Pore-Axis (A)')
	plt.ylabel('Potential (kT/e)')
	plt.ylim(-20,10)
	plt.savefig(args.outname,dpi=600)

elif (args.choice == "apbs") and (args.exp == None):
	for i in range(0,len(Final),1):
		plt.plot(Pore_Axes,Final[i],color=Color_Scheme[i])
	plt.title(input("Title? "))
	plt.xlabel('Pore-Axis (A)')
	plt.ylabel('Potential (kT/e)')
	plt.ylim(-20,10)
	plt.savefig(args.outname,dpi=600)

elif (args.choice == "hole") and (args.exp != None):
	for i in range(0,len(Final),1):
		plt.plot(Pore_Axes,Final[i],color=Color_Scheme[i])
	plt.plot(Pore_Axes,experiment[0],color='black')
	plt.title(input("Title? "))
	plt.xlabel('Pore-Axis (A)')
	plt.ylabel('Pore Radii (A)')
	plt.ylim(0,15)
	plt.savefig(args.outname,dpi=600)

elif (args.choice == "hole") and (args.exp == None):
	for i in range(0,len(Final),1):
		plt.plot(Pore_Axes,Final[i],color=Color_Scheme[i])
	plt.title(input("Title? "))
	plt.xlabel('Pore-Axis (A)')
	plt.ylabel('Pore Radii (A)')
	plt.ylim(0,15)
	plt.savefig(args.outname,dpi=600)
