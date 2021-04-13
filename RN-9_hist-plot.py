# plot_data.py
# by Umair Khan
# Reichow Lab

# Plot histograms of interactions from data files.

# Usage:
#  -t: title of plot
#  -p: plot command specified "-p [path to file] [column] [bin start] [bin end] [bin size] [label] [color]"
#      colors from https://matplotlib.org/2.0.2/examples/color/named_colors.html
#  -s: optional -- filepath to save figure to (instead of showing)
#  -nl: optional -- remove the legend

# Imports
import sys
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

# Styling
sns.set()
sns.set_context("talk")
matplotlib.rcParams["font.family"] = "HK Grotesk"
matplotlib.rcParams["font.weight"] = "semibold"
matplotlib.rcParams["axes.titleweight"] = "bold"
matplotlib.rcParams["axes.labelweight"] = "medium"
matplotlib.rcParams["axes.titlesize"] = "xx-large"
matplotlib.rcParams["axes.labelsize"] = "large"
matplotlib.rcParams['xtick.labelsize'] = "large"
matplotlib.rcParams['ytick.labelsize'] = "large"
matplotlib.rcParams['text.color'] = "#333333"
matplotlib.rcParams['axes.labelcolor'] = "#333333"
matplotlib.rcParams['xtick.color'] = "#333333"
matplotlib.rcParams['ytick.color'] = "#333333"
matplotlib.rcParams["axes.facecolor"] = "white"
matplotlib.rcParams["axes.edgecolor"] = "#333333"
matplotlib.rcParams["lines.linewidth"] = 5.0
matplotlib.rcParams["axes.spines.top"] = False
matplotlib.rcParams["axes.spines.right"] = False

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-t", "--title", dest = "title", action = "store")
parser.add_argument("-p", "--plot", dest = "plot", nargs = "*", action = "store")
parser.add_argument("-s", "--save", dest = "save", action = "store", default = "show")
parser.add_argument("-nl", "--no-legend", dest = "no_legend", action = "store_true", default = False)
args = parser.parse_args()

if args.title == None:
	sys.exit("No title specified. Exiting...")

if args.plot == None:
	sys.exit("Please specify at least one histogram to plot. Exiting...")

params = {"filepath": args.plot[0],
		  "col": int(args.plot[1]),
		  "start": float(args.plot[2]),
		  "end": float(args.plot[3]),
		  "step": float(args.plot[4]),
		  "label": args.plot[5],
		  "color": args.plot[6]}

# Grab data from file
data = [float(line.split()[params["col"]]) for line in open(params["filepath"]).readlines()[1:]]

# Make list of bin edges
bin_edges = list(np.arange(params["start"], params["end"], params["step"])) + [params["end"]]

# Create figure and plot histogram
fig = plt.figure(figsize = (16, 10))
plt.hist(data, bins = bin_edges, density = True, color = params["color"], label = params["label"])
if params["col"] in [2, 4]:
	plt.xlabel(r"distance (${\rm \AA}$)", labelpad = 15)
else:
	plt.xlabel("angle (°)", labelpad = 15)
ttl = plt.title(args.title, pad = 20)
ttl.set_position([0.5, 1.03])
plt.tight_layout(pad = 2)
lgd = plt.legend(loc = "upper right")
if args.no_legend:
	lgd.set_visible(False)

# Save or show
if args.save == "show":
	plt.show()
else:
	plt.savefig(args.save, dpi = fig.dpi)