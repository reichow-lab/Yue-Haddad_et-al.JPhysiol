#
#	Program:	process_ion.py
#	Author:		Bassam Haddad
#
#		It is not always possible to load MD simulation trajectories in a single VMD session, and is therefore not possible to
#	pull entire ion trajectories at once. For example, if there are 100 ions in the simulation, but I need to load the simulation in
#	two halves, then I will end up with 200 ion trajectory files. such that ion_firsthalf_13 and ion_secondhalf_13 are the same ion
#	but their trajectories are split into two separate files. This convention makes future applications of DetBa.py problematic
#	specifically the ability to track the permeation of a single ion. Thus it would help to have a 1 file per ion. This program does
#	just that, it appends the trajectories of the "second half" to the "first half" and adjusts the time value such that it is a
#	continuous process.
# process_ion.py
# by Umair Khan
# Reichow Lab

# Slide together ion trajectories with ease.
# Usage: python process_ion.py [number of files] [output prefix] [list of input prefixes]

# Imports
import sys
import os
from tqdm import tqdm

# Get input arguments
num_files = int(sys.argv[1])
output_prefix = sys.argv[2]
input_prefixes = sys.argv[3:]
num_inputs = len(input_prefixes)

# For each input file, create a list keyed with the filename
# that contains each line of the input file as a list, with
# the first entry being the frame (as an int) and the second
# entry being the coordinate (as a str).
#
print("\nReading input files...")
input_dict = {}
for prefix in input_prefixes:
    for suffix in tqdm(range(num_files), ncols = 100, desc = prefix + "[n]"):
        input_dict[prefix + str(suffix)] = []
        lines = [line.rstrip().split() for line in open(prefix + str(suffix), "r")]
        for line in lines:
            input_dict[prefix + str(suffix)].append([int(line[0]), line[1]])

# For each prefix and suffix, get the final frame value of the
# previous prefix with the same suffix and add that value
# to each frame of the current file.
#
print("\nSliding...")
for i in range(1, num_inputs):
    for suffix in tqdm(range(num_files), ncols = 100, desc = input_prefixes[i] + "[n]"):
        offset = input_dict[input_prefixes[i - 1] + str(suffix)][-1][0] + 1
        for line in input_dict[input_prefixes[i] + str(suffix)]:
            line[0] += offset

# For each suffix, combine all the lines from the input prefixes
# as strings and then combine all of those strings, then write
# the result to the output prefix with the suffix.
#
#	where 'xxxx' is the adjustment needed to make to the time column of the "second half" data ... for example, them time in
#	'ion_firsthalf_13' goes from 0 - 10000, and the time in 'ion_secondhalf_13' also goes from 0 - 10000, the slide amount (xxxx)
#	will be 10000, such that your resulting file, 'ion_total_13', goes from 0 - 20000.
#
#	The prefix refers to the string up to the number, using our prior examples, the prefix for 'ion_firsthalf_13' would be
#	'ion_firsthalf_'.

from tqdm	import tqdm
from sys        import argv
import os

script, slide_amount = argv

num_files               =       int(input("How many input files are there? "))

prefix			=	str(input("What is the prefix of the files you want to slide? "))

for c in tqdm(range(0,num_files,1)):

	dat = str(prefix) + str(c)

	Data = open(dat, 'r')

	new = 'temp_' + str(c)

	newD = open(new, 'w')

	for line in Data:

		values = line.split()

		pos = float(values[0]) + float(slide_amount)

		newD.write(str(pos) + ' ' + str(values[1]) + '\n')

	Data.close()
	newD.close()

prefix_init		=	str(input("What is the prefix of the leading file? "))
prefix_new		=	str(input("What is the prefix of your new file? "))

for n in tqdm(range(0,num_files,1)):

	outfile = str(prefix_init) + str(n)

	outnew	= str(prefix_new) + str(n)

	oldfile = str(prefix) + str(n)

	infile	= 'temp_' + str(n)

	fout	= open(outfile, 'a')

	fin	= open(infile, 'r')

	data	= fin.read()

	fin.close()

	fout.write(data)

	fout.close()

	os.rename(outfile, outnew)

	os.remove(infile)

	os.remove(oldfile)

print("\nWriting output files...")
for suffix in tqdm(range(num_files), ncols = 100, desc = output_prefix + "[n]"):
    split_files = []
    for prefix in input_prefixes:
        split_files.append("\n".join(["\t".join([str(line[0]), line[1]]) for line in input_dict[prefix + str(suffix)]]))
    with open(output_prefix + str(suffix), "w") as f:
        f.write("\n".join(split_files))

# For each input prefix and suffix combination, delete the file.
#
rm_input = input("\nWould you like to delete the input files? (y/n) ")
if rm_input == "y":
    print("\nDeleting input files...")
    for prefix in input_prefixes:
        for suffix in tqdm(range(num_files), ncols = 100, desc = prefix + "[n]"):
            os.remove(prefix + str(suffix))

# Line break for consistency
print("")
