# pull_data.tcl
# by Umair Khan
# Reichow Lab

# Monitor the following interactions throughout the given trajectory:
#  - inter-atomic distance between CZ/CG (R9/N9) of chain i and CZ/CG (R9/N9) of chain i+1
#  - dihedral angle for C-CA-CB-CG (R9/N9)
#  - inter-atomic distance between CZ/CG (R9/N9) of chain i and CD (E12) of chain i+1

# Command-line arguments (in order):
#  - filename to write data to
#  - chain configuration (0 for Cx46 base, 1 for Cx50 base)
#  - ninth residue identifier (0 for R9, 1 for N9)

# Each line of the file is formatted as follows:
# "frame from_chain [R/N]9_C[Z/G]_i+1 dihed E12_CD_i+1
#  - frame	           -- the frame of the trajectory
#  - from_chain        -- the chain being measured from
#  - [R/N]9_C[Z/G]_i+1 -- the distance to the CZ/CG (R9/N9) of the next chain
#  - dihed             -- the C-CA-CB-CG (R9/N9) dihedral angle
#  - E12_CD_i+1        -- the distance to the CD (E12) of the next chain


# Pull the data from the list of given chains.
proc get_data {frame outfile chain_indices chains ninth_res} {

	# Go through given chains
	foreach index $chain_indices {

		# Get chains under consideration
		set this_chain [lindex $chains $index]
		set next_chain [lindex $chains [expr ($index + 1) % 6]]

		# Get atom indices
		set c [[atomselect top "protein and chain $this_chain and resid 9 and name C"] get index]
		set ca [[atomselect top "protein and chain $this_chain and resid 9 and name CA"] get index]
		set cb [[atomselect top "protein and chain $this_chain and resid 9 and name CB"] get index]
		set cg [[atomselect top "protein and chain $this_chain and resid 9 and name CG"] get index]
		set cz [[atomselect top "protein and chain $this_chain and resid 9 and name CZ"] get index]
		set next_cg [[atomselect top "protein and chain $next_chain and resid 9 and name CG"] get index]
		set next_cz [[atomselect top "protein and chain $next_chain and resid 9 and name CZ"] get index]
		set next_e12_cd [[atomselect top "protein and chain $next_chain and resid 12 and name CD"] get index]

		# Get data
		if {$ninth_res == 0} {
			set ninth_dist [measure bond [list $cz $next_cz]]
			set e12_dist [measure bond [list $cz $next_e12_cd]]
		} else {
			set ninth_dist [measure bond [list $cg $next_cg]]
			set e12_dist [measure bond [list $cg $next_e12_cd]]
		}
		set dihed [measure dihed [list $c $ca $cb $cg]]

		# Write data
		puts $outfile "$frame\t$this_chain\t$ninth_dist\t$dihed\t$e12_dist"

	}

}


# Main running procedure.
proc run {filename chain_conf ninth_res} {

	# Get the correct chain configuration
	if {$chain_conf == 0} {
		puts "Selected Cx46 base chain configuration."
		set chain_top {A B C D E F}
		set chain_bot {G H I J K L}
	} elseif {$chain_conf == 1} {
		puts "Selected Cx50 base chain configuration."
		set chain_top {A B C D E H}
   		set chain_bot {G I J K L F}
	} else {
		puts "Invalid chain configuration."
		return
	}

	# Check the residue identifier
	if {$ninth_res == 0} {
		puts "Ninth residue is R."
	} elseif {$ninth_res == 1} {
		puts "Ninth residue is N."
	} else {
		puts "Invalid ninth residue identifier."
		return
	}

	# Confirm with user
	puts -nonewline "Proceed? "
	flush stdout
	set confirm [gets stdin]
	if {"$confirm" ne "y"} { return }

	# Set up basic variables
	set outfile [open $filename w]
	set num_frames [molinfo top get numframes]
	set chain_indices {0 1 2 3 4 5}

	# Write fields to first line of file
	if {$ninth_res == 0} {
		puts $outfile "frame\tfrom_chain\tR9_CZ_i+1\tdihed\tE12_CD_i+1"
	} else {
		puts $outfile "frame\tfrom_chain\tN9_CG_i+1\tdihed\tE12_CD_i+1"
	}

	# Go through frames and get data from each hemisphere
	for {set frame 0} {$frame < $num_frames} {incr frame} {
		animate goto $frame
		get_data $frame $outfile $chain_indices $chain_top $ninth_res
		get_data $frame $outfile $chain_indices $chain_bot $ninth_res
		if {$frame % 100 == 0} { puts "Finished frame $frame." }
	}

	# Cleanup
	close $outfile

}
