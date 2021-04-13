#
#			Program	: Edge_Eraser.py
#			Author	: Bassam Haddad
#

import numpy as np

def edge_erase(pop_mat,bin_size=1,cutoff=10):
	"""
	Due to the periodic boundaries of a molecular simulation, one sees ion transition events going from the top of the box to the bottom of the box. When calculating transition
	probabilities, we find that there are some spurious probabilities for an ion moving "220" A in a single step (2 ps). This program takes the transition-population matrix and
	erases the counts involving transitions across the periodic boundaries, which are manifested near the edges of the transition matrix which should be nearly diagonal.
	"""
	cutoff		= cutoff / bin_size
	num_row,num_col	= pop_mat.shape
	for i in range (num_row):
		for j in range(num_col):
			if abs(i - j) >= cutoff:
				pop_mat[i,j] = 0
			else:
				pass
	return pop_mat
	"""
	The purpose of tri_diag is to only save the tri-diagonal of the transition
	matrix, such that we are only looking at nearest neighbors. This is more si-
	milar to the data that is in the rates method. The purpose is to compare to
	the rates method to understand why the Pss method gives a symmetric PMF.
	"""
def tri_diag(pop_mat,bin_size):
	num_row,num_col = pop_mat.shape
	for i in range(num_col):
		for j in range(num_row):
			if abs(i - j) > 1:
				# keep the connection between periodic boudaries
				if abs(i - j) == (num_row - 1):
					pass
				else:
					pop_mat[j,i] = 0
			else:
				pass
	return pop_mat
