#
#	Program	: Diffusion_Calc.py
#	Author	: Bassam Haddad
#
#	This program takes the transition matrix that is populated in the Propogator program, row-normalizes it, then uses it to calculate a diffusion coefficient.
#

import numpy as np

def normalize(pop_mat):
	"""
	This function takes the populated transition matrix (i.e. matrix with raw numbers of transitions) and converts it to a 'functionally' probability transition matrix
	by row normalizing the matrix ... in other words, the sum of each row should add to a probability of 1. Consider this, if each row describes the starting position of
	an ion, then it has 3 possible moves: move backwards, stay, move forwards. What is the probability that you do any of the three possibile actions given that you began
	in row i ... the answer is 1. This new matrix (trans_mat) is full of effective transition probabilities.
	"""
  
	column_sum = pop_mat.sum(axis=0)
	colnorm = pop_mat / column_sum
	return colnorm

def Diff_Calc(trans_mat,bin_size,sim_min=-100,sim_max=100,pore_min=-60,pore_max=60):
	"""
	The first set of loops populate a matrix known as 'bin_centers' which gives the center of the bins in units A, and allows for us to use transition matrices with arbitrary
	bin sizes. The next set of loops are performing the actual calculation of the average Diffusion coefficient within the pore. The bin_centers matrix should corrospond to
	the bins in the transition matrix. The values: -110, 110, -60, and 60 are chosen based off of all the PMFs I have made during this project, and are specific to Cx46/50
	or really any connexin gap junction I suppose. Point is, if you have an explicit reason to use a different value, edit it here.
	"""
	num_row,num_col		= trans_mat.shape
	bin_centers		= np.zeros((num_row,1))
	bin_centers[0,0]	= sim_min + (bin_size / 2)
	for i in range(1,num_row,1):
		bin_centers[i,0]	= bin_centers[i-1,0] + bin_size
	Sum,counts,hold,tau = 0,0,0,2e-12									# tau = 2e-12 seconds (2 picoseconds) (each MD frame is written every tau)
	for i in range(0,num_row,1):
		if bin_centers[i,0] >= pore_min and bin_centers[i,0] <= pore_max:
			for j in range(0,num_row,1):
				Sum = Sum + (((bin_centers[i,0] - bin_centers[j,0])**2)	* trans_mat[i,j])	# The sum will equal <dX^2> for bin i
			counts	+=	1
			hold	=	hold + Sum								# This is the sum of all the bin-specific <dX^2> within the pore
			Sum	=	0
		else:
			pass
	Diff_pore	=	(hold) / (2*tau*counts)
	return	Diff_pore
