#
#	Program	: Current_Calculator.py
#	Author	: Bassam Haddad
#		This module of DetBa.py takes a PMF (currently just potassium) and then predicts the Current-Voltage relationship (I-V curve) for the given channel. The work is based primarily of the methods
#	from (Zonta et al. 2014). It calculates a mean first passage time (Tau) for an ion traversing the pore, both 'forward' and 'backward'. Taking the inverse of Tau, one gets a rate of transition (1/s),
#	which can then be used to calculate a current:
#
#		I = q[Tau_f - Tau_r] : I has units picoamperes (pA), which is defined as picocoulombs per second (pC/s), and q is the charge that is moving in units pC.
#
#	Then, assuming that a potential difference is a linear function added to a PMF (two terms that do work on the ion), one can calculate the Taus for the channel under a theoretically applied voltage:
#
#		U(x) = W(x) + q(dV)/L * (L-x) : U(x) is the PMF with voltage (dV) applied linearly ((L-x)/L) to the original PMF W(x). Here, q has units 'e' and just gives the sign of the charge. dV is in units
#						Kcal/(mol*e).
#
#	Thus the current (I) can be calculated for each applied voltage (mV).
#
#	Constants Used:
#
#		DifC	= 1.957e11	: (Angstrom)^2/s				: This value was experimentally determined for K+ in bulk water (Samson et al. 2003)
#		q	= 1.60217663e-7	: elementary charge of proton in 'picocoulombs'
#
#	Methods:
#
#		Current	:
#
#			Inputs	: PMF, num_bins, bin_size, [Diffusion coefficient (DifC), and charge of ion are defaulted]	: Eventually DifC will be calculated from the transition matrix.
#			Outputs	: Current (I) in picoamperes, Tau_forward/reverse in seconds
#
#		VoltPMF	:
#
#			Inputs	: PMF (equilibrium), dV (voltage in mV), [ion_sign, and charge of ion are defaulted]		: Will need it to be negative for Cl- PMFs, however they are a minor species.
#			Outputs	: PMF (voltage applied)
#
#		Text2PMF:
#
#			Inputs	: Pre-prepared PMF in a text file.
#			Outputs	: PMF in an array, num_bins

import numpy as np
R	= 0.00198588		# Ideal gas constant in units Kcal/molK
T	= 310			# Temperature in units K

def Current(PMF,num_bins,bin_s,DifC=1.957e11,q=1.60217663e-7):
	tau_forward = 0
	tau_reverse = 0
	for y in range(0,num_bins,1):
		for z in range(0,y+1,1):
			hold		= np.exp((PMF[y,1] - PMF[z,1])/(R*T))*(bin_s**2)
			tau_forward	= tau_forward + hold
	for y in range(num_bins-1,-1,-1):
		for z in range(num_bins-1,y-1,-1):
			hold		= np.exp((PMF[y,1] - PMF[z,1])/(R*T))*(bin_s**2)
			tau_reverse	= tau_reverse + hold
	tau_forward	= tau_forward/DifC
	tau_reverse	= tau_reverse/DifC
	K_forward	= 1/tau_forward
	K_reverse	= 1/tau_reverse
	Current		= q*(K_forward - K_reverse)
	return (Current,tau_forward,tau_reverse)
def Text2PMF(PMF,n=0):
	num_bins	= sum(1 for line in open(PMF))
	PMF_z		= np.zeros((num_bins,2))
	pmfile		= open(PMF)
	for dataline in pmfile:
		val		= dataline.split()
		if val[0] == "Pore":
			pass
		else:
			PMF_z[n,1]	= float(val[1])
			PMF_z[n,0]	= float(val[0])
		n		+= 1
	return (PMF_z,num_bins)
def VoltPMF(PMF,dV,num_bins,ion_sign=1,q=1.60217663e-7):
	if ion_sign == 1:
		q	= q
	elif ion_sign == -1:
		q	= -q
	U_x	= np.zeros((num_bins,2))	# Initialize new voltage applied PMF array
	dV_e	= float(dV) * (0.001) * 23.061	# Conversion of mV (dV) to Kcal/mol (dV_e)
	for x in range(0,num_bins,1):
		U_x[x,1]	= PMF[x,1] + ((dV_e)/(num_bins))*(num_bins - x)
		U_x[x,0]	= PMF[x,0]
	return U_x
