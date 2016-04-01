import pyfits
import numpy as np
import matplotlib.pyplot as plt
from lya_functions import *
from mean_spectrum import *
import time
from line_profiler import LineProfiler

def do_profile(follow=[]):
	def inner(func):
		def profiled_func(*args, **kwargs):
			try:
				profiler = LineProfiler()
				profiler.add_function(func)
				for f in follow:
					profiler.add_function(f)
				profiler.enable_by_count()
				return func(*args, **kwargs)
			finally:
				profiler.print_stats()
		return profiled_func
	return inner

def chisq_minimization(flux,template,weight):
	# This gives a result equivalent to the following quantity:
	# np.sum((flux/template)*weight)/np.sum(weight)
	# note that we need to divide the weight by the template**2
	# because the proper weights are the variance about F, NOT the variance about F/C
	# As a result, using chisq minimization is EXACTLY equivalent to the weighted means
	# that we use in the other parts of the iteration
	# I write this weighted mean as a chisq minimization to easily allow for generalization
	# to a multi-parameter continuum model
	weight = np.diag(np.sqrt(weight/(template**2)))
	flux = np.reshape(flux,(len(flux),1))
	template = np.reshape(template,(len(template),1))
	return np.linalg.lstsq(np.dot(weight,template), np.dot(weight,flux))

def weighted_mean(flux,template,weight):
	# this is equivalent to chisq minimization, using the following 
	return np.sum((flux/template)*weight)/np.sum(weight)

#@do_profile()
def find_fz(all_qso, a_qso, mean_spec, test_fz, forest_low_limit, forest_high_limit, lya, zmax):
	fz = test_fz
	all_fz = {}
	for z in fz[:,0]:
		all_fz[z] = [[],[]]
	all_chisq_dof = np.zeros(len(all_qso.values()))
	all_aqso = np.zeros(len(all_qso.values()))

	for i in range(len(all_qso.values())):
		qso = all_qso.values()[i]
		spectrum = qso[0]
		z = qso[1]
		spectrum = mask_data(spectrum)
		wavelength = 10**spectrum[:,0]
		flux = spectrum[:,1]
		var = 1./spectrum[:,2]
		rest_wavelength = wavelength/(1.+z)
		template = mean_spec[np.argmin(np.abs(rest_wavelength[:,np.newaxis]-mean_spec[:,0]),axis=1),1]
	
		# First step: find A_qso (possibly other parameters) for each quasar
		#print len(var), len(template), len(rest_wavelength)
		#print 1./(var/(template**2))
		#print a_qso[i]**2*np.array(map(lambda x: intrinsic_variance(x, lya), rest_wavelength))
		weight = 1./((var/(template**2)) + a_qso[i]**2*np.array(map(lambda x: intrinsic_variance(x, lya), rest_wavelength)))
		# weight is the inverse variance of (F/C)
		flux_red_of_lya = flux[rest_wavelength > lya]
		template_red_of_lya = template[rest_wavelength > lya]
		weight_red_of_lya = weight[rest_wavelength > lya]
		#a_qso = weighted_mean(flux_red_of_lya,template_red_of_lya,weight_red_of_lya)
		chisq_output = chisq_minimization(flux_red_of_lya,template_red_of_lya,weight_red_of_lya)
		a_qso_ind = chisq_output[0][0][0]
		chisq = chisq_output[1][0]
		all_aqso[i] = a_qso_ind
		all_chisq_dof[i] = chisq/len(flux_red_of_lya)
	
		# Second step: find <F(z)>
		#zabs = wavelength/lya - 1.
		#fz_regrid = fz[np.argmin(np.abs(zabs[:,np.newaxis]-fz[:,0]),axis=1),1]
		# For z > zmax, we have no measurement of the continuum.  use the Faucher-Giguere estimate.
		#fz_regrid[zabs > zmax] = np.exp(-0.0018*(1+zabs[zabs > zmax])**3.92)
		# Technically fz_regrid should have some kind of extrapolation to z > zmax.  since we have to include quasars
		# with z > zmax (in order to measure the forest up to zmax) but then we're also using these quasars to compute the continuum
		#fz_regrid[rest_wavelength > lya] = 1.
		fz_regrid = create_fz_grid(wavelength, rest_wavelength, fz, lya, zmax)
		weight = 1./((var/(a_qso_ind**2*template**2))+fz_regrid**2*map(lambda x: intrinsic_variance(x, lya), rest_wavelength))
	
		for i in range(len(wavelength)):
			wv = wavelength[i]
			zabs_ind = wv/lya - 1.
			if rest_wavelength[i] < forest_high_limit and rest_wavelength[i] > forest_low_limit and zabs_ind < zmax:
				ind = np.argmin(np.abs(wv-lya*(1.+fz[:,0])))
				all_fz[fz[ind,0]][0].append(flux[i]/(a_qso_ind*template[i]))
				all_fz[fz[ind,0]][1].append(weight[i])

	for i in range(len(fz[:,0])):
		z = fz[i,0]
		ind_fz = all_fz[z]
		#print np.mean(ind_fz[0])
		fz[i,1] = np.sum(np.array(ind_fz[0])*np.array(ind_fz[1]))/np.sum(ind_fz[1])
		fz[i,2] = np.sqrt(1./np.sum(ind_fz[1]))
	
	return [fz, all_aqso, all_chisq_dof]

def recalc_mean_spec(mean_spec, fz, zabs, all_beta):
	all_spec = {}
	for wv in mean_spec[:,0]:
		all_spec[wv] = [[],[]]
	for i in range(len(all_qso.values())):
		qso = all_qso.values()[i]
		spectrum = qso[0]
		z = qso[1]
		spectrum = mask_data(spectrum)
		wavelength = 10**spectrum[:,0]
		flux = spectrum[:,1]
		var = 1./spectrum[:,2]
		rest_wavelength = wavelength/(1.+z)
	
		template = mean_spec[np.argmin(np.abs(rest_wavelength[:,np.newaxis]-mean_spec[:,0]),axis=1),1]

		weight = 1./((var/((all_beta[i][0][0]*template)**2))+map(lambda x: intrinsic_variance(x, lya), rest_wavelength))
		norm_flux = normalize(rest_wavelength,flux)
		#print weight

		if np.all(np.isnan(norm_flux)):
			continue

		forest_low_limit = 1041.
		forest_high_limit = 1185.
		lya = 1215.24
		inds = np.argmin(np.abs(rest_wavelength[:,np.newaxis]-mean_spec[:,0]),axis=1)
		for i in range(len(rest_wavelength)):
			wv = rest_wavelength[i]
			obs_wv = wavelength[i]
			ind = inds[i]
			if wv > forest_low_limit and wv < forest_high_limit:
				zind = obs_wv/lya-1.
				fz_ind = fz[np.argmin(np.abs(zind-zabs))]
				all_spec[mean_spec[:,0][ind]][0].append(norm_flux[i]/fz_ind)
				all_spec[mean_spec[:,0][ind]][1].append(weight[i])
	#print all_spec
	for i in range(len(mean_spec[:,0])):
		wv = mean_spec[:,0][i]
		spec = all_spec[wv]
		if wv > forest_low_limit and wv < forest_high_limit:
			#print wv
			mean_spec[i,1] = np.sum(np.array(spec[0])*np.array(spec[1]))/np.sum(spec[1])
			
	return mean_spec
