import pyfits
import numpy as np
from lya_functions import *
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

#@do_profile()
def mean_spectrum(all_qso, mean_spec, lya, iteration=0, fz=None, all_beta=None, 
	forest_low_limit=None, forest_high_limit=None):
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
		norm = normalize(rest_wavelength, flux, var, lya)
		norm_flux = flux/norm
		
		if np.all(np.isnan(norm_flux)):
			continue
		
		if iteration > 0:
			template = mean_spec[np.argmin(np.abs(rest_wavelength[:,np.newaxis]-mean_spec[:,0]),axis=1),1]
			weight = 1./((var/((all_beta[i]*template)**2))+map(lambda x: intrinsic_variance(x, lya), rest_wavelength))
		else:
			weight = 1./((var/(norm**2))+map(lambda x: intrinsic_variance(x, lya), rest_wavelength))

		inds = np.argmin(np.abs(rest_wavelength[:,np.newaxis]-mean_spec[:,0]),axis=1)
		for j in range(len(rest_wavelength)):
			wv = rest_wavelength[j]
			obs_wv = wavelength[j]
			ind = inds[j]
			if iteration > 0:
				if wv > forest_low_limit and wv < forest_high_limit:
					zind = obs_wv/lya-1.
					fz_ind = fz[np.argmin(np.abs(zind-fz[:,0])),1]
					all_spec[mean_spec[:,0][ind]][0].append(norm_flux[j]/fz_ind)
					all_spec[mean_spec[:,0][ind]][1].append(weight[j])
			else:
				all_spec[mean_spec[:,0][ind]][0].append(norm_flux[j])
				all_spec[mean_spec[:,0][ind]][1].append(weight[j])

	for i in range(len(mean_spec[:,0])):
		wv = mean_spec[:,0][i]
		spec = all_spec[wv]
		if iteration > 0:
			if wv > forest_low_limit and wv < forest_high_limit:
				mean_spec[i,1] = np.sum(np.array(spec[0])*np.array(spec[1]))/np.sum(spec[1])
		else:
			mean_spec[i,1] = np.sum(np.array(spec[0])*np.array(spec[1]))/np.sum(spec[1])
		
	return mean_spec