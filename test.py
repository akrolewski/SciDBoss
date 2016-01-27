mean_spec = continuum

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
	zabs = wavelength/lya - 1.
	fz_regrid = fz[np.argmin(np.abs(zabs[:,np.newaxis]-fz[:,0]),axis=1),1]
	# Technically fz_regrid should have some kind of extrapolation to z > zmax.  since we have to include quasars
	# with z > zmax (in order to measure the forest up to zmax) but then we're also using these quasars to compute the continuum
	fz_regrid[rest_wavelength > lya] = 1.
	weight = 1./((var/(a_qso_ind**2*template**2))+fz_regrid**2*map(lambda x: intrinsic_variance(x, lya), rest_wavelength))
	
	for i in range(len(wavelength)):
		wv = wavelength[i]
		if rest_wavelength[i] < forest_high_limit and rest_wavelength[i] > forest_low_limit:
			ind = np.argmin(np.abs(wv-lya*(1.+fz[:,0])))
			all_fz[fz[ind,0]][0].append(flux[i]/(a_qso_ind*template[i]))
			all_fz[fz[ind,0]][1].append(weight[i])

for i in range(len(fz[:,0])):
	z = fz[i,0]
	ind_fz = all_fz[z]
	#print np.mean(ind_fz[0])
	fz[i,1] = np.sum(np.array(ind_fz[0])*np.array(ind_fz[1]))/np.sum(ind_fz[1])
	fz[i,2] = np.sqrt(1./np.sum(ind_fz[1]))