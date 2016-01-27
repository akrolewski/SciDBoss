all_spec = {}
mean_spec = continuum

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
	a = a_qso[i]
	template = mean_spec[np.argmin(np.abs(rest_wavelength[:,np.newaxis]-mean_spec[:,0]),axis=1),1]
	zabs = wavelength/lya - 1.
	fz_regrid = fz[np.argmin(np.abs(zabs[:,np.newaxis]-fz[:,0]),axis=1),1]
	# Technically fz_regrid should have some kind of extrapolation to z > zmax.  since we have to include quasars
	# with z > zmax (in order to measure the forest up to zmax) but then we're also using these quasars to compute the continuum
	fz_regrid[rest_wavelength > lya] = 1.
	weight = 1./((var/(a**2*fz_regrid**2))+template**2*map(lambda x: intrinsic_variance(x, lya), rest_wavelength))
	
	inds = np.argmin(np.abs(rest_wavelength[:,np.newaxis]-mean_spec[:,0]),axis=1)
	for j in range(len(rest_wavelength)):
		wv = rest_wavelength[j]
		obs_wv = wavelength[j]
		#print wv, flux[j]/(a*fz_regrid[j]), weight[j]
		if np.isnan(weight[j]):
			print 5/0
		ind = inds[j]
		all_spec[mean_spec[:,0][ind]][0].append(flux[j]/(a*fz_regrid[j]))
		all_spec[mean_spec[:,0][ind]][1].append(weight[j])
	#print 5/0

for i in range(len(mean_spec[:,0])):
	wv = mean_spec[:,0][i]
	spec = all_spec[wv]
	mean_spec[i,1] = np.sum(np.array(spec[0])*np.array(spec[1]))/np.sum(spec[1])