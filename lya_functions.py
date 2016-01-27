import numpy as np
import pyfits
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


def intrinsic_variance(wavelength, ly_alpha):
	# Returns the intrinsic variance for the quasar: see eqn. 3.3
	# in doi:10.1088/1475-7516/2012/11/059 (Font-Ribera DLA paper)
	var_red_of_ly_alpha = 0.01
	var_blue_of_ly_alpha = 0.1
	# These are rough values that can be adjusted later
	# ly_alpha = 1215.24
	if wavelength > ly_alpha:
		return var_red_of_ly_alpha
	else:
		return var_blue_of_ly_alpha

#@do_profile()
def load_data(plate, fiber, mjd, qso_dir):
	strplate = str(plate)
	strfiber = str(fiber)
	while len(strplate) < 4:
		strplate = '0' + strplate
	while len(strfiber) < 4:
		strfiber = '0' + strfiber
	strmjd = str(mjd)
	spectrum = np.loadtxt(qso_dir + '%s-%s-%s.dat' % (strplate, strmjd, strfiber))
	return spectrum


def load_all_data(qso_list, qso_dir):
	hdulist = pyfits.open(qso_list)
	boss_test_sample = hdulist[1].data
	all_qso = {}

	for qso in boss_test_sample:
		plate = qso['plate']
		fiber = qso['fiberid']
		mjd = qso['mjd']
		z = qso['z_vi']
		strplate = str(plate)
		strfiber = str(fiber)
		while len(strplate) < 4:
			strplate = '0' + strplate
		while len(strfiber) < 4:
			strfiber = '0' + strfiber
		strmjd = str(mjd)
		spectrum = np.loadtxt(qso_dir + '%s-%s-%s.dat' % (strplate, strmjd, strfiber))
		all_qso[(plate, fiber, mjd)] = [spectrum, z]

	return all_qso


def mask_data(spectrum):
	# Mask data where any maskbit is set, and where ivar == 0
	loglam = spectrum[:,0][(spectrum[:,3] == 0) & (spectrum[:,2] != 0)]
	flux = spectrum[:,1][(spectrum[:,3] == 0) & (spectrum[:,2] != 0)]
	ivar = spectrum[:,2][(spectrum[:,3] == 0) & (spectrum[:,2] != 0)]

	# Remove regions around particularly prominent contaminating lines, which often 
	# aren't properly subtracted by SDSS
	dx = 6 # Exclude 6 angstroms around each emission line
	bkg_lines = [5577, 5890, 6300, 6363] 
	# night sky line, sodium D (interstellar medium), night sky line, night sky line
	for line in bkg_lines:
		flux2 = flux[(10**loglam < (line-dx)) | (10**loglam > (line+dx))]
		ivar2 = ivar[(10**loglam < (line-dx)) | (10**loglam > (line+dx))]
		loglam2 = loglam[(10**loglam < (line-dx)) | (10**loglam > (line+dx))]
		flux = flux2
		ivar = ivar2
		loglam = loglam2
	#print 1/0
	return np.transpose(np.array([loglam2, flux2, ivar2]))
	
def fnorm(all_qso, norm_range_low, norm_range_high, lya):
	fnorm = np.zeros(len(all_qso.values()))
	for i in range(len(all_qso.values())):
		qso = all_qso.values()[i]
		spectrum = qso[0]
		z = qso[1]
		spectrum = mask_data(spectrum)
		wavelength = 10**spectrum[:,0]
		flux = spectrum[:,1]
		var = 1./spectrum[:,2]
		rest_wavelength = wavelength/(1.+z)
		weight = 1./(var+flux**2*map(lambda x: intrinsic_variance(x, lya), rest_wavelength))
		fnorm[i] = np.sum(flux[(rest_wavelength > norm_range_low) 
		& (rest_wavelength < norm_range_high)]*weight[(rest_wavelength > norm_range_low) 
		& (rest_wavelength < norm_range_high)])/np.sum(weight[(rest_wavelength > norm_range_low) 
		& (rest_wavelength < norm_range_high)])
	return fnorm
		
def create_fz_grid(wavelength, rest_wavelength, fz, lya, zmax):
	zabs = wavelength/lya - 1.
	fz_regrid = fz[np.argmin(np.abs(zabs[:,np.newaxis]-fz[:,0]),axis=1),1]
	# For z > zmax, we have no measurement of the continuum.  use the Faucher-Giguere estimate.
	fz_regrid[zabs > zmax] = np.exp(-0.0018*(1+zabs[zabs > zmax])**3.92)
	# Technically fz_regrid should have some kind of extrapolation to z > zmax.  since we have to include quasars
	# with z > zmax (in order to measure the forest up to zmax) but then we're also using these quasars to compute the continuum
	fz_regrid[rest_wavelength > lya] = 1.	
	return fz_regrid


def normalize(rest_wavelength, flux, var, lya):
	norm_range_low = 1275.
	norm_range_high = 1285.
	weight = 1./((var/(flux**2))+map(lambda x: intrinsic_variance(x, lya), rest_wavelength))
	normalization = np.sum(flux[(rest_wavelength > norm_range_low) 
		& (rest_wavelength < norm_range_high)]*weight[(rest_wavelength > norm_range_low) 
		& (rest_wavelength < norm_range_high)])/np.sum(weight[(rest_wavelength > norm_range_low) 
		& (rest_wavelength < norm_range_high)])
	# normalization = np.nanmean(flux[(rest_wavelength > norm_range_low) 
	#	& (rest_wavelength < norm_range_high)])
	# Right now, if there is no data in the 1275 to 1285 A range, the spectrum is 
	# ignored for computing the mean
	# if we want to fold in information from lower-redshift quasars, 
	# then we will have to change this step
	if np.isnan(normalization):
		return np.nan
	else:
		return normalization
