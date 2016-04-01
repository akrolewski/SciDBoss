import pyfits
import numpy as np
import matplotlib.pyplot as plt
from lya_functions import *
from mean_continuum import *
from find_fz import *
from prettyplot import *
import time
import copy
import os
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

t0 = time.time()
qso_list = '/Users/ALEX/Berkeley/Lyman_alpha_forest/boss_test_sample.fits'
qso_dir = '/Users/ALEX/Berkeley/Lyman_alpha_forest/boss_test_sample/'

all_qso = load_all_data(qso_list, qso_dir)
print time.time()-t0

lya = 1215.67

deltaz = 0.03

begin_bin = 3566.97 # currently where the spectra begin.  this may change if we decide
# to mask some pixels at the very blue-most end of the spectrograph
zmin = begin_bin/lya - 1. + deltaz/2.
# increase by deltaz/2 because I want the bin EDGE to be at (3566.97/1215.24)-1.

# this "begin_bin" is not exactly correct.  the "begin_bin" is different for different plates
# so I really should just mask a bunch of pixels at the beginning of the spectrum

min_wv = 600.
max_wv = 3600.
wv_bin = 1.
forest_low_limit = 1041.
forest_high_limit = 1185.

norm_range_low = 1275.
norm_range_high = 1285.

zwidth = 1.5
zmax = zmin + zwidth

iteration = 0
stat = 2.*(zwidth/deltaz)

#z_pivot = 2.2
znorm_low = 2.2
znorm_high = 2.6

blue_lambda_low = 1045.
blue_lambda_high = 1055.
blue_amp = 1.25

continuum = np.transpose([np.linspace(min_wv,max_wv,1+int((max_wv-min_wv)/wv_bin)),np.ones(1+int((max_wv-min_wv)/wv_bin))])
a_qso = fnorm(all_qso, norm_range_low, norm_range_high, lya)

# Ignore quasars with no data in 1275-1285 A range
inds = [i for i in range(len(a_qso)) if not np.isnan(a_qso[i])]
new_all_qso = {}
keys = np.array(all_qso.keys())[inds]
vals = np.array(all_qso.values())[inds]
for i in range(len(vals)):
	new_all_qso[tuple(keys[i])] = list(vals[i])
all_qso = new_all_qso
a_qso = a_qso[inds]

#fz = np.transpose([np.arange(zmin,zmax,deltaz), np.zeros(1+int((zmax-zmin)/deltaz)), np.zeros(1+int((zmax-zmin)/deltaz))])
fz = np.transpose([np.linspace(zmin,zmax,1+int((zmax-zmin)/deltaz)),np.zeros(1+int((zmax-zmin)/deltaz)),np.zeros(1+int((zmax-zmin)/deltaz))])
fz[:,1] = np.exp(-0.0018*(1+fz[:,0])**3.92)

def normalize_amplitude(wv, fl, norms):
	flux_norm = np.mean(fl[(wv > norms[0]) & (wv < norms[1])])
	return fl/flux_norm
	
def normalize_spectral_slope_blue_side(wv, fl, lya, norms, lambda_norms, blue_amp):
	blue_amp_meas = np.mean(fl[(wv > lambda_norms[0]) & (wv < lambda_norms[1])])
	fix_slope = np.log(blue_amp)/np.log(np.mean(lambda_norms)/np.mean(norms))
	meas_slope = np.log(blue_amp_meas)/np.log(np.mean(lambda_norms)/np.mean(norms))
	newfl = copy.copy(fl)
	newfl[wv < np.mean(norms)] = fl[wv < np.mean(norms)]*(np.mean(norms)**(meas_slope-fix_slope))*(wv[wv < np.mean(norms)]**(fix_slope-meas_slope))
	return newfl

while stat > ((zmax-zmin)/deltaz):
	continuum = mean_continuum(all_qso, a_qso, continuum, fz, lya, zmax)
	
	continuum[:,1] = normalize_amplitude(continuum[:,0], continuum[:,1], [norm_range_low,norm_range_high])
	continuum[:,1] = normalize_spectral_slope_blue_side(continuum[:,0], continuum[:,1],lya, [norm_range_low,norm_range_high],[blue_lambda_low,blue_lambda_high],blue_amp)
	
	
	#flux_norm = np.mean(continuum[(continuum[:,0] > norm_range_low) & (continuum[:,0] < norm_range_high),1])
	#continuum[:,1] = continuum[:,1]/flux_norm
	
	output = find_fz(all_qso, a_qso, continuum, fz, forest_low_limit, forest_high_limit, lya, zmax+deltaz/2.)
	fz = output[0]
	a_qso = output[1]
	chisq_dof = output[2]
	
	fz_faucher_giguere = np.exp(-0.0018*(1+fz[:,0])**3.92)
	fz_norm = np.mean(fz_faucher_giguere[(fz[:,0] > znorm_low) & (fz[:,0] < znorm_high)])
	fz_val = np.mean(fz[(fz[:,0] > znorm_low ) & (fz[:,0] < znorm_high),1])
	fz[:,1] = fz[:,1]*(fz_norm/fz_val)
	
	if iteration > 0:
		stat = np.nansum(((fz[:,1]-orig_fz[:,1])/fz[:,2])**2)
	
	# plots
	prettyplot()
	plt.plot(continuum[:,0],continuum[:,1],color='b')
	plt.savefig("continuum_%i.pdf" % iteration)
	
	prettyplot()
	plt.errorbar(fz[:,0],fz[:,1],yerr=fz[:,2],color='b')
	plt.savefig("fz_%i.pdf" % iteration)
	
	# write to text file
	mean_spec_file = open('continuum%i.dat' % iteration,'w')
	for i in range(len(continuum)):
		mean_spec_file.write('%f %f\n' % (continuum[i,0],continuum[i,1]))
	mean_spec_file.close()
		
	fz_file = open('fz%i.dat' % iteration,'w')
	fz_file.write('stat: %f\n' % stat)
	for i in range(len(fz)):
		fz_file.write('%f %f %f\n' % (fz[i,0],fz[i,1],fz[i,2]))
	fz_file.close()
	
	orig_fz = copy.copy(fz)
	iteration = iteration + 1
	
	# Stuff to do:
	# plots for each quasar
	# Add in functionality for restricted wavelength ranges


# prettyplot()
# plt.plot(fz[:,0],np.exp(-0.0018*(1+fz[:,0])**3.92),color='r',label='Faucher-Giguere 2008')
# plt.plot(fz[:,0],fz[:,1],color='b',label='data')
# plt.ylabel("<F(z)>")
# plt.xlabel("z")
# plt.legend()
# plt.savefig("fz_simple_continuum_16_01_19.pdf")
# plt.clf()
# 
# prettyplot()
# plt.plot(mean_spec[:,0],mean_spec[:,1])
# plt.xlim([600.,3000.])
# plt.ylim([0.,4.])
# plt.savefig("mean_spec_16_01_19.pdf")
# 
# def chisq_minimization(flux,template,weight):
# 	flux = np.reshape(flux,(len(flux),1))
# 	template = np.reshape(template,(len(template),1))
# 	weight = np.diag(np.sqrt(weight))
# 	return np.linalg.lstsq(np.dot(weight,template), np.dot(weight,flux))[0][0][0]
# 
# def weighted_mean(flux,template,weight):
# 	return np.sum((flux/template)*weight)/np.sum(weight)
# 	
# def weighted_mean2(flux,template,weight):
# 	return np.sum(template**2*(flux/template)*weight)/np.sum(template**2*weight)
# 
dir = "quasars/"
os.mkdir(dir)
mean_spec = continuum
#test_amps = np.arange(1,5,0.01)
#chisqs = np.zeros(len(test_amps))
inds = chisq_dof.argsort()
chisq_dof_sorted = chisq_dof[inds]
a_qso_sorted = a_qso[inds]
allval = np.array(all_qso.values())
data = allval[:,0]
redshift = allval[:,1]
keys = np.array(all_qso.keys())
data = data[inds]
redshift = redshift[inds]
keys = keys[inds]
for i in range(len(all_qso.values())):
	#if i == 2000 or i == 3000:
	#if i == 3000:
	spectrum = data[i]
	z = redshift[i]
	beta = a_qso_sorted[i]
	spectrum = mask_data(spectrum)
	wavelength = 10**spectrum[:,0]
	flux = spectrum[:,1]
	var = 1./spectrum[:,2]
	rest_wavelength = wavelength/(1.+z)
	#weight = 1./((var/(flux**2))+map(lambda x: intrinsic_variance(x, lya), rest_wavelength))
	template = mean_spec[np.argmin(np.abs(rest_wavelength[:,np.newaxis]-mean_spec[:,0]),axis=1),1]
	zabs = wavelength/lya - 1.
	fz_regrid = fz[np.argmin(np.abs(zabs[:,np.newaxis]-fz[:,0]),axis=1),1]
	fz_regrid[zabs > zmax] = np.exp(-0.0018*(1+zabs[zabs > zmax])**3.92)
	fz_regrid[rest_wavelength > lya] = 1.
	#flux_red_of_lya = flux[rest_wavelength > lya]
	#template_red_of_lya = template[rest_wavelength > lya]
	#weight = 1./(var+flux**2*map(lambda x: intrinsic_variance(x, lya), rest_wavelength))
	#weight_red_of_lya = weight[rest_wavelength > lya]
	#weight = 1./(var/(template**2)+flux**2*map(lambda x: intrinsic_variance(x, lya), rest_wavelength)/(template**2))
	#weight_red_of_lya = weight[rest_wavelength > lya]
	plt.clf()
	prettyplot()
	plt.plot(rest_wavelength,flux,color='b')
	plt.plot(rest_wavelength,beta*template,color='r')
	plt.plot(rest_wavelength,fz_regrid*beta*template,'r--')
	print chisq_dof_sorted[i]
	plt.savefig(dir+'%i_%i-%i-%i.pdf' % (i, keys[i][0],keys[i][1],keys[i][2]))
	#plt.xlim([1000,2200])
	#plt.ylim([-2,6])
	#print all_qso.keys()[i][0]
	#plt.savefig("ex_2_16_01_19.pdf")
	#print all_chisq_dof[i], all_aqso[i]
		#for i in range(len(test_amps)):
		#	amp = test_amps[i]
		#	chisqs[i] = np.sum((flux_red_of_lya-amp*template_red_of_lya)**2*weight_red_of_lya)
		
prettyplot()
plt.plot(fz[:,0],fz[:,1],color='b')
plt.plot(fz[:,0],np.exp(-0.0018*(1+fz[:,0])**3.92),color='r')
plt.savefig("fz_with_fg.pdf")

aqso_chisq_file = open('aqso_chisq.dat','w')
aqso_chisq_file.write('A_qso chisq\n')
for i in range(len(a_qso)):
	aqso_chisq_file.write('%f %f\n' % (a_qso[i],chisq_dof[i]))

aqso_chisq_file.close()

prettyplot()
plt.plot(continuum[:,0],continuum[:,1],color='b')
plt.xlim([1000,1200])
plt.ylim([0,3])
plt.savefig("continuum_zoom_in.pdf")
# 
# print chisq_minimization(flux_red_of_lya,template_red_of_lya,weight_red_of_lya)
# print test_amps[np.argmin(chisqs)]
# print weighted_mean(flux_red_of_lya,template_red_of_lya,weight_red_of_lya)
# # 
# # plt.hist(all_aqso,bins=100)
# # plt.xlim([0,10])
# # plt.xlabel("Aqso")
# # plt.ylabel("Counts")
# # plt.savefig("aqso_hist.pdf")
# # 
# # redshifts = map(lambda x: x[1], all_qso.values())
# # plt.hist(redshifts,bins=100)
# # plt.xlabel("z")
# # plt.ylabel("counts")
# # plt.savefig("redshift_hist.pdf")
# # 
# # plt.scatter(redshifts,all_aqso)
# # plt.xlabel("z")
# # plt.ylabel("Aqso")
# # plt.savefig("z_aqso.pdf")
# # 
# # plt.hist(all_chisq_dof,bins=500)
# # plt.xlim([0,10])
# # plt.xlabel("chisq per dof")
# # plt.savefig("redchisq.pdf")
# # 
# # mean_spec_file = open('final_mean_spec.dat','w')
# # for i in range(len(mean_spec)):
# # 	mean_spec_file.write('%f %f\n' % (mean_spec[i,0],mean_spec[i,1]))
# # 
# # mean_spec_file.close()
# # 
# # amplitude_chisq_file = open('amplitude_chisq.txt','w')
# # amplitude_chisq_file.write('Aqso Chisq/dof\n')
# # for i in range(len(all_aqso)):
# # 	amplitude_chisq_file.write('%f %f\n' % (all_aqso[i], all_chisq_dof[i]))
# # 
# # amplitude_chisq_file.close()
