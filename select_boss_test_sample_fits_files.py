import pyfits
import random
import urllib
import os
import sys
# This script selects a test sample from BOSS DR12 that I can test the fitting with
# It creates a .fits file summarizing the test sample, and downloads the quasar spectra

# Don't run the script if we've already downloaded the test sample
#if os.path.exists('boss_test_sample.fits') or os.path.exists('boss_test_sample/'):
#	print "please don't overwrite the existing files"
#	sys.exit()

file = pyfits.open('DR12Q.fits')
all_qso = file[1].data
z_type = 'Z_VI' # this could be Z_VI (visual inspection redshift), Z_PIPE (pipeline redshift),
# or Z_PCA (pca redshift).  think more about which one to use later!
z_low = 2.15 # set by the red end of spectrograph
z_high = 5 # no good reason for upper bound...

# Select only quasars within the redshift range, and quasars that don't have a BAL
# may have other cuts later...
qso_selection = all_qso[(all_qso[z_type] > z_low) & (all_qso[z_type] < z_high) & (all_qso['BAL_FLAG_VI'] == 0)]

# randomly select 5000 quasars
n_test = 5000
inds = random.sample(xrange(len(qso_selection)),n_test)
qso_test_sample = qso_selection[inds]

#hdulist = pyfits.HDUList([file[1].header, qso_selection])
#pyfits.writeto('boss_test_sample.fits',qso_test_sample,file[1].header)

qso_dir = 'boss_test_sample_fits'
os.mkdir(qso_dir)

cnt = 0
for qso in qso_test_sample:
	plate = qso['plate']
	fiber = qso['fiberid']
	mjd = qso['mjd']
	strplate = str(plate)
	strfiber = str(fiber)
	while len(strplate) < 4:
		strplate = '0' + strplate
	while len(strfiber) < 4:
		strfiber = '0' + strfiber
	strmjd = str(mjd)
	name = 'spec' +'-' + strplate + '-' + strmjd + '-' + strfiber + '.fits'
	opener = urllib.URLopener()
	#rsync = 'rsync://data.sdss3.org/dr12/boss/spectro/redux/v5_7_0/spectra/' + strplate + '/' + name
	url = 'http://data.sdss3.org/sas/dr12/boss/spectro/redux/v5_7_0/spectra/' + strplate + '/' + name
	#spec = pyfits.open(url)
	#flux = spec[1].data['flux']
	#loglam = spec[1].data['loglam']
	#ivar = spec[1].data['ivar']
	#mask = spec[1].data['and_mask']
	#qsofile = open(qso_dir + '/' + strplate + '-' + strmjd + '-' + strfiber + '.dat','w')
	#print url
	asshole = "wget %s %s/" % (url,qso_dir)
	os.system(asshole)
	cnt = cnt + 1
	if cnt == 50:
		break
	#for i in range(len(flux)):
	#	qsofile.write('%.8f %.8f %.8f %i\n' % (loglam[i],flux[i],ivar[i],mask[i]))