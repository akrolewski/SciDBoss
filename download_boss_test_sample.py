import pyfits
import random
import urllib
import os
import sys
# Downloads the BOSS test sample

file = pyfits.open('boss_test_sample.fits')
qso_test_sample = file[1].data

qso_dir = 'boss_test_sample'
#os.mkdir(qso_dir)

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
	if os.path.exists(qso_dir + '/' + strplate + '-' + strmjd + '-' + strfiber + '.dat'):
		continue
	else:
		name = 'spec' +'-' + strplate + '-' + strmjd + '-' + strfiber + '.fits'
		opener = urllib.URLopener()
		try:
			url = 'http://data.sdss3.org/sas/dr12/boss/spectro/redux/v5_7_0/spectra/' + strplate + '/' + name
			spec = pyfits.open(url)
		except IOError:
			url = 'http://data.sdss3.org/sas/dr12/boss/spectro/redux/v5_7_2/spectra/' + strplate + '/' + name
			spec = pyfits.open(url)
		print url
		flux = spec[1].data['flux']
		loglam = spec[1].data['loglam']
		ivar = spec[1].data['ivar']
		mask = spec[1].data['and_mask']
		qsofile = open(qso_dir + '/' + strplate + '-' + strmjd + '-' + strfiber + '.dat','w')
		for i in range(len(flux)):
			qsofile.write('%.8f %.8f %.8f %i\n' % (loglam[i],flux[i],ivar[i],mask[i]))