from scidbpy import interface, SciDBQueryError, SciDBArray
import numpy as np
import copy
import time

# Connect to the database
sdb = interface.SciDBShimInterface('http://mndlscidb07-ib:32700')

# Define the lyman alpha wavelength
lya = 1215.67

# Night sky lines
ns1 = 5577
ns2 = 5890
ns3 = 6300
ns4 = 6363
# Remove 6 A on either side around each night sky line
nswidth = 6

# Adjustable parameters for the fitting
# Define the continuum between these points
min_wv = 600.
max_wv = 3600.
wv_bin = 1.
# The intrinsic fractional variance
var_blue_side = 0.1
var_red_side = 0.01
# The range over which we normalize the continuum
norm_range_low = 1275.
norm_range_high = 1285.
# The range that we use to normalize the blue-side continuum slope
blue_lambda_low = 1045.
blue_lambda_high = 1055.
# The amplitude of the continuum at 1050 A (so the continuum is 25% greater here than at 1280 A)
blue_amp = 1.25
# The redshift range for fitting the mean transmission
zmin = 1.5
zmax = 8.1
deltaz = 0.03
# The limits of the Lyman-alpha forest that we use
forest_low_limit = 1041.
forest_high_limit = 1185.
# Normalize to the Faucher-Giguere fit over this range
# so that the mean of <F(z)>in this range matches the F-G mean
z_norm_low = 2.2
z_norm_high = 2.6

t0 = time.time()
# Create the mean continuum array
sdb.query("create array tempcont <lambdar:float> [i=0:%i,?,?]" % int((max_wv-min_wv)/wv_bin))
sdb.query("store(build(tempcont,i*%f+%f),tempcont)" % (wv_bin,min_wv))
sdb.query("store(apply(tempcont,cont_flux,1.),continuum)")
continuum = sdb.wrap_array("continuum")

# Create the mean transmitted flux array
sdb.query("create array tempfz <z:float> [i=0:%i,?,?]" % int((zmax-zmin)/deltaz))
sdb.query("store(build(tempfz,i*%f+%f),tempfz)" % (deltaz,zmin))
sdb.query("store(apply(tempfz,transmission,exp(-0.0018*pow(1+z,3.92))),fz)")
fz = sdb.wrap_array("fz")

print "created fz and continuum arrays"
print time.time()-t0

# Define the weights (eqn. 5 in documentation)
sdb.query("store(apply(boss_test_sample,lambda,pow(10,loglambda),lambdar,pow(10,loglambda)/(1.+redshift),"
"weight,1./(1./ivar+pow(flux,2)*iif(pow(10,loglambda)/(1.+redshift)<%f,%f,%f))),bts2)" % (lya, var_blue_side, var_red_side))
sdb.query("store(apply(bts2,fw,flux*weight),bts3)")

# Remove masked pixels, pixels with ivar = 0,
# and pixels falling inside any of 4 night sky lines
sdb.query("store(filter(bts3,(andmask=0) and (ivar<>0) and"
"((lambdar > %f) and (lambdar < %f)) and"
"((pow(10,loglambda) < %f) or (pow(10,loglambda) > %f)) and"
"((pow(10,loglambda) < %f) or (pow(10,loglambda) > %f)) and"
"((pow(10,loglambda) < %f) or (pow(10,loglambda) > %f)) and"
"((pow(10,loglambda) < %f) or (pow(10,loglambda) > %f))), bts_good_data)"
% (min_wv, max_wv,
ns1-nswidth, ns1+nswidth, ns2-nswidth, ns2+nswidth, ns3-nswidth, ns3+nswidth, ns4-nswidth, ns4+nswidth))

print "removed bad pixels"
print time.time()-t0

# Compute initial guess for A_qso using the average flux in the 1275-1285 A range (eqn. 4 in documentation)
initial_a_qso_chunk_size = 30
sdb.query("store(filter(bts_good_data, (pow(10,loglambda)/(1.+redshift) > %f) and (pow(10,loglambda)/(1.+redshift) < %f)),normrange)"
% (norm_range_low, norm_range_high))
sdb.query("store(uniq(project(normrange,name)),namenum)")
sdb.query("store(redimension(index_lookup(normrange,namenum,normrange.name,numname),<loglambda:float,flux:float,"
"ivar:float,andmask:int64,name:string,redshift:float,weight:double,fw:double> [i=0:*,1048576,0,numname=0:*,%i,0]),fancynormrange)" % initial_a_qso_chunk_size)
sdb.query("store(apply(aggregate(fancynormrange,sum(weight) as w_sum,sum(fw) as fw_sum, numname as name,numname),a_qso,fw_sum/w_sum),fancy_a_qso)")

print "computed initial guess for aqso"
print time.time()-t0

# Mean continuum

# Create attributes telling you the wavelength (lambda_cont) and the index (lambda_index) corresponding to the nearest continuum point
sdb.query("store(apply(bts_good_data, lambda_cont, iif((lambdar-%f)/%f - floor((lambdar-%f)/%f) > 0.5, %f + float(ceil((lambdar-%f)/%f))*%f, %f + float(floor((lambdar-%f)/%f))*%f), "
"z_abs, iif(((lambda/%f-1.)-%f)/%f - floor((((lambda/%f-1.)-%f))/%f) > 0.5, %f + float(ceil(((lambda/%f-1.)-%f)/%f))*%f, %f + float(floor(((lambda/%f-1.)-%f)/%f))*%f)),bts_good_data_mod)"
% (min_wv, wv_bin, min_wv, wv_bin, min_wv, min_wv, wv_bin, wv_bin, min_wv, min_wv, wv_bin, wv_bin,
lya, zmin, deltaz, lya, zmin, deltaz, zmin, lya, zmin, deltaz, deltaz, zmin, lya, zmin, deltaz, deltaz))

sdb.query("store(apply(bts_good_data_mod, lambda_index, iif((lambdar-%f)/%f - floor((lambdar-%f)/%f) > 0.5, ceil((lambdar-%f)/%f), floor((lambdar-%f)/%f)),"
"z_abs_index, iif(((lambda/%f-1.)-%f)/%f - floor((((lambda/%f-1.)-%f))/%f) > 0.5, ceil(((lambda/%f-1.)-%f)/%f), floor(((lambda/%f-1.)-%f)/%f))), bts_test)"
% (min_wv, wv_bin, min_wv, wv_bin, min_wv, wv_bin, min_wv, wv_bin, 
lya, zmin, deltaz, lya, zmin, deltaz, lya, zmin, deltaz, lya, zmin, deltaz))
sdb.query("remove(bts_good_data_mod)")
sdb.query("rename(bts_test, bts_good_data_mod)")

# Add the continuum flux as an attribute
# Use an initial guess of 1 for the continuum flux
sdb.query("store(apply(bts_good_data_mod, cont_flux, 1.), bts_test)")
sdb.query("remove(bts_good_data_mod)")
sdb.query("rename(bts_test, bts_good_data_mod)")


# Add the mean transmission as an attribute
# Use the Faucher-Giguere mean transmission as an initial guess (eqn. 3 in documentation)
sdb.query("store(apply(bts_good_data_mod, real_transmission, iif(lambdar > %f, 1., exp(-0.0018*pow(1+(lambda/%f-1.),3.92)))), bts_test)" % (lya, lya))
sdb.query("remove(bts_good_data_mod)")
sdb.query("rename(bts_test, bts_good_data_mod)")

print "add transmission and continuum as attributes, and add their indices"
print time.time()-t0

# Add the A_qso guess (lines 83-90) as an attribute.
# Begin by redimensioning the array along the index corresponding to the quasar name
sdb.query("store(redimension(index_lookup(bts_good_data_mod,namenum,bts_good_data_mod.name,numname),<lambda:double, flux:float, ivar:float, andmask:int64,"
"name:string, redshift:float, lambdar:double, lambda_cont:double, lambda_index:int64, cont_flux:double, z_abs:double, z_abs_index:int64, real_transmission:double> [i=0:*,1048576,0, numname=0:*,%i,0]),bts_good_data_numname_redim)" % initial_a_qso_chunk_size) #23 s

print "add amplitude as attribute (redimension)"
print time.time()-t0

# Then cross-join with the A_qso guess array to add the amplitude back to the original data
sdb.query("store(redimension(cross_join(bts_good_data_numname_redim, fancy_a_qso, bts_good_data_numname_redim.numname,fancy_a_qso.numname),<lambda:double, flux:float, ivar:float,"
"andmask:int64, name:string, redshift:float, lambdar:double, lambda_cont:double, lambda_index:int64, cont_flux:double, z_abs:double, z_abs_index:int64, real_transmission:double, numname:int64, a_qso:double NULL> [i=0:*,1048576,0]),bts_test)") #~15 s
sdb.query("remove(bts_good_data_mod)")
sdb.query("rename(bts_test, bts_good_data_mod)")

print "add ampliutde as attribute (cross join)"
print time.time()-t0

# Create weights for the continuum flux (eqn. 8 in documentation)
sdb.query("store(apply(bts_good_data_mod, weight_cont, 1./(((1./ivar)/(pow(a_qso,2)*pow(real_transmission,2)))+pow(cont_flux,2)*iif(lambdar <%f,%f,%f)),"
"adj_flux, flux/(a_qso*real_transmission)), bts_test)" % (lya, var_blue_side, var_red_side))
sdb.query("store(apply(bts_test, adj_fw, adj_flux*weight_cont), bts_test2)")
sdb.query("remove(bts_good_data_mod)")
sdb.query("remove(bts_test)")
sdb.query("rename(bts_test2, bts_good_data_mod)")

# Create new continuum flux (eqn. 7 in documentation)
continuum_chunk_size = 1000
sdb.query("store(redimension(bts_good_data_mod,<lambda:double, flux:float, ivar:float, andmask:int64, name:string,"
"redshift:float, lambdar:double, lambda_cont:double, cont_flux:double, z_abs:double, z_abs_index:int64, real_transmission:double, numname:int64,"
"a_qso:double NULL DEFAULT null, weight_cont:double NULL DEFAULT null, adj_flux:double NULL DEFAULT null, adj_fw:double NULL DEFAULT null> [i=0:*,1048576,0,lambda_index=0:*,%i,0]),bts_good_data_lambda_redim)" % continuum_chunk_size)
sdb.query("store(apply(aggregate(bts_good_data_lambda_redim,sum(weight_cont) as w_sum,sum(adj_fw) as fw_sum, lambda_index as index, lambda_index),flux,fw_sum/w_sum,lambda_cont,%f+float(lambda_index)*%f),continuum2)" % (min_wv, wv_bin))

print "create new continuum (redimension plus aggregate)"
print time.time()-t0

# Normalize the continuum to 1 between 1275 and 1285 A
sdb.query("store(aggregate(filter(continuum2, (lambda_cont > %f) and (lambda_cont < %f)), avg(flux) as norm), cont_norm)" % (norm_range_low, norm_range_high))
cont_norm = sdb.wrap_array("cont_norm").toarray()[0]
sdb.query("store(apply(continuum2, flux_norm, flux/%f), continuum2_norm)" % cont_norm)

# Normalize the continuum slope in the blue region to match the continuum slope found from the Paris PCA fits
# This removes the degeneracy between blue-side continuum slope and slope of the mean transmission
sdb.query("store(aggregate(filter(continuum2_norm, (lambda_cont > %f) and (lambda_cont < %f)), avg(flux_norm) as norm), blue_amp_meas)" % (blue_lambda_low, blue_lambda_high))
blue_amp_meas = sdb.wrap_array("blue_amp_meas").toarray()[0]
lambda_norms = [norm_range_low, norm_range_high]
blue_lambda_norms = [blue_lambda_low, blue_lambda_high]
fix_slope = np.log(blue_amp)/np.log(np.mean(blue_lambda_norms)/np.mean(lambda_norms))
meas_slope = np.log(blue_amp_meas)/np.log(np.mean(blue_lambda_norms)/np.mean(lambda_norms))
sdb.query("store(apply(continuum2_norm, cont_flux2, iif(lambda_cont < %f, flux_norm*pow(%f, %f)*pow(lambda_cont, %f),flux_norm)), continuum2_slope_norm)"
% (np.mean(lambda_norms), np.mean(lambda_norms), meas_slope-fix_slope, fix_slope-meas_slope))

# Load the normalized continuum back into the main array
sdb.query("store(redimension(cross_join(bts_good_data_lambda_redim, continuum2_slope_norm, bts_good_data_lambda_redim.lambda_index, continuum2_slope_norm.lambda_index),<lambda:double, flux:float,"
"ivar:float, andmask:int64, name:string, redshift:float, lambdar:double, lambda_cont:double, lambda_index:int64, cont_flux2:double NULL DEFAULT null, z_abs:double, z_abs_index:int64,"
"real_transmission:double, numname:int64, a_qso:double NULL DEFAULT null> [i=0:*,1048576,0]),bts_test)")
sdb.query("remove(bts_good_data_mod)")
sdb.query("rename(bts_test, bts_good_data_mod)")

print "add new continuum (redimension and cross join)"
print time.time()-t0

# Find A_qso in the first iteration

second_a_qso_chunk_size = 30
# Select only points red of lyalpha
# and compute the weight (equation 9 in documentation)
sdb.query("store(apply(bts_good_data_mod, weight, 1./(((1./ivar)/pow(cont_flux2,2))+pow(a_qso,2)*iif(lambdar < %f, %f, %f))), bts_test)"
% (lya, var_blue_side, var_red_side))
sdb.query("store(apply(bts_test, fw, (flux/cont_flux2)*weight), bts_test2)")
sdb.query("remove(bts_test)")
sdb.query("remove(bts_good_data_mod)")
sdb.query("rename(bts_test2, bts_good_data_mod)")

print "compute new weights"
print time.time()-t0

# Compute A_qso (using the weighted mean method which is exactly equivalent to the least-squares method
# explained in and around eqn. 10 of the documentation)
sdb.query("store(redimension(bts_good_data_mod, <lambda:double, flux:float,"
"ivar:float, andmask:int64, name:string, redshift:float, lambdar:double, lambda_cont:double, lambda_index:int64, cont_flux2:double NULL DEFAULT null, z_abs:double, z_abs_index:int64,"
"real_transmission:double, weight:double NULL DEFAULT null,fw:double NULL DEFAULT null> [i=0:*,1048576,0,numname=0:*,%i,0]),bts_good_data_numname_redim2)" % second_a_qso_chunk_size)
sdb.query("store(apply(aggregate(filter(bts_good_data_numname_redim2, lambdar > %f), sum(weight) as w_sum,sum(fw) as fw_sum, numname as name,numname),a_qso2,fw_sum/w_sum),all_a_qso)" % lya)

print "find a_qso (redimension plus aggregate)"
print time.time()-t0

# Cross join to load A_qso back into the main array
sdb.query("store(cross_join(bts_good_data_numname_redim2, all_a_qso, bts_good_data_numname_redim2.numname, all_a_qso.numname),bts_good_data_cj)")

print "load a_qso (cross-join)"
print time.time()-t0

# Compute the chisq and reduced chisq for the A_qso fit
sdb.query("store(apply(aggregate(apply(bts_good_data_cj, ominusesq, ivar*pow(flux-a_qso2*cont_flux2,2)), sum(ominusesq) as chisq, count(ominusesq) as npts, numname as name, numname), red_chisq, chisq/(npts-1.)), all_chisq)")

print "compute chisq"
print time.time()-t0

# Redimension the main array to finish loading A_qso
sdb.query("store(redimension(bts_good_data_cj,<lambda:double, flux:float,"
"ivar:float, andmask:int64, name:string, redshift:float, lambdar:double, lambda_cont:double, lambda_index:int64, cont_flux2:double NULL DEFAULT null, z_abs:double, z_abs_index:int64,"
"real_transmission:double, numname:int64, a_qso2: double NULL DEFAULT null> [i=0:*,1048576,0]),bts_test)") #~15 s
sdb.query("remove(bts_good_data_mod)")
sdb.query("rename(bts_test, bts_good_data_mod)")

print "load a_qso (redimension)"
print time.time()-t0

# Select only the pixels in the Lyalpha forest region, and compute the weight for <F(z)> (eqn. 12 in documentation)
mean_flux_chunk_size = 10000000
sdb.query("store(apply(filter(bts_good_data_mod, (lambdar > %f) and (lambdar < %f)), weight, 1./((1./ivar)/(pow(a_qso2,2)*pow(cont_flux2,2))+pow(real_transmission,2)*iif(lambdar < %f, %f, %f))), bts_forest_region)" % 
(forest_low_limit, forest_high_limit, lya, var_blue_side, var_red_side))
sdb.query("store(apply(bts_forest_region, fw, (flux/(a_qso2*cont_flux2))*weight), bts_forest_region2)")

# Compute <F(z)> (eqn. 11 in documentation)
sdb.query("store(redimension(bts_forest_region2, <lambda:double, flux:float, ivar:float, andmask:int64, name:string, redshift:float, lambdar:double,"
"lambda_cont:double, lambda_index:int64, cont_flux2:double NULL DEFAULT null, z_abs:double, real_transmission:double, numname:int64,"
"a_qso2:double NULL DEFAULT null, weight:double NULL DEFAULT null, fw:double NULL DEFAULT null> [i=0:*,1048576,0,z_abs_index=0:*,%i,0]), bts_forest_region3)" % mean_flux_chunk_size)
sdb.query("store(apply(aggregate(bts_forest_region3, sum(weight) as w_sum, sum(fw) as fw_sum, z_abs_index as name, z_abs_index), transmission2, fw_sum/w_sum), fz2)")

# Normalize <F(z)> to the average value of the Faucher-Giguere <F(z)> from z = 2.2 to z = 2.6
z_norm_low_ind = int((z_norm_low-zmin)/deltaz)
z_norm_high_ind = int((z_norm_high-zmin)/deltaz)
trans2 = sdb.wrap_array("project(fz2,transmission2)").toarray()
mean_trans = np.mean(trans2[z_norm_low_ind:z_norm_high_ind])

fz_faucher_giguere = np.exp(-0.0018*(1+np.arange(zmin,zmin+deltaz*len(trans2),deltaz))**3.92)
fz_norm = np.mean(fz_faucher_giguere[z_norm_low_ind:z_norm_high_ind])

trans2 = trans2*(fz_norm/mean_trans)

transmission = sdb.wrap_array("project(fz,transmission)").toarray()
z = sdb.wrap_array("project(fz,z)").toarray()

trans2_mod = np.append(trans2,np.zeros(len(transmission)-len(trans2)))
trans2_mod = np.where(trans2_mod != 0, trans2_mod, transmission)

trans2_norm_scidb = sdb.from_array(trans2_mod)

sdb.query("rename(%s,fz2_norm)" % trans2_norm_scidb.name)

print "find mean flux transmission (redimension and aggregagte)"
print time.time()-t0

# Cross join <F(z)> back into the big array
sdb.query("store(redimension(bts_good_data_mod, <lambda:double, flux:float, ivar:float, andmask:int64,"
"name:string, redshift:float, lambdar:double, lambda_cont:double, lambda_index:int64, cont_flux2:double NULL DEFAULT null,"
"z_abs:double, numname:int64, a_qso2:double NULL DEFAULT null> [i=0:*,1048576,0,z_abs_index=0:*,%i,0]), bts_good_data_z_redim)" % mean_flux_chunk_size)

print "find mean flux transmission (redimension)"
print time.time()-t0

sdb.query("store(redimension(cross_join(bts_good_data_z_redim, fz2_norm, bts_good_data_z_redim.z_abs_index, fz2_norm.i0), <lambda:double, flux:float, ivar:float, andmask:int64,"
"name:string, redshift:float, lambdar:double, lambda_cont:double, lambda_index:int64, cont_flux2:double NULL DEFAULT null,"
"z_abs:double, z_abs_index:int64, f0:double, numname:int64, a_qso2:double NULL DEFAULT null> [i=0:*,1048576,0]), bts_test)")

# Require <F(z)> to be 1 red of Ly alpha
sdb.query("store(apply(bts_test,transmission,iif(lambdar > %f, 1, f0)),bts_test2)" % lya)
sdb.query("remove(bts_good_data_mod)")
sdb.query("remove(bts_test)")
sdb.query("rename(bts_test2, bts_good_data_mod)")

print "load mean flux back into array (cross join)"
print time.time()-t0

# Remove temporary arrays
sdb.query("remove(tempcont)")
sdb.query("remove(continuum)")
sdb.query("remove(tempfz)")
sdb.query("remove(fz)")
sdb.query("remove(bts2)")
sdb.query("remove(bts3)")
sdb.query("remove(bts_good_data)")
sdb.query("remove(normrange)")
sdb.query("remove(namenum)")
sdb.query("remove(fancynormrange)")
sdb.query("remove(fancy_a_qso)")
sdb.query("remove(bts_good_data_mod)")
sdb.query("remove(bts_good_data_lambda_redim)")
sdb.query("remove(continuum2)")
sdb.query("remove(continuum2_norm)")
sdb.query("remove(blue_amp_meas)")
sdb.query("remove(continuum2_slope_norm)")
sdb.query("remove(bts_good_data_numname_redim)")
sdb.query("remove(bts_good_data_numname_redim2)")
sdb.query("remove(bts_good_data_cj)")
sdb.query("remove(all_chisq)")
sdb.query("remove(bts_forest_region)")
sdb.query("remove(bts_forest_region2)")
sdb.query("remove(bts_forest_region3)")
sdb.query("remove(fz2)")
sdb.query("remove(fz2_norm)")
sdb.query("remove(bts_good_data_z_redim)")