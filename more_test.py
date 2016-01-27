a_qso = all_aqso
chisq_dof = all_chisq_dof

fz_faucher_giguere = np.exp(-0.0018*(1+fz[:,0])**3.92)
fz_norm = np.mean(fz_faucher_giguere[(fz[:,0] > znorm_low) & (fz[:,0] < znorm_high)])
# from eqn. 21 in Faucher-Giguere et al.
fz_val = np.mean(fz[(fz[:,0] > znorm_low ) & (fz[:,0] < znorm_high),1])
fz[:,1] = fz[:,1]*(fz_norm/fz_val)

stat = np.nansum(((fz[:,1]-orig_fz[:,1])/fz[:,2])**2)

# plots
prettyplot()
plt.plot(continuum[:,0],continuum[:,1])
plt.savefig("continuum_%i.pdf" % iteration)

prettyplot()
plt.plot(fz[:,0],fz[:,1])
plt.plot(fz[:,0],fz[:,1])
plt.savefig("fz_%i.pdf" % iteration)

# write to text file
mean_spec_file = open('continuum%i.dat' % iteration,'w')
for i in range(len(continuum)):
	mean_spec_file.write('%f %f\n' % (continuum[i,0],continuum[i,1]))
mean_spec_file.close()
	
fz_file = open('fz%i.dat' % iteration,'w')
for i in range(len(fz)):
	fz_file.write('%f %f\n' % (fz[i,0],fz[i,1]))
fz_file.close()
