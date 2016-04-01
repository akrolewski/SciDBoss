import matplotlib.pyplot as plt
import matplotlib
from matplotlib import font_manager

# This file was created by Alex Krolewski on November 12, 2015
# It makes all matplotlib plots prettier!

def prettyplot():
	ticks_font = font_manager.FontProperties(family='Helvetica', style='normal',
		size=16, weight='normal', stretch='normal')

	font = {'family': 'Helvetica', 'size': 10}
	matplotlib.rc('font',**font)
	#matplotlib.rc('ylabel',fontweight='bold',fontsize=18,labelpad=20)
	matplotlib.rcParams['axes.labelsize'] = 18
	matplotlib.rcParams['axes.labelweight'] = 'bold'
	matplotlib.rcParams['axes.titlesize'] = 20
	#matplotlib.rcParams['axes.titleweight'] = 'bold'

	plt.figure()
	ax = plt.axes()
	
	for label in ax.get_xticklabels():
		#print label.get_text()
		label.set_fontproperties(ticks_font)
	for label in ax.get_yticklabels():
		label.set_fontproperties(ticks_font)
	
	plt.minorticks_on()
	plt.tick_params(axis='both', which='major', labelsize=12)
	plt.gcf().subplots_adjust(bottom=0.15)
	plt.gcf().subplots_adjust(left=0.15)
	
	t = plt.title('')
	t.set_y(1.05)
	t.set_fontweight('bold')
	
	x = ax.set_xlabel('',labelpad=20)
	y = ax.set_ylabel('',labelpad=20)