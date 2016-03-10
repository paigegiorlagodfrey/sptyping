from __future__ import division
from BDNYCdb import BDdb
# db=BDdb.get_db('/Users/paigegiorla/Desktop/PG_DB_2_16_15.db')
db=BDdb.get_db('/Users/paigegiorla/Dropbox/BDNYCdb/BDNYC.db')
from matplotlib import pyplot as plt 
import numpy as np
from BDNYCdb import utilities as u
import astropy.units as q
import scipy.stats as s
import small_functions as m
from matplotlib import gridspec
from itertools import cycle
import astrotools as a	
import pickle	
	
def showme(short_name,extraction,band,chisquare,object,output_sorted):	
	'''The plotting code for the spectral typing reduced chi squared routine. This code makes the figure used in the paper
	output_sorted=[chilist,namelist,specidlist,sptlist,templist]
	'''
	plt.figure(figsize=(8,10))
	gs = gridspec.GridSpec(2, 1, height_ratios=[1,2]) 

# plot the chi square values vs spec type with polynomial fit
	ax1 = plt.subplot(gs[0])
   	Ts=['L9','T0','T1','T2','T3','T4','T5','T6','T7','T8','T9']
	plt.scatter(chisquare[0],chisquare[1])
	chisquare=zip(*sorted(zip(*chisquare)))
	pfit = np.polyfit(chisquare[0],chisquare[1], 3)   # Fit a 3rd order polynomial to (x, y) data
	yfit = np.polyval(pfit, chisquare[0])   # Evaluate the polynomial at x
	plt.plot(chisquare[0], yfit)
	ax1.set_xlim(19.8,29.2)
 	plt.locator_params(axis = 'x', nbins = 11)
 	ax1.set_xticklabels(Ts)
	plt.xlabel('Spectral Type',fontsize="x-large")
	plt.ylabel('Reduced Chi Squared',fontsize="x-large")

# plot the object in question
	ax2 = plt.subplot(gs[1])
	lamb=np.gradient(object[0])
	xer=lamb/2
	plt.errorbar(object[0],object[1],xerr=xer,yerr=object[2],fmt=None,ecolor='k',marker='o',label='{}'.format(short_name))



# plot a spectrum from half and one plus/minus the spectral type of the best fitting template and one earlier and later types for show
	half, one, two = [output_sorted[0][3] + 0.5, output_sorted[0][3] - 0.5], [output_sorted[0][3] + 1, output_sorted[0][3] - 1], [output_sorted[0][3] + 2.0, output_sorted[0][3] - 2.0]
	
	half_plus=output_sorted[next(n for n in range(len(output_sorted)) if output_sorted[n][3]==half[0])]
	designation_hp=db.list("select names from sources where id={}".format(half_plus[1])).fetchone()
	half_minus=output_sorted[next(n for n in range(len(output_sorted)) if output_sorted[n][3]==half[1])]
	designation_hm=db.list("select names from sources where id={}".format(half_minus[1])).fetchone()
	
	one_plus=output_sorted[next(n for n in range(len(output_sorted)) if output_sorted[n][3]==one[0])]
	designation_op=db.list("select names from sources where id={}".format(one_plus[1])).fetchone()
	one_minus=output_sorted[next(n for n in range(len(output_sorted)) if output_sorted[n][3]==one[1])]
	designation_om=db.list("select names from sources where id={}".format(one_minus[1])).fetchone()
	
	two_plus=output_sorted[next(n for n in range(len(output_sorted)) if output_sorted[n][3]==two[0])]
	designation_tp=db.list("select names from sources where id={}".format(two_plus[1])).fetchone()
	two_minus=output_sorted[next(n for n in range(len(output_sorted)) if output_sorted[n][3]==two[1])]
	designation_tm=db.list("select names from sources where id={}".format(two_minus[1])).fetchone()

	spectype=u.specType(two_minus[3])
	plt.plot(two_minus[4][0],two_minus[4][1],linestyle=':',color='#8080FF',linewidth=1.85,label='{}, {}'.format(spectype,designation_tm[0].split(',')[0]))
	
	spectype=u.specType(one_minus[3])
	plt.plot(one_minus[4][0],one_minus[4][1],linestyle='-.',color='#0000FF',linewidth=1.5,label='{}, {}'.format(spectype,designation_om[0].split(',')[0]))
	
	spectype=u.specType(half_minus[3])
	plt.plot(half_minus[4][0],half_minus[4][1],linestyle='--',color='#000080',linewidth=1.5,label='{}, {}'.format(spectype,designation_hm[0].split(',')[0]))

# plot the best fitting spectrum template 	
	designation=db.list("select names from sources where id={}".format(output_sorted[0][1])).fetchone()
	spectype=u.specType(output_sorted[0][3])	
	plt.plot(output_sorted[0][4][0],output_sorted[0][4][1],color='r', linewidth=2,label='{}, {}'.format(spectype,designation[0].split(',')[0]))

	spectype=u.specType(half_plus[3])
	plt.plot(half_plus[4][0],half_plus[4][1],linestyle='--',color='#004C00',linewidth=1.5,label='{}, {}'.format(spectype,designation_hp[0].split(',')[0]))
	
	spectype=u.specType(one_plus[3])
	plt.plot(one_plus[4][0],one_plus[4][1],linestyle='-.',color='#009900',linewidth=1.5,label='{}, {}'.format(spectype,designation_op[0].split(',')[0]))
	
	spectype=u.specType(two_plus[3])
	plt.plot(two_plus[4][0],two_plus[4][1],linestyle=':',color='#80CC80',linewidth=1.85,label='{}, {}'.format(spectype,designation_tp[0].split(',')[0]))
	
	plt.legend( loc='upper right',fontsize=2,prop={'size':9.8}, ncol=1, numpoints=1)
	plt.xlabel(r'Wavelength ($\mu$m)', fontsize="x-large")
	plt.ylabel(r' Normalized Flux Density', fontsize="x-large")
	plt.tick_params(labelsize="large")
 	plt.ylim(0.001,1.1)
 	plt.xlim(min(object[0])-0.1,max(object[0])+0.1)

	plt.savefig('/Users/paigegiorla/Code/Python_modules/Publications/'+'{}'.format(short_name)+'/Images/{}_'.format(short_name)+'{}'.format(extraction)+'{}'.format(band)+'.eps')