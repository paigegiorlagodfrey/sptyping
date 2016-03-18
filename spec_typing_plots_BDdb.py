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
	
def showme(short_name,extraction,band,chisquare,object,output_table,spt_range,constraint_on_pfit,order,avg_arr=[]):	
	'''The plotting code for the spectral typing reduced chi squared routine. This code makes the figure used in the paper
	output_table=[chilist,namelist,specidlist,sptlist,templist] astropy table
	If you want to plot the averages of the chisquare values, then run chisquare through m.average_chisq() first and provide here
	Order can be integer or list of two different orders. If list, constraint_on_pfit needs to be two lists also. 
	'''
	plt.figure(figsize=(8,10))
	gs = gridspec.GridSpec(2, 1, height_ratios=[1,2]) 

# plot the chi square values vs spec type with polynomial fit
	spt_ticks = []
	ax1 = plt.subplot(gs[0])
	for spt in np.arange(spt_range[0],spt_range[1]+0.5,1):
		spt_tick = a.specType(spt)
		spt_ticks.append(spt_tick)

	if avg_arr:
		plt.scatter(chisquare[0],chisquare[1], color='gray')
		for i in range(len(chisquare[0])):
			if output_table[i][5] == 'VL-G':
				plt.scatter(chisquare[0][i]+np.random.random()/10,chisquare[1][i], color='g')
		plt.scatter(avg_arr[0],avg_arr[1],color='k')
	else:
		for i in range(len(chisquare[0])):
			if output_table[i][5] == None or output_table[i][5] == '':
				plt.scatter(chisquare[0][i]+np.random.random()/10,chisquare[1][i], color='k')
			else:	
				plt.scatter(chisquare[0][i]+np.random.random()/10,chisquare[1][i], color='g')

	chisquare=zip(*sorted(zip(*chisquare)))
	
	if isinstance(order,list):
		pfit_1,yfit_1,chisquare_1 = m.polynomialfit(chisquare,constraint_on_pfit[0],order[0])
		pfit_2,yfit_2,chisquare_2 = m.polynomialfit(chisquare,constraint_on_pfit[1],order[1])
 		plt.plot(chisquare_1[0], yfit_1,'b')
 		plt.plot(chisquare_2[0], yfit_2,'r')
	else:
		pfit,yfit,chisquare = m.polynomialfit(chisquare,constraint_on_pfit,order)
 		plt.plot(chisquare[0], yfit,'b')
	
	ax1.set_xlim(spt_range[0]-0.2,spt_range[1]+0.2)
 	ax1.set_xticks(np.arange(spt_range[0],spt_range[1]+0.5,1))
 	ax1.set_xticklabels(spt_ticks)
	plt.xlabel('Spectral Type',fontsize="x-large")
	plt.ylabel('Reduced Chi Squared',fontsize="x-large")

# plot the object in question
	ax2 = plt.subplot(gs[1])
	lamb=np.gradient(object[0])
	xer=lamb/2
	plt.errorbar(object[0],object[1],xerr=xer,yerr=object[2],fmt=None,ecolor='k',marker='o',label='{}'.format(short_name))



# plot a spectrum from half and one plus/minus the spectral type of the best fitting template and one earlier and later types for show
	half, one, two = [output_table[0][3] + 0.5, output_table[0][3] - 0.5], [output_table[0][3] + 1, output_table[0][3] - 1], [output_table[0][3] + 2.0, output_table[0][3] - 2.0]
	
	half_plus=output_table[next(n for n in range(len(output_table)) if output_table[n][3]==half[0])]
	designation_hp=db.list("select names from sources where id={}".format(half_plus[1])).fetchone()
	half_minus=output_table[next(n for n in range(len(output_table)) if output_table[n][3]==half[1])]
	designation_hm=db.list("select names from sources where id={}".format(half_minus[1])).fetchone()
	
	one_plus=output_table[next(n for n in range(len(output_table)) if output_table[n][3]==one[0])]
	designation_op=db.list("select names from sources where id={}".format(one_plus[1])).fetchone()
	one_minus=output_table[next(n for n in range(len(output_table)) if output_table[n][3]==one[1])]
	designation_om=db.list("select names from sources where id={}".format(one_minus[1])).fetchone()
	
	two_plus=output_table[next(n for n in range(len(output_table)) if output_table[n][3]==two[0])]
	designation_tp=db.list("select names from sources where id={}".format(two_plus[1])).fetchone()
	two_minus=output_table[next(n for n in range(len(output_table)) if output_table[n][3]==two[1])]
	designation_tm=db.list("select names from sources where id={}".format(two_minus[1])).fetchone()

	spectype=u.specType(two_minus[3])
	plt.plot(two_minus[4][0],two_minus[4][1],linestyle=':',color='#8080FF',linewidth=1.85,label='{}, {}'.format(spectype,designation_tm[0].split(',')[0]))
	
	spectype=u.specType(one_minus[3])
	plt.plot(one_minus[4][0],one_minus[4][1],linestyle='-.',color='#0000FF',linewidth=1.5,label='{}, {}'.format(spectype,designation_om[0].split(',')[0]))
	
	spectype=u.specType(half_minus[3])
	plt.plot(half_minus[4][0],half_minus[4][1],linestyle='--',color='#000080',linewidth=1.5,label='{}, {}'.format(spectype,designation_hm[0].split(',')[0]))

# plot the best fitting spectrum template 	
	designation=db.list("select names from sources where id={}".format(output_table[0][1])).fetchone()
	spectype=u.specType(output_table[0][3])	
	plt.plot(output_table[0][4][0],output_table[0][4][1],color='r', linewidth=2,label='{}, {}'.format(spectype,designation[0].split(',')[0]))

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
 	plt.xlim(min(object[0])-0.01,max(object[0])+0.01)

	plt.savefig('/Users/paigegiorla/Publications/'+'{}'.format(short_name)+'/Images/{}_'.format(short_name)+'{}'.format(extraction)+'{}'.format(band)+'.eps')
	
	if isinstance(order,int):
	 	return pfit,yfit