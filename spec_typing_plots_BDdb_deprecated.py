from __future__ import division
from BDNYCdb import BDdb
db=BDdb.get_db('/Users/paigegiorla/Dropbox/BDNYCdb/BDNYCdeprecated.db')
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
	
def showme(short_name,extraction,band,chisquare,object,output_table,spt_range,constraint_on_pfit,order,best_spt,next_best,avg_arr=[],plot_polyfit_avg=True,plot_polyfit=False):	
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

	chisquare = zip(*sorted(zip(*chisquare)))
	polyfit_dict = {}


	for i in range(len(chisquare[0])):
		if output_table[i][5] == None or output_table[i][5] == '':
			plt.scatter(chisquare[0][i],chisquare[1][i], color='gray')
		else:	
			plt.scatter(chisquare[0][i],chisquare[1][i], color='g')
	if plot_polyfit==True:	
		if isinstance(order,list):
			pfit_1,yfit_1,chisquare_1 = m.polynomialfit(chisquare,constraint_on_pfit[0],order[0])
			pfit_2,yfit_2,chisquare_2 = m.polynomialfit(chisquare,constraint_on_pfit[1],order[1])
			plt.plot(chisquare_1[0], yfit_1,'b',linewidth=1.85)
			plt.plot(chisquare_2[0], yfit_2,'r',linewidth=1.85)
			polyfit_dict['pfit_{}'.format(order[0])] = pfit_1
			polyfit_dict['pfit_{}'.format(order[1])] = pfit_2
			polyfit_dict['yfit_{}'.format(order[0])] = yfit_1
			polyfit_dict['yfit_{}'.format(order[1])] = yfit_2
		else:
			pfit,yfit,chisquare = m.polynomialfit(chisquare,constraint_on_pfit,order)
			plt.plot(chisquare[0], yfit,'b',linewidth=1.85)
			polyfit_dict['pfit_{}'.format(order)] = pfit
			polyfit_dict['yfit_{}'.format(order)] = yfit

	if avg_arr:
		avg_arr = zip(*sorted(zip(*avg_arr)))
		plt.scatter(avg_arr[0],avg_arr[1],color='k')
		if plot_polyfit_avg==True:
			if isinstance(order,list):
				pfit_1,yfit_1,avg_arr_1 = m.polynomialfit(avg_arr,constraint_on_pfit[0],order[0])
				pfit_2,yfit_2,avg_arr_2 = m.polynomialfit(avg_arr,constraint_on_pfit[1],order[1])
				plt.plot(avg_arr_1[0], yfit_1,'c',linestyle='--',linewidth=1.85)
				plt.plot(avg_arr_2[0], yfit_2,'m',linestyle='--',linewidth=1.85)
				polyfit_dict['pfit_avg_{}'.format(order[0])] = pfit_1
				polyfit_dict['pfit_avg_{}'.format(order[1])] = pfit_2
				polyfit_dict['yfit_avg_{}'.format(order[0])] = yfit_1
				polyfit_dict['yfit_avg_{}'.format(order[1])] = yfit_2
			else:
				pfit,yfit,avg_arr = m.polynomialfit(avg_arr,constraint_on_pfit,order)
				plt.plot(avg_arr[0], yfit,'c',linestyle='--',linewidth=1.85)
				polyfit_dict['pfit_avg_{}'.format(order)] = pfit
				polyfit_dict['yfit_avg_{}'.format(order)] = yfit
	
	ax1.set_xlim(spt_range[0]-0.2,spt_range[1]+0.2)
 	ax1.set_xticks(np.arange(spt_range[0],spt_range[1]+0.5,1))
 	ax1.set_xticklabels(spt_ticks)
	plt.xlabel('Spectral Type',fontsize="x-large")
	plt.ylabel('Reduced Chi Squared',fontsize="x-large")

# plot the object in question
	ax2 = plt.subplot(gs[1])
	lamb=np.gradient(object[0])
	xer=lamb[0]/2
	plt.errorbar(object[0],object[1],xerr=xer,yerr=object[2],fmt=None,ecolor='k',marker='o',label='{}'.format(short_name))



# plot a spectrum from half and one plus/minus the spectral type of the best fitting template and one earlier and later types for show
	
	half_plus = output_table[next(n for n in range(len(output_table)) if output_table[n][3]==next_best[0])]
	designation_hp = db.list("select names from sources where id={}".format(half_plus[1])).fetchone()
	half_minus = output_table[next(n for n in range(len(output_table)) if output_table[n][3]==next_best[1])]
	designation_hm = db.list("select names from sources where id={}".format(half_minus[1])).fetchone()
	
	one_plus = output_table[next(n for n in range(len(output_table)) if output_table[n][3]==next_best[2])]
	designation_op = db.list("select names from sources where id={}".format(one_plus[1])).fetchone()
	one_minus = output_table[next(n for n in range(len(output_table)) if output_table[n][3]==next_best[3])]
	designation_om = db.list("select names from sources where id={}".format(one_minus[1])).fetchone()
	
	two_plus = output_table[next(n for n in range(len(output_table)) if output_table[n][3]==next_best[4])]
	designation_tp = db.list("select names from sources where id={}".format(two_plus[1])).fetchone()
	two_minus = output_table[next(n for n in range(len(output_table)) if output_table[n][3]==next_best[5])]
	designation_tm = db.list("select names from sources where id={}".format(two_minus[1])).fetchone()	
	
	spectype = u.specType(half_plus[3])
	plt.plot(half_plus[4][0],half_plus[4][1],linestyle=':',color='#80CC80',linewidth=1.5,label='{}, {}'.format(spectype,half_plus[6]))

	spectype = u.specType(half_minus[3])
	plt.plot(half_minus[4][0],half_minus[4][1],linestyle='-.',color='#009900',linewidth=1.5,label='{}, {}'.format(spectype,half_minus[6]))
	
	spectype = u.specType(one_plus[3])
	plt.plot(one_plus[4][0],one_plus[4][1],linestyle='--',color='#004C00',linewidth=1.5,label='{}, {}'.format(spectype,one_plus[6]))

	spectype = u.specType(one_minus[3])
	plt.plot(one_minus[4][0],one_minus[4][1],linestyle='--',color='#000080',linewidth=1.5,label='{}, {}'.format(spectype,one_minus[6]))

	spectype = u.specType(two_plus[3])
	plt.plot(two_plus[4][0],two_plus[4][1],linestyle='-.',color='#0000FF',linewidth=1.85,label='{}, {}'.format(spectype,two_plus[6]))

	spectype = u.specType(two_minus[3])
	plt.plot(two_minus[4][0],two_minus[4][1],linestyle=':',color='#8080FF',linewidth=1.85,label='{}, {}'.format(spectype,two_minus[6]))


	# plot the best fitting spectrum template 	
	best = output_table[next(n for n in range(len(output_table)) if output_table[n][3]==best_spt)]
	designation = db.list("select names from sources where id={}".format(best[1])).fetchone()
	spectype = u.specType(best[3])	
	plt.plot(best[4][0],best[4][1],color='r', linewidth=2,label='{}, {}'.format(spectype,best[6]))

	plt.legend( loc='best',fontsize=2,prop={'size':9.8}, ncol=1, numpoints=1)
	plt.xlabel(r'Wavelength ($\mu$m)', fontsize="x-large")
	plt.ylabel(r' Normalized Flux Density', fontsize="x-large")
	plt.tick_params(labelsize="large")
 	plt.ylim(0.001,1.1)
 	plt.xlim(min(object[0])-0.01,max(object[0])+0.01)
 	
	pickle.dump(polyfit_dict,open('/Users/paigegiorla/Publications/'+'{}'.format(short_name)+'/Results/polynomial_{}'.format(extraction)+'.pkl','wb'))

	plt.savefig('/Users/paigegiorla/Publications/'+'{}'.format(short_name)+'/Images/{}_'.format(short_name)+'{}'.format(extraction)+'{}'.format(band)+'.eps')
	designations = [designation[0].split(',')[0],designation_hp[0].split(',')[0],designation_hm[0].split(',')[0],designation_op[0].split(',')[0],designation_om[0].split(',')[0],designation_tp[0].split(',')[0],designation_tm[0].split(',')[0]]
	plt.clf()
	return polyfit_dict, designations