from __future__ import division
from astrodbkit import astrodb
db=astrodb.Database('/Users/paigegiorla/Dropbox/BDNYCdb/BDNYCdev.db')
Jdb=astrodb.Database('/Users/paigegiorla/Code/BDNYCv1.0.db')
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
from astropy.table import Table
	
def showme(short_name,extraction,band,chisquare,object,output_table,spt_range,constraint_on_pfit,order,best_spts,avg_arr=[],plot_polyfit_avg=True,plot_polyfit=False):	
	'''The plotting code for the spectral typing reduced chi squared routine. This code makes the figure used in the paper
	output_table=[chilist,namelist,specidlist,sptlist,templist] astropy table
	If you want to plot the averages of the chisquare values, then run chisquare through m.average_chisq() first and provide here
	Order can be integer or list of two different orders. If list, constraint_on_pfit needs to be two lists also. 
	Half, one, two are the spectral types above/below the object's projected type that you choose to plot against it for visual comparison 
	'''
	plt.figure(figsize=(8,10))
	gs = gridspec.GridSpec(2, 1, height_ratios=[1,2]) 

# plot the chi square values vs spec type with polynomial fit
	ax1 = plt.subplot(gs[0])
	spt_ticks = m.xticks(spt_range)
	
	chisquare = zip(*sorted(zip(*chisquare)))
	polyfit_dict = {}
	output_table.sort('spectral_types')

	for i in range(len(chisquare[0])):
		if output_table[i][5] == 'VL-G' or output_table[i][5] == 'b' or output_table[i][5] == 'g':
			plt.scatter(chisquare[0][i],chisquare[1][i], color='orange')
		else:	
			plt.scatter(chisquare[0][i],chisquare[1][i], color='gray')
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
	
	xmin = spt_range[0]
	xmax = spt_range[1]
	ax1.set_xlim(xmin-0.2,xmax+0.2)
	c=[]
	dup = chisquare
	dup_sort = zip(*sorted(zip(*dup)))
	for i in range(len(dup_sort[0])):
		if dup_sort[0][i] >=xmin and dup_sort[0][i]<=xmax:
			c.append(dup_sort[1][i])
	c = sorted(c)	
		
	ax1.set_ylim(0, c[-1]+1)
 	ax1.set_xticks(np.arange(spt_range[0],spt_range[1]+0.5,1))
 	ax1.set_xticklabels(spt_ticks)
	plt.xlabel('Spectral Type',fontsize="x-large")
	plt.ylabel('Chi Squared',fontsize="x-large")

# plot the object in question
	ax2 = plt.subplot(gs[1])
	lamb=np.gradient(object[0])
	xer=lamb/2
	plt.errorbar(object[0],object[1],xerr=xer,yerr=object[2],fmt=None,ecolor='k',marker='o',label='{}'.format(short_name))



# plot a spectrum from half and one plus/minus the spectral type of the best fitting template and one earlier and later types for show
	
	colors = ['#8080FF','#0000FF','#000080','r','#004C00','#009900','#80CC80']
	linewidths = [1.85,1.5,1.5,1.5,1.5,1.5,1.85]
	linestyles = [':','-.','--','-','--','-.',':']
	designations = []
	for i,spt in zip(range(len(best_spts)),best_spts):
		print i,spt
		obj = output_table[next(n for n in range(len(output_table)) if output_table[n][3]==spt)]
		print obj[6]
		designation = db.query("select names from sources where id={}".format(obj[1]))
		spectype = u.specType(obj[3])
		plt.plot(obj[4][0],obj[4][1],linestyle=linestyles[i],color=colors[i],linewidth=linewidths[i],label='{}, {}'.format(spectype,obj[6]))
		designations.append(designation[0][0].split(',')[0])
		
# 	half_plus = output_table[next(n for n in range(len(output_table)) if output_table[n][3]==next_best[0])]
# 	designation_hp = db.query("select names from sources where id={}".format(half_plus[1]))
# 	half_minus = output_table[next(n for n in range(len(output_table)) if output_table[n][3]==next_best[1])]
# 	designation_hm = db.query("select names from sources where id={}".format(half_minus[1]))
# 	
# 	one_plus = output_table[next(n for n in range(len(output_table)) if output_table[n][3]==next_best[2])]
# 	designation_op = db.query("select names from sources where id={}".format(one_plus[1]))
# 	one_minus = output_table[next(n for n in range(len(output_table)) if output_table[n][3]==next_best[3])]
# 	designation_om = db.query("select names from sources where id={}".format(one_minus[1]))
# 	
# 	two_plus = output_table[next(n for n in range(len(output_table)) if output_table[n][3]==next_best[4])]
# 	designation_tp = db.query("select names from sources where id={}".format(two_plus[1]))
# 	two_minus = output_table[next(n for n in range(len(output_table)) if output_table[n][3]==next_best[5])]
# 	designation_tm = db.query("select names from sources where id={}".format(two_minus[1]))
# 	spectype = u.specType(two_minus[3])
# 	plt.plot(two_minus[4][0],two_minus[4][1],linestyle=':',color='#8080FF',linewidth=1.85,label='{}, {}'.format(spectype,two_minus[6]))
# 	
# 	spectype = u.specType(one_minus[3])
# 	plt.plot(one_minus[4][0],one_minus[4][1],linestyle='-.',color='#0000FF',linewidth=1.5,label='{}, {}'.format(spectype,one_minus[6]))
# 	
# 	spectype = u.specType(half_minus[3])
# 	plt.plot(half_minus[4][0],half_minus[4][1],linestyle='--',color='#000080',linewidth=1.5,label='{}, {}'.format(spectype,half_minus[6]))
# 
# # plot the best fitting spectrum template 	
# 	best = output_table[next(n for n in range(len(output_table)) if output_table[n][3]==best_spt)]
# 	designation = db.query("select names from sources where id={}".format(best[1]))
# 	spectype = u.specType(best[3])	
# 	plt.plot(best[4][0],best[4][1],color='r', linewidth=2,label='{}, {}'.format(spectype,best[6]))
# 
# 	spectype = u.specType(half_plus[3])
# 	plt.plot(half_plus[4][0],half_plus[4][1],linestyle='--',color='#004C00',linewidth=1.5,label='{}, {}'.format(spectype,half_plus[6]))
# 	
# 	spectype = u.specType(one_plus[3])
# 	plt.plot(one_plus[4][0],one_plus[4][1],linestyle='-.',color='#009900',linewidth=1.5,label='{}, {}'.format(spectype,one_plus[6]))
# 	
# 	spectype = u.specType(two_plus[3])
# 	plt.plot(two_plus[4][0],two_plus[4][1],linestyle=':',color='#80CC80',linewidth=1.85,label='{}, {}'.format(spectype,two_plus[6]))
	
	plt.legend( loc='best',fontsize=0.5,prop={'size':7.8}, ncol=1, numpoints=1)
	plt.xlabel(r'Wavelength ($\mu$m)', fontsize="x-large")
	plt.ylabel(r' Normalized Flux Density', fontsize="x-large")
	plt.tick_params(labelsize="large")
 	plt.ylim(0.001,max(object[1])+0.5)
 	plt.xlim(min(object[0])-0.01,max(object[0])+0.01)

	pickle.dump(polyfit_dict,open('/Users/paigegiorla/Publications/'+'{}'.format(short_name)+'/Results/polynomial_{}'.format(extraction)+'.pkl','wb'))

	plt.savefig('/Users/paigegiorla/Publications/'+'{}'.format(short_name)+'/Images/{}_'.format(short_name)+'{}'.format(extraction)+'{}'.format(band)+'.eps')
	plt.clf()
	return polyfit_dict, designations