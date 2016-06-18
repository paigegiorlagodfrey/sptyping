from __future__ import division
from astrodbkit import astrodb
db=astrodb.Database('/Users/paigegiorla/Dropbox/BDNYCdb/BDNYCdev.db')
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
	
def showme(short_name,extraction,band,chisquare,object,output_table,spt_range,constraint_on_pfit,order,avg_arr=[],plot_polyfit_avg=True,plot_polyfit=False):	
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
		if output_table[i][3] == 'f':
			print 'f:',output_table[i][3], output_table[i][1],chisquare[0][i]
			plt.scatter(chisquare[0][i],chisquare[1][i], color='gray')
		else:	
			print 'else:',output_table[i][3], output_table[i][1],chisquare[0][i]
			plt.scatter(chisquare[0][i],chisquare[1][i], color='orange')
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
	ax1.set_ylim(min(c)-1, max(c)+1)
 	ax1.set_xticks(np.arange(spt_range[0],spt_range[1]+0.5,1))
 	ax1.set_xticklabels(spt_ticks)
	plt.xlabel('Spectral Type',fontsize="x-large")
	plt.ylabel('Chi Squared',fontsize="x-large")
	plt.tight_layout()

# plot the object in question
	ax2 = plt.subplot(gs[1])
	lamb=np.gradient(object[0])
	xer=lamb/2
	plt.errorbar(object[0],object[1],xerr=xer,yerr=object[2],fmt=None,ecolor='k',marker='o',label='{}'.format(short_name))


# plot L5f,L6f,L7f,L8f

	five = output_table[next(n for n in range(len(output_table)) if (output_table[n][1]==15.0 and output_table[n][3]=='g'))]
	spectype = a.specType(five[1])
	plt.plot(five[2][0],five[2][1],color='r',linewidth=1.85,label='{}, {}'.format(spectype,five[3]))
	five_2 = output_table[next(n for n in range(len(output_table)) if (output_table[n][1]==15.0 and output_table[n][3]=='f'))]
	spectype = a.specType(five_2[1])
	plt.plot(five_2[2][0],five_2[2][1],color='m',linewidth=1.85,label='{}, {}'.format(spectype,five_2[3]))
	six = output_table[next(n for n in range(len(output_table)) if output_table[n][1]==16.0)]
	spectype = a.specType(six[1])
	plt.plot(six[2][0],six[2][1],linestyle='--',color='#004C00',linewidth=1.5,label='{}, {}'.format(spectype,six[3]))
	seven = output_table[next(n for n in range(len(output_table)) if output_table[n][1]==17.0)]
	spectype = a.specType(seven[1])
	plt.plot(seven[2][0],seven[2][1],linestyle='-.',color='#009900',linewidth=1.5,label='{}, {}'.format(spectype,seven[3]))
	eight = output_table[next(n for n in range(len(output_table)) if output_table[n][1]==18.0)]
	spectype = a.specType(eight[1])
	plt.plot(eight[2][0],eight[2][1],linestyle=':',color='#8080FF',linewidth=1.85,label='{}, {}'.format(spectype,eight[3]))
 	
	plt.legend( loc='best',fontsize=2,prop={'size':9.8}, ncol=1, numpoints=1)
	plt.xlabel(r'Wavelength ($\mu$m)', fontsize="x-large")
	plt.ylabel(r' Normalized Flux Density', fontsize="x-large")
	plt.tick_params(labelsize="large")
 	plt.ylim(0.001,1.1)
 	plt.xlim(min(object[0])-0.01,max(object[0])+0.01)

	pickle.dump(polyfit_dict,open('/Users/paigegiorla/Publications/'+'{}'.format(short_name)+'/Results/polynomial_{}'.format(extraction)+'.pkl','wb'))

	plt.savefig('/Users/paigegiorla/Publications/'+'{}'.format(short_name)+'/Images/{}_'.format(short_name)+'{}'.format(extraction)+'{}'.format(band)+'.eps')
	plt.clf()
	return polyfit_dict