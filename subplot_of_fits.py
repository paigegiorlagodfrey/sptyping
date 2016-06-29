import pickle
from matplotlib import pyplot as plt
import astropy.units as q
from pylab import *
from matplotlib import colors
import matplotlib as mpl
import numpy as np
from astrodbkit import astrodb
db=astrodb.Database('/Users/paigegiorla/Dropbox/BDNYCdb/BDNYCdev.db')
from BDNYCdb import utilities as u
import astropy.table as table

def scrub(data):
  '''
  For input data [w,f,e] or [w,f] returns the list with NaN, negative, and zero flux (and corresponsing wavelengths and errors) removed. 
  '''
  data = [i*q.Unit('') for i in data]
  data = [i[np.where(np.logical_and(data[1].value>0,~np.isnan(data[1].value)))] for i in data]
  data = [i[np.unique(data[0], return_index=True)[1]] for i in data]
  return [i[np.lexsort([data[0]])] for i in data]

def showme(best_spts=False,plot_these=['1256-1257b','2114-2251','1741-4642'],individual_objs=['1256-1257b', '0437-0228b', '2307+2108b', '2114-2251', '1300+1221', '0112+1703','2144+1446','1741-4642']):
	if individual_objs:
		io_colors = ['b', 'c', 'm', 'purple', 'orange', 'fuchsia','lime','c']
		io_names = ['VHS 1256B', '51 Eri b', 'HR 8799b', 'PSO 318', 'Ross 458C', 'GU Psc b','HN Peg b','1741-4642']

	files = [('J','/Users/paigegiorla/Publications/HD2562B/Results/J_062816.pkl','/Users/paigegiorla/Code/Data/GOI28/J_Spectrum_GOI28_JHK1K2_v5.ascii'),('H','/Users/paigegiorla/Publications/HD2562B/Results/H_062816.pkl','/Users/paigegiorla/Code/Data/GOI28/H_Spectrum_GOI28_JHK1K2_v5.ascii'),('K','/Users/paigegiorla/Publications/HD2562B/Results/K_062816.pkl','/Users/paigegiorla/Code/Data/GOI28/K1K2_Spectrum_GOI28_JHK1K2_v5.ascii'),('full','/Users/paigegiorla/Publications/HD2562B/Results/full_notfull_062816.pkl','/Users/paigegiorla/Code/Data/GOI28/edited_Spectrum_GOI28_JHK1K2_v5.ascii')]
	lims = [[1.143, 1.375],[1.375,2.000],[1.937,2.4],[0.9,2.4]]

	fontsize = ['large','large','large','x-large']
	axJ = plt.subplot2grid((2,3), (0,0))
	axH = plt.subplot2grid((2,3), (0,1))
	axK = plt.subplot2grid((2,3), (0,2))
	axf = plt.subplot2grid((2,3), (1,0), colspan=3)
	axJ.set_xticklabels([])
	axH.set_xticklabels([])
	axK.set_xticklabels([])
	subplots = [axJ, axH, axK, axf]
	bestspts_colors = ['#8080FF','#0000FF','#000080','r','green','#009900','#80CC80']
	linewidths = [1.85,1.5,1.5,1.5,1.5,1.5,1.85]
	linestyles = [':','-.','--','-','--','-.',':']
	colors = ['orchid','orange','darkorange','palevioletred']
	fontsize = ['large','large','large','large','x-large']
	leg1 = Rectangle((0, 0), 0, 0, alpha=0.0)

	for k in range(len(files)):
		object_file = open(files[k][2], 'rb').readlines()
		W_obj,F_obj,U_obj = [],[],[]
		for line in object_file:
			columns=line.split()
			wav=float(columns[0])
			flu=float(columns[1])
			un=float(columns[2])
			W_obj.append(wav) ; F_obj.append(flu) ; U_obj.append(un)

		object = [W_obj,F_obj,U_obj]
		object = scrub(object)	

		fb = open(files[k][1],'rb')
		output = pickle.load(fb)
		fb.close()
		output_table = table.Table(output,names=('chi_values','source_id','spectra_id','spectral_types','template','gravity','shortnames'), meta={'name':'results from GOF'})

		output_table.sort('chi_values')

		subplots[k].errorbar(object[0],object[1],yerr=object[2],fmt=None,ecolor='gray',marker='o',label='HD2562B')
		designations = []

		if best_spts:
			for j,spt in zip(range(len(best_spts)),best_spts):
				obj = output_table[next(n for n in range(len(output_table)) if output_table[n][3]==spt)]
				if best_spts[j]==best_spts[j-1]:
					obj = output_table[next(n+1 for n in range(len(output_table)) if output_table[n][3]==spt)]
				designation = db.query("select names from sources where id={}".format(obj[1]))
				spectype = u.specType(obj[3])
				grav = ('' if obj[5]==None else obj[5])
				name = (io_names[individual_objs.index(obj[6])] if obj[6] in individual_objs else obj[6])
				subplots[k].plot(obj[4][0],obj[4][1],linestyle=linestyles[j],color=bestspts_colors[j],linewidth=linewidths[j],label='{}{}, {}'.format(spectype,grav,name))
				designations.append(designation[0][0].split(',')[0])
		else:
			obj = output_table[0]
			designation = db.query("select names from sources where id={}".format(obj[1]))
			spectype = u.specType(obj[3])
			grav = ('' if obj[5]==None else obj[5])
			name = (io_names[individual_objs.index(obj[6])] if obj[6] in individual_objs else obj[6])
			subplots[k].plot(obj[4][0],obj[4][1],linestyle='-',color=colors[k],linewidth=2,label='{}{}, {}'.format(spectype,grav[:4],name))
				
			if plot_these:
				pcolors = ['blue','red','mediumturquoise']
				for p in range(len(plot_these)):	
					obj_p = output_table[[n for n in range(len(output_table)) if output_table[n][6]==plot_these[p]]][0]
					if obj[6]==obj_p[6]:
						print 'best fit already plotted'
					else:
						spectype = u.specType(obj_p[3])
						grav = ('' if obj_p[5]==None else obj_p[5])
						name = (io_names[individual_objs.index(obj_p[6])] if obj_p[6] in individual_objs else obj_p[6])
						subplots[k].plot(obj_p[4][0],obj_p[4][1],linestyle='--',color=pcolors[p],linewidth=1.5,label='{}{}, {}'.format(spectype,grav[:4],name))

				designations.append(designation[0][0].split(',')[0])
				legend = subplots[k].legend(loc='best',frameon=False,fontsize = 'xx-small',handlelength=0)

		subplots[k].set_xlim(min(object[0])-0.01,max(object[0])+0.01)
		min_lim, max_lim = [],[]
		for i in range(len(object[1])):
			min_lim.append(object[1][i]-object[2][i])
			max_lim.append(object[1][i]+object[2][i])
		y_lims = (min(min_lim),max(max_lim))	
		print y_lims
 		subplots[k].set_ylim(y_lims[0],y_lims[1])
 		subplots[k].set_yticklabels([])


	plt.subplots_adjust(wspace=0,hspace=0)	
	plt.xlabel('Wavelength ($\mu$m)',fontsize="xx-large")
	plt.ylabel('Flux',fontsize="xx-large")
	plt.savefig('/Users/paigegiorla/Publications/HD2562B/Images/subplot.eps')
	plt.clf()	
	return designations