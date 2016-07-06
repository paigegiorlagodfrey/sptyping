import pickle
from matplotlib import pyplot as plt
import astropy.units as q
from pylab import *
from matplotlib import colors
import matplotlib as mpl
import numpy as np
from astrodbkit import astrodb
import small_functions as m
from matplotlib import gridspec
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

def showme(spt_range=[(10,25),(10,25),(10,25),(0,0)],plot_polyfit=False,best_spts=False,plot_these=['1256-1257b','2114-2251','1741-4642','1516+3053'],individual_objs=['1256-1257b', '2114-2251', '1741-4642','1516+3053']):

	files = [('J','/Users/paigegiorla/Publications/HD2562B/Results/J_070216.pkl','/Users/paigegiorla/Code/Data/GOI28/J_Spectrum_GOI28_JHK1K2_v5.ascii'),('H','/Users/paigegiorla/Publications/HD2562B/Results/H_070216.pkl','/Users/paigegiorla/Code/Data/GOI28/H_Spectrum_GOI28_JHK1K2_v5.ascii'),('K','/Users/paigegiorla/Publications/HD2562B/Results/K_070216.pkl','/Users/paigegiorla/Code/Data/GOI28/K1K2_Spectrum_GOI28_JHK1K2_v5.ascii'),('full','/Users/paigegiorla/Publications/HD2562B/Results/full_notfull_070216.pkl','/Users/paigegiorla/Code/Data/GOI28/edited_Spectrum_GOI28_JHK1K2_v5.ascii')]
	lims = [[1.143, 1.375],[1.375,2.000],[1.937,2.4],[0.9,2.4]]

	fontsize = ['large','large','large','x-large']
	axJ1 = plt.subplot2grid((3,3), (0,0))
	axJ2 = plt.subplot2grid((3,3), (1,0))
	axH1 = plt.subplot2grid((3,3), (0,1))
	axH2 = plt.subplot2grid((3,3), (1,1))
	axK1 = plt.subplot2grid((3,3), (0,2))
	axK2 = plt.subplot2grid((3,3), (1,2))
	axf = plt.subplot2grid((3,3), (2,0), colspan=3)

	axJ2.set_xticklabels([])
	axH2.set_xticklabels([])
	axK2.set_xticklabels([])
	axf.set_yticklabels([])
	axH1.get_yaxis().set_visible(False)
	axK1.get_yaxis().set_visible(False)
	axJ1.tick_params(labelbottom='off',labeltop='on',labelsize=6)
	axH1.tick_params(labelbottom='off',labeltop='on',labelsize=6)
	axK1.tick_params(labelbottom='off',labeltop='on',labelsize=6)
	subplots = [(axJ1, axJ2), (axH1, axH2), (axK1, axK2), (0,axf)]
	bestspts_colors = ['#8080FF','#0000FF','#000080','r','green','#009900','#80CC80']
	linewidths = [1.85,1.5,1.5,1.5,1.5,1.5,1.85]
	linestyles = [':','-.','--','-','--','-.',':']
	colors = ['orchid','orange','darkorange','mediumturquoise']
	fontsize = ['large','large','large','large','x-large']

	designations = []
	for k in range(len(files)):
		fb = open(files[k][1],'rb')
		output = pickle.load(fb)
		fb.close()
		output_table = table.Table(output,names=('chi_values','source_id','spectra_id','spectral_types','template','gravity','shortnames'), meta={'name':'results from GOF'})

		if isinstance(subplots[k][0],int):
			print 'full spectrum time'
		else:	
			output_table.sort('spectral_types')
			chisquare = [output_table['spectral_types'],output_table['chi_values']]
			chisquare = zip(*sorted(zip(*chisquare)))
			spt_ticks = m.xticks([spt_range[k][0],spt_range[k][1]])[1::2]
			polyfit_dict = {}

			for i in range(len(chisquare[0])):
				rand = np.random.random()/8
				if output_table[i][5] == 'VL-G' or output_table[i][5] == 'b' or output_table[i][5] == 'g':
					subplots[k][0].scatter(chisquare[0][i]+rand,chisquare[1][i], color='orange') 		
				else:	
					subplots[k][0].scatter(chisquare[0][i]+rand,chisquare[1][i], color='gray')
				if plot_these:
					io_colors = ['blue','red','mediumturquoise','lime']
					io_names = ['VHS 1256B', 'PSO 318','1741-4642','1516+3053']
					for io in range(len(plot_these)):
						if output_table[i][6] == plot_these[io]:
							subplots[k][0].plot(chisquare[0][i]+rand,chisquare[1][i],linestyle='None',marker='D',markeredgewidth=2,markeredgecolor=io_colors[io], markerfacecolor='None',label=io_names[io])

			if plot_polyfit==True:	
					pfit,yfit,chisquare = m.polynomialfit(chisquare,[10,25],2)
					subplots[k][0].plot(chisquare[0], yfit,'lime',linewidth=1.85)

			xmin = spt_range[k][0]
			xmax = spt_range[k][1]
 			subplots[k][0].set_xlim(xmin-0.2,xmax+0.2)
			c=[]
			dup = chisquare
			dup_sort = zip(*sorted(zip(*dup)))
			for i in range(len(dup_sort[0])):
				if dup_sort[0][i] >=xmin and dup_sort[0][i]<=xmax:
					c.append(dup_sort[1][i])
			c = sorted(c)	
			avg = sum(c)/len(c)	
 	 		subplots[k][0].set_ylim(-10, 400)
			subplots[k][0].set_xticks(np.arange(spt_range[k][0],spt_range[k][1]+0.5,1)[1::2])
			subplots[k][0].set_xticklabels(spt_ticks)
			subplots[k][0].set_xlabel('Spectral Type',fontsize="x-large")
			subplots[k][0].set_ylabel('Chi Squared',fontsize="large")		
		
##-------------------------		
		
		object_file = open(files[k][2], 'rb').readlines()
		W_obj,F_obj,U_obj = [],[],[]
		for line in object_file:
			columns=line.split()
			wav=float(columns[0])
			flu=float(columns[1])
			un=float(columns[2])
			W_obj.append(wav) ; F_obj.append(flu) ; U_obj.append(un)
		
		if k==2 or k==3:
			object = [W_obj[:-1],F_obj[:-1],U_obj[:-1]]
		else:
			object = [W_obj,F_obj,U_obj]
		object = scrub(object)	
		subplots[k][1].scatter(object[0],object[1],color='gray',s=3)
		subplots[k][1].errorbar(object[0],object[1],yerr=object[2],fmt=None,ecolor='gray',marker='o',label='HR 2562B')
		
		output_table.sort('chi_values')
		if best_spts:
			for j,spt in zip(range(len(best_spts)),best_spts):
				obj = output_table[next(n for n in range(len(output_table)) if output_table[n][3]==spt)]
				if best_spts[j]==best_spts[j-1]:
					obj = output_table[next(n+1 for n in range(len(output_table)) if output_table[n][3]==spt)]
				designation = db.query("select names from sources where id={}".format(obj[1]))
				spectype = u.specType(obj[3])
				grav = ('' if (obj[5]==None or obj[5]=='None') else obj[5])
				name = (io_names[individual_objs.index(obj[6])] if obj[6] in individual_objs else obj[6])
				subplots[k][1].plot(obj[4][0],obj[4][1],linestyle=linestyles[j],color=bestspts_colors[j],linewidth=linewidths[j],label='{}{}, {}'.format(spectype,grav,name))
				designations.append(designation[0][0].split(',')[0])
		else:
			obj = output_table[0]
			designation = db.query("select names from sources where id={}".format(obj[1]))
			spectype = u.specType(obj[3])
			grav = ('' if (obj[5]==None or obj[5]=='None') else obj[5])
			name = (io_names[individual_objs.index(obj[6])] if obj[6] in individual_objs else obj[6])
			subplots[k][1].plot(obj[4][0],obj[4][1],linestyle='-',color=colors[k],linewidth=2,label='{}{}, {}'.format(spectype,grav[:4],name))
				
			if plot_these:
				pcolors = ['blue','red','mediumturquoise','lime']
				for p in range(len(plot_these)):	
					obj_p = output_table[[n for n in range(len(output_table)) if output_table[n][6]==plot_these[p]]][0]
					if obj[6]==obj_p[6]:
						print 'best fit already plotted'
					else:
						spectype = u.specType(obj_p[3])
						grav = ('' if (obj_p[5]==None or obj_p[5]=='None') else obj_p[5])
						name = (io_names[individual_objs.index(obj_p[6])] if obj_p[6] in individual_objs else obj_p[6])
						subplots[k][1].plot(obj_p[4][0],obj_p[4][1],linestyle='--',color=pcolors[p],linewidth=1.5,label='{}{}, {}'.format(spectype,grav[:4],name))

				designations.append(designation[0][0].split(',')[0])
				subplots[k][1].legend(loc='best',frameon=False,fontsize = 'xx-small',handlelength=0,prop={'size':6})

		subplots[k][1].set_xlim(min(object[0])-0.01,max(object[0])+0.01)
		min_lim, max_lim = [],[]
		for i in range(len(object[1])):
			min_lim.append(object[1][i]-object[2][i])
			max_lim.append(object[1][i]+object[2][i])
		y_lims = (min(min_lim),max(max_lim))	
 		subplots[k][1].set_ylim(y_lims[0],y_lims[1])
 		subplots[k][1].set_yticklabels([])

	plt.subplots_adjust(wspace=0,hspace=0)	
	plt.xlabel('Wavelength ($\mu$m)',fontsize="large")
	plt.ylabel('Normalized Flux Density',fontsize="large",x=1,y=1)
	u.manual_legend(['Spectral Type'],['None'],['None'],fontsize=12,bbox_to_anchor=(0.36,3.3))
	u.manual_legend(['VHS 1256B', 'PSO 318','1741-4642','1516+3053'],['None','None','None','None'],['D','D','D','D'],edges=['blue','red','mediumturquoise','lime'],sizes=[4,4,4,4],fontsize=7,bbox_to_anchor=(.1,2.95),ncol=1)
	u.manual_legend(['field','young'],['gray','orange'],['o','o'],edges=['None','None'],sizes=[4,4],fontsize=7,bbox_to_anchor=(.01,2.95),ncol=1)
	subplots[k][1].legend(loc='bottom right',frameon=False,fontsize = 'xx-small',handlelength=0,prop={'size':6})

	plt.savefig('/Users/paigegiorla/Publications/HD2562B/Images/subplot.eps')
	plt.clf()	
	return designations