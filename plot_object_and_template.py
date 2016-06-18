from astrodbkit import astrodb
db=astrodb.Database('/Users/paigegiorla/Dropbox/BDNYCdb/BDNYCdev.db')
from matplotlib import pyplot as plt
from BDNYCdb import utilities as u
import numpy as np 
import small_functions as m
import astropy.units as q

def scrub(data):
  '''
  For input data [w,f,e] or [w,f] returns the list with NaN, negative, and zero flux (and corresponsing wavelengths and errors) removed. 
  '''
  data = [i*q.Unit('') for i in data]
  data = [i[np.where(np.logical_and(data[1].value>0,~np.isnan(data[1].value)))] for i in data]
  data = [i[np.unique(data[0], return_index=True)[1]] for i in data]
  return [i[np.lexsort([data[0]])] for i in data]

def showme(obj_name,band,filepath,db_spectra_ids,colors,facecolors):
	'''obj_spectrum: [ name, wavelength, flux]
		 db_spectra_ids: list of ids from astrodb to plot against obj
	'''
	object_file=open('{}'.format(filepath), 'rb').readlines()
	W_obj,F_obj,U_obj,F_obj_norm,U_obj_norm,W_obj_microns=[],[],[],[],[],[]
	for line in object_file:
		columns=line.split()
		wav=float(columns[0])
		flu=float(columns[1])
		un=float(columns[2])
		W_obj.append(wav) ; F_obj.append(flu) ; U_obj.append(un)
		
	maxFobj=max(F_obj)
	F_obj=[F_obj[i]/maxFobj for i in range(len(F_obj))]
	U_obj=[U_obj[i]/maxFobj for i in range(len(U_obj))]
	object=[W_obj,F_obj,U_obj]
 	object=scrub(object)	

	i=0
	for id in db_spectra_ids:
		comp_spectrum = db.query("select s.spectrum,so.shortname from spectra as s join sources as so on s.source_id=so.id where s.id={}".format(id),fmt='dict')
		W,F,U = comp_spectrum[0]['s.spectrum'].data
		W,F,U = scrub([W,F,U])
		maxF = max(F)
		F = F/maxF
		U = U/maxF
		temp = [W,F,U]
		template = u.rebin_spec(temp,object[0])
		template = [k.value for k in template]
		a=float(sum(object[1]*template[1]/((template[2]**2)+(object[2]**2))))
		b=float(sum(template[1]*template[1]/((template[2]**2)+(object[2]**2))))
		c=a/b
		template[1]=template[1]*c
		template[2]=template[2]*c
	 	plt.plot(template[0],template[1],color=colors[i],label=comp_spectrum[0]['so.shortname'])
	 	plt.fill_between(template[0], template[1]-template[2], template[1]+template[2],linewidth=0,facecolor=facecolors[i])
	 	i+=1
	 
	plt.errorbar(object[0],object[1],object[2],color='k',label=obj_name)
	plt.legend(loc='best')
	plt.ylim(0.01,1.01)
	plt.xlim(object[0][0]-0.01,object[0][-1]+0.01)

	plt.savefig('/Users/paigegiorla/Publications/GOI28/Images/comp_to_Ls_{}'.format(band)+'.eps')
	plt.clf()
	
def bands_subplot(obj_name,db_spectra_ids):
	'''obj_spectrum: [ name, wavelength, flux]
		 db_spectra_ids: list of ids from astrodb to plot against obj
	'''
	filepaths = ['/Users/paigegiorla/Code/Data/GOI28/J_Spectrum_GOI28_JHK1K2_v5.ascii','/Users/paigegiorla/Code/Data/GOI28/H_Spectrum_GOI28_JHK1K2_v5.ascii','/Users/paigegiorla/Code/Data/GOI28/K1K2_Spectrum_GOI28_JHK1K2_v5.ascii']
	i=0
	bands = []
	for filepath in filepaths:
		object_file=open('{}'.format(filepath), 'rb').readlines()
		W_obj,F_obj,U_obj,F_obj_norm,U_obj_norm,W_obj_microns=[],[],[],[],[],[]
		for line in object_file:
			columns=line.split()
			wav=float(columns[0])
			flu=float(columns[1])
			un=float(columns[2])
			W_obj.append(wav) ; F_obj.append(flu) ; U_obj.append(un)
		maxFobj=max(F_obj)
		F_obj=[F_obj[i]/maxFobj for i in range(len(F_obj))]
		U_obj=[U_obj[i]/maxFobj for i in range(len(U_obj))]
		object=[W_obj,F_obj,U_obj]
 		object=scrub(object)	
 		bands.append(object)
		i+=1
	fig, axes = plt.subplots(nrows=1, ncols=3,sharex=True,sharey=True)    # hide frame	
	colors = ['b','r']
	facecolors = ['#089FFF','#FF9848']

	for k,ax in zip(range(len(bands)),axes.flat):

		for m in range(len(db_spectra_ids)):
			comp_spectrum = db.query("select s.spectrum,so.shortname from spectra as s join sources as so on s.source_id=so.id where s.id={}".format(db_spectra_ids[m]),fmt='dict')
			W,F,U = comp_spectrum[0]['s.spectrum'].data
			W,F,U = scrub([W,F,U])
			maxF = max(F)
			F = F/maxF
			U = U/maxF
			temp = [W,F,U]
			template = u.rebin_spec(temp,bands[k][0])
			template = [n.value for n in template]
			a=float(sum(bands[k][1]*template[1]/((template[2]**2)+(bands[k][2]**2))))
			b=float(sum(template[1]*template[1]/((template[2]**2)+(bands[k][2]**2))))
			c=a/b
			template[1]=template[1]*c
			template[2]=template[2]*c
			ax.scatter(template[0],template[1],color=colors[m])
			ax.fill_between(template[0], template[1]-template[2], template[1]+template[2],linewidth=0,facecolor=facecolors[m])
			ax.annotate(str(comp_spectrum[0]['so.shortname']),xy=(2, 0.2+(m/2)),color=colors[m])

		ax.errorbar(bands[k][0],bands[k][1],bands[k][2],color='k')
		ax.annotate(str(obj_name),xy=(2, 0.15),color='k')
		ax.set_ylim(0.01,1.01)
		ax.set_xlim(bands[k][0][0]-0.01,bands[k][0][-1]+0.01)

	plt.savefig('/Users/paigegiorla/Publications/GOI28/Images/comp_to_Ls_subplots.eps')
	plt.clf()
	
