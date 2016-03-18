from __future__ import division
from astrodbkit import astrodb
db=astrodb.Database('/Users/paigegiorla/Dropbox/BDNYCdb/BDNYCv1.0.db')
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

def scrub(data):
  '''
  For input data [w,f,e] or [w,f] returns the list with NaN, negative, and zero flux (and corresponsing wavelengths and errors) removed. 
  '''
  data = [i*q.Unit('') for i in data]
  data = [i[np.where(np.logical_and(data[1].value>0,~np.isnan(data[1].value)))] for i in data]
  data = [i[np.unique(data[0], return_index=True)[1]] for i in data]
  return [i[np.lexsort([data[0]])] for i in data]
  
def showme(filepath,short_name,extraction,spt_range,plot=False):
	'''
		This code will bin down all T dwarfs in database to the wavelength range and resolution of the object called in to question. 
		It returns the list of the spectral types of all of the T dwarf templates, and their corresponding reduced chi squared fit value.
	'''

#!	open file and read in the arrays for wavelength, flux, and uncertainty
	object_file=open('{}'.format(filepath), 'rb').readlines()
	print short_name
	W_obj,F_obj,U_obj,F_obj_norm,U_obj_norm,W_obj_microns=[],[],[],[],[],[]
	for line in object_file:
		columns=line.split()
		wav=float(columns[0])
		flu=float(columns[1])
		un=float(columns[2])
		W_obj.append(wav) ; F_obj.append(flu) ; U_obj.append(un)
	if W_obj[0]>10:
		W_obj = [W_obj[i]/1000. for i in range(len(W_obj))]

	maxFobj=max(F_obj)
	F_obj=[F_obj[i]/maxFobj for i in range(len(F_obj))]
	U_obj=[U_obj[i]/maxFobj for i in range(len(U_obj))]
	object=[W_obj,F_obj,U_obj]

	meta_data = db.query("select s.id,s.source_id,c.shortname,t.spectral_type,t.gravity,s.spectrum from spectra s join spectral_types t on s.source_id=t.source_id join sources c on c.id=s.source_id where t.spectral_type between ? and ? and t.adopted=1 and s.instrument_id=6 and s.telescope_id=7 and s.mode_id=1",params=(spt_range[0],spt_range[1]),fmt='dict')	
	dupes = [i for pos, i in enumerate(meta_data) if i['s.source_id'] in [j['s.source_id'] for j in meta_data[pos+1:]]]
	for dupe in dupes:
		meta_data.remove(dupe)
	print "objects in db = ", len(meta_data)	
#!	loop through each object
	chisquare, chilist, namelist, specidlist, sptlist,templist,removelist,gravlist,output,shortnames = [],[],[],[],[],[],[],[],[],[]

	starting_wl=object[0][0]
	jump=np.abs((object[0][0]-object[0][1])/2)
	num_steps=len(object[0])
	print starting_wl, jump, num_steps

	object=scrub(object)	
	dof=len(object[0])-2-1	
	print "dof=", dof																#DOF is num of wavelength points minus free params (teff and logg) minus 1
	
	for item in meta_data:
		if item['s.spectrum']==None:
# 			print item['s.source_id']
			continue
		W,F,U = item['s.spectrum'].data
		W,F,U = u.scrub([W,F,U])

		maxF = max(F)
		F = F/maxF
		U = U/maxF
		temp = [W,F,U]
		template = u.rebin_spec(temp,object[0])
		template = [k.value for k in template]
		if template[1][0]==0.0:
 			print 'temp', item['s.source_id']
			continue
#!		plot template over template's original flux
		if plot==True:
			plt.plot(temp[0],temp[1])
			plt.plot(template[0],template[1])			
			plt.savefig('/Users/paigegiorla/Publications/'+'{}'.format(short_name)+'/Images/original_over_template/{}'.format(extraction)+'/{}'.format(j)+'.pdf')
			plt.clf()		
		a=float(sum(object[1]*template[1]/((template[2]**2)+(object[2]**2))))
		b=float(sum(template[1]*template[1]/((template[2]**2)+(object[2]**2))))
		c=a/b
		template[1]=template[1]*c
		template[2]=template[2]*c

#!		plot template over object
		if plot==True:
			plt.scatter(template[0],template[1], color='red')
			plt.scatter(object[0],object[1], color='blue')
			plt.xlim(1.14,1.8)
			plt.savefig('/Users/paigegiorla/Publications/'+'{}'.format(short_name)+'/Images/template_over_object/{}'.format(extraction)+'/{}'.format(j)+'.pdf')
			plt.clf()

		chi=sum(((template[1]-object[1])**2)/((template[2]**2)+(object[2]**2)))/dof     	#both data sets have uncertainties
		namelist.append(item['s.source_id'])
		sptlist.append(item['t.spectral_type'])
		chilist.append(float(chi))
		specidlist.append(item['s.id'])
		templist.append(template)
		gravlist.append(item['t.gravity'])
		shortnames.append(item['c.shortname'])
# 		output.append([chi,item['s.source_id'],item['t.spectral_type'],template,item['t.gravity']])
	output = [chilist,namelist,specidlist,sptlist,templist,gravlist]
	output_table = Table([chilist,namelist,specidlist,sptlist,templist,gravlist,shortnames],names=('chi_values','source_id','spectra_id','spectral_types','template','gravity','shortnames'), meta={'name':'results from GOF'})
	output_table.sort('chi_values')
	pickle.dump(output_table,open('/Users/paigegiorla/Publications/'+'{}'.format(short_name)+'/Results/{}'.format(extraction)+'.pkl','wb'))
	chisquare=[sptlist,chilist]	
	print "# of objects fit = ", len(output_table)
	print "min chi value = ", output_table[0][0]
	print "spectral type = ", output_table[0][3]
	print "source id = ", output_table[0][1]
	print "spectrum id = ", output_table[0][2]
	print "young?", output_table[0][5]
	
	return chisquare, object, output_table
