from __future__ import division
from astrodbkit import astrodb
db=astrodb.Database('/Users/paigegiorla/Dropbox/BDNYCdb/BDNYCdev.db')
from matplotlib import pyplot as plt 
import numpy as np
from BDNYCdb import utilities as ut
import astropy.units as q
import scipy.stats as s
import small_functions as m
from matplotlib import gridspec
from itertools import cycle
import astrotools
import pickle
from astropy.table import Table
import glob

def scrub(data):
  '''
  For input data [w,f,e] or [w,f] returns the list with NaN, negative, and zero flux (and corresponsing wavelengths and errors) removed. 
  '''
  data = [i*q.Unit('') for i in data]
  data = [i[np.where(np.logical_and(data[1].value>0,~np.isnan(data[1].value)))] for i in data]
  data = [i[np.unique(data[0], return_index=True)[1]] for i in data]
  return [i[np.lexsort([data[0]])] for i in data]
  
def showme(filefolder,filepath,short_name,extraction,spt_range,plot=False,reduced=False):
	'''
		This code will bin down all T dwarfs in database to the wavelength range and resolution of the object called in to question. 
		It returns the list of the spectral types of all of the T dwarf templates, and their corresponding reduced chi squared fit value.
		
		filefolder: path to folder with text files of templates
		filepath: path to text file of object in question
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
 	

	template_folder = glob.glob(filefolder)
	template_list = []
	for file in template_folder:
		w_temp,f_temp,u_temp = [],[],[]
		open_file = open(file).readlines()
		for line in open_file:
			columns=line.split()
			w = float(columns[0])
			f = float(columns[1])
			un = float(columns[2])
			w_temp.append(w),f_temp.append(f),u_temp.append(un)			
		dict={'spectral_type':file.split('/')[-1].split('.')[0].split('_')[0],'gravity':file.split('/')[-1].split('.')[0].split('_')[2], 'spectrum':[w_temp,f_temp,u_temp]}
		template_list.append(dict)	
	
	
#!	loop through each object
	chisquare, chilist, namelist, specidlist, sptlist,templist,removelist,gravlist,output,shortnames = [],[],[],[],[],[],[],[],[],[]

	starting_wl=object[0][0]
	jump=np.abs((object[0][0]-object[0][1])/2)
	num_steps=len(object[0])
	print starting_wl, jump, num_steps

	object=scrub(object)	
	dof=len(object[0])-2-1	
	print "dof=", dof																#DOF is num of wavelength points minus free params (teff and logg) minus 1
	for item in template_list:
		try:
			W,F,U = item['spectrum']
			W,F,U = scrub([W,F,U])

			maxF = max(F)
			F = F/maxF
			U = U/maxF
			temp = [W,F,U]
			template = ut.rebin_spec(temp,object[0])
			template = [k.value for k in template]
	#!		plot template over template's original flux
			if plot==True:
				plt.plot(temp[0],temp[1],label='original')
				plt.plot(template[0],template[1],label='template')
				plt.legend()			
				plt.savefig('/Users/paigegiorla/Publications/'+'{}'.format(short_name)+'/Images/original_over_template/{}'.format(extraction)+'/{}_{}'.format(item['spectral_type'],item['gravity'])+'.pdf')
				plt.clf()		
			a=float(sum(object[1]*template[1]/((template[2]**2)+(object[2]**2))))
			b=float(sum(template[1]*template[1]/((template[2]**2)+(object[2]**2))))
			c=a/b
			template[1]=template[1]*c
			template[2]=template[2]*c

	#!		plot template over object
			if plot==True:
				plt.scatter(template[0],template[1], color='red',label='template')
				plt.scatter(object[0],object[1], color='blue',label='object')
				plt.xlim(1.14,1.8)
				plt.legend()
				plt.savefig('/Users/paigegiorla/Publications/'+'{}'.format(short_name)+'/Images/template_over_object/{}'.format(extraction)+'/{}_{}'.format(item['spectral_type'],item['gravity'])+'.pdf')
				plt.clf()

			if reduced:
				chi=sum(((template[1]-object[1])**2)/((template[2]**2)+(object[2]**2)))/dof     	#both data sets have uncertainties
			else:
				chi=sum(((template[1]-object[1])**2)/((template[2]**2)+(object[2]**2)))     	#both data sets have uncertainties
			spt = astrotools.specType(item['spectral_type'])
			sptlist.append(spt)
			chilist.append(float(chi))
			templist.append(template)
			gravlist.append(item['gravity'])
		except AttributeError: pass
	output = [chilist,sptlist,templist,gravlist]
	output_table = Table([chilist,sptlist,templist,gravlist],names=('chi_values','spectral_types','template','gravity'), meta={'name':'results from GOF'})
	output_table.sort('chi_values')
	pickle.dump(output_table,open('/Users/paigegiorla/Publications/'+'{}'.format(short_name)+'/Results/{}'.format(extraction)+'.pkl','wb'))
	chisquare=[sptlist,chilist]	
	print "# of objects fit = ", len(output_table)
	print "min chi value = ", output_table[0][0]
	print "spectral type = ", output_table[0][1]
	print "young?", output_table[0][3]
	
	return chisquare, object, output_table
