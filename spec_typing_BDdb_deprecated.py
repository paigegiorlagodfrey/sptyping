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
from astropy.table import Table

def scrub(data):
  '''
  For input data [w,f,e] or [w,f] returns the list with NaN, negative, and zero flux (and corresponsing wavelengths and errors) removed. 
  '''
  data = [i*q.Unit('') for i in data]
  data = [i[np.where(np.logical_and(data[1].value>0,~np.isnan(data[1].value)))] for i in data]
  data = [i[np.unique(data[0], return_index=True)[1]] for i in data]
  return [i[np.lexsort([data[0]])] for i in data]
  
def showme(filepath,short_name,extraction,plot=False, M_type= False, L_type=False, T_type=True, Y_type=False,):
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


#!	get source_id for all T dwarfs from database, turn it into a list and use list to get spectra_id for each source's spex data	
 	M, L, T, Y = [],[],[],[]
	if M_type==True:
			M = db.list("select distinct source_id from spectral_types where spectral_type>=0 and spectral_type<=9 and regime='IR' and adopted=1").fetchall()
	if L_type==True:
			L = db.list("select distinct source_id from spectral_types where spectral_type>=10 and spectral_type<=19 and regime='IR' and adopted=1").fetchall()
	if T_type==True:
			T = db.list("select distinct source_id from spectral_types where spectral_type>=20 and spectral_type<=29 and regime='IR' and adopted=1").fetchall()
	if Y_type==True:
			Y = db.list("select distinct source_id from spectral_types where spectral_type>=30 and spectral_type<=39 and regime='IR' and adopted=1").fetchall()
	sources = M + L + T +  Y

# 	sources=db.list("select distinct source_id from spectral_types where spectral_type>=20 and spectral_type<=30 and regime='IR' and adopted=1").fetchall()
	sourcelist, ids, idlist=[],[],[]
	index=0
	while index<len(sources):
		sourcelist.append(sources[index][0])
		index=index+1
	for i in sourcelist:
		id=db.list("SELECT distinct id FROM spectra WHERE source_id='{}' AND instrument_id=6 and mode_id=1".format(i)).fetchone()
		ids.append(id)

	ids=filter(None,ids)	
	index=0
	while index<len(ids):
		idlist.append(ids[index][0])
		index=index+1		

#!	loop through each object
	chisquare, chilist, namelist, specidlist, sptlist,templist,removelist,gravlist,shortnames=[],[],[],[],[],[],[],[],[]

	starting_wl=object[0][0]
	jump=np.abs((object[0][0]-object[0][1])/2)
	num_steps=len(object[0])
	print starting_wl, jump, num_steps

	object=scrub(object)	
	dof=len(object[0])-2-1	
	print "dof=", dof																#DOF is num of wavelength points minus free params (teff and logg) minus 1
	print "# of objects", len(idlist)
	for j in idlist:
		name,W,F,U=db.list("SELECT source_id,wavelength,flux,unc from spectra where id={}".format(j)).fetchone()
		W,F,U=u.scrub([W,F,U])
		maxF=max(F)
		F=F/maxF
		U=U/maxF
		temp=[W, F, U]
		template=u.rebin_spec(temp,object[0])
  		template=[k.value for k in template]

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

# 		source,spectype,gravity=db.list("SELECT spectral_types.source_id, spectral_types.spectral_type,spectral_types.gravity from spectral_types join spectra on spectral_types.source_id=spectra.source_id where spectral_types.regime='IR' and spectral_types.adopted=1 and spectral_types.spectral_type>=20 and spectral_types.spectral_type<=30 and spectra.id='{}'".format(j)).fetchone()
		source,spectype,gravity,shortname = db.list("SELECT spectral_types.source_id, spectral_types.spectral_type,spectral_types.gravity,sources.shortname from spectral_types join spectra on spectral_types.source_id=spectra.source_id join sources on sources.id=spectral_types.source_id where spectral_types.adopted=1 and spectra.id='{}'".format(j)).fetchone()
#!		plot template over object
		if plot==True:
			plt.scatter(template[0],template[1], color='red')
			plt.scatter(object[0],object[1], color='blue')
			plt.xlim(1.14,1.8)
			plt.savefig('/Users/paigegiorla/Publications/'+'{}'.format(short_name)+'/Images/template_over_object/{}'.format(extraction)+'/{}'.format(j)+'.pdf')
			plt.clf()

		chi=sum(((template[1]-object[1])**2)/((template[2]**2)+(object[2]**2)))/dof     	#both data sets have uncertainties
		namelist.append(source)
		sptlist.append(spectype)
		chilist.append(float(chi))
		specidlist.append(j)
		templist.append(template)
		gravlist.append(gravity)
		shortnames.append(shortname)
 	output = [chilist,namelist,specidlist,sptlist,templist,gravlist,shortnames]
# 	return output
	output_sorted = zip(*sorted(zip(*output)))
	output_table = Table([output_sorted[0],output_sorted[1],output_sorted[2],output_sorted[3],output_sorted[4],output_sorted[5],output_sorted[6]],names=('chi_values','source_id','spectra_id','spectral_types','template','gravity','shortnames'), meta={'name':'results from GOF'})
	pickle.dump(output,open('/Users/paigegiorla/Publications/'+'{}'.format(short_name)+'/Results/{}'.format(extraction)+'.pkl','wb'))
	chisquare=[sptlist,chilist]	
	print "min chi value = ", output_table[0][0]
	print "spectral type = ", output_table[0][3]
	print "source id = ", output_table[0][1]
	print "spectrum id = ", output_table[0][2]
	print "young?", output_table[0][5]

	
	return chisquare, object, output_table
