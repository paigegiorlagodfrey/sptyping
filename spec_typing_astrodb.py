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
import astropy.table as table

def scrub(data):
  '''
  For input data [w,f,e] or [w,f] returns the list with NaN, negative, and zero flux (and corresponsing wavelengths and errors) removed. 
  '''
  data = [i*q.Unit('') for i in data]
  data = [i[np.where(np.logical_and(data[1].value>0,~np.isnan(data[1].value)))] for i in data]
  data = [i[np.unique(data[0], return_index=True)[1]] for i in data]
  return [i[np.lexsort([data[0]])] for i in data]


def txtfile_obj(shortname=['1741-4642','1516+3053'],spt=[17.0,20.5],source_id=[1628,1423],g=['candidate','None'],filepath=['/Users/paigegiorla/Code/Data/1741-4642.txt','/Users/paigegiorla/Code/Data/SDSSJ1516+30.txt']):
	txt=[]
	for i in range(len(filepath)):
		object_file = open('{}'.format(filepath[i]), 'rb').readlines()
		W_obj,F_obj,U_obj = [],[],[]
		for line in object_file:
			columns=line.split()
			wav=float(columns[0])
			flu=float(columns[1])
			un=float(columns[2])
			W_obj.append(wav) ; F_obj.append(flu) ; U_obj.append(un)
		spec = [W_obj,F_obj,U_obj]
		txt.append({'s.id':'txt','s.source_id':source_id[i],'c.shortname':shortname[i],'t.spectral_type':spt[i],'t.gravity':g[i],'s.spectrum':spec})
	return txt

  
def showme(filepath,short_name,extraction,spt_range,norm_to,plot=False,reduced=False, M_type= False, L_type=False, T_type=True, Y_type=False,txtfile=False):
	'''
		This code will bin down all T dwarfs in database to the wavelength range and resolution of the object called in to question. 
		It returns the list of the spectral types of all of the T dwarf templates, and their corresponding reduced chi squared fit value.
	'''

## open file and read in the arrays for wavelength, flux, and uncertainty
	object_file=open('{}'.format(filepath), 'rb').readlines()
	W_obj,F_obj,U_obj,F_obj_norm,U_obj_norm,W_obj_microns=[],[],[],[],[],[]
	for line in object_file:
		columns=line.split()
		wav=float(columns[0])
		flu=float(columns[1])
		un=float(columns[2])
		W_obj.append(wav) ; F_obj.append(flu) ; U_obj.append(un)
	if W_obj[0]>10:
		W_obj = [W_obj[i]/1000. for i in range(len(W_obj))]

	object=[W_obj,F_obj,U_obj]
# 	object = m.normalize_to_band(norm_to[0],object)	
 	M, L, L_young, T, Y,txt = [],[],[],[],[],[]

	if txtfile==True:
			txt = txtfile_obj()
	if M_type==True:
			M = db.query("select s.id,s.source_id,c.shortname,t.spectral_type,t.gravity,s.spectrum from spectra s join spectral_types t on s.source_id=t.source_id join sources c on c.id=s.source_id where spectral_type>=0.0 and t.spectral_type<=9.0 and t.regime='IR' and t.adopted=1",fmt='dict')
	if L_type==True:
 			L = Jdb.query("select s.id,s.source_id,c.shortname,t.spectral_type,t.gravity,s.spectrum from spectra s join spectral_types t on s.source_id=t.source_id join sources c on c.id=s.source_id where t.spectral_type between 10.0 and 20.0 and t.adopted=1 and s.regime='NIR'",fmt='dict')
			source_ids = []
			for i in L:
				source_ids.append(i['s.source_id'])
 			L_y = Jdb.query("select s.id,s.source_id,c.shortname,t.spectral_type,t.gravity,s.spectrum from spectra s join spectral_types t on s.source_id=t.source_id join sources c on c.id=s.source_id where t.spectral_type between 10.0 and 25.0 and s.regime='NIR' and s.source_id not in ({})".format(','.join(map(str,source_ids))),fmt='dict')
			vhs = db.query("select s.id,s.source_id,c.shortname,t.spectral_type,t.gravity,s.spectrum from spectra s join spectral_types t on s.source_id=t.source_id join sources c on c.id=s.source_id where s.id=2240",fmt='dict')
# 			wise1741 = db.query("select s.id,s.source_id,c.shortname,t.spectral_type,t.gravity,s.spectrum from spectra s join spectral_types t on s.source_id=t.source_id join sources c on c.id=s.source_id where s.id=2703",fmt='dict')
			L_young = L_y + vhs 
	if T_type==True:
			T = db.query("select s.id,s.source_id,c.shortname,t.spectral_type,t.gravity,s.spectrum from spectra s join spectral_types t on s.source_id=t.source_id join sources c on c.id=s.source_id where t.spectral_type between 20.0 and 30.0 and t.regime='IR' and t.adopted=1 and s.instrument_id=6 and s.mode_id=1",fmt='dict')
	if Y_type==True:
			Y = db.query("select s.id,s.source_id,c.shortname,t.spectral_type,t.gravity,s.spectrum from spectra s join spectral_types t on s.source_id=t.source_id join sources c on c.id=s.source_id where t.spectral_type between 30.0 and 39.0 and t.regime='IR' and s.telescope_id=5",fmt='dict')

	meta_data = M + L + L_young + T + Y	+ txt

	dupes = [i for pos, i in enumerate(meta_data) if i['s.source_id'] in [j['s.source_id'] for j in meta_data[pos+1:]]]
	for dupe in dupes:
		meta_data.remove(dupe)
	print "objects in db = ", len(meta_data)	
	chisquare, chilist, namelist, specidlist, sptlist,templist,removelist,gravlist,output,shortnames = [],[],[],[],[],[],[],[],[],[]

	starting_wl=object[0][0]
	jump=np.abs((object[0][0]-object[0][1])/2)
	num_steps=len(object[0])
	print starting_wl, jump, num_steps

	object=scrub(object)	
	dof=len(object[0])-2-1	
	print "dof=", dof																#DOF is num of wavelength points minus free params (teff and logg) minus 1

## loop through each object
	for item in meta_data:
		if item['s.spectrum']==None:
			continue	
		try:
			if isinstance(item['s.spectrum'], astrodb.Spectrum):
				if len(item['s.spectrum'].data)==2:
					wave,flux =  item['s.spectrum'].data
					e = flux/10
					item['s.spectrum'].data = np.vstack((item['s.spectrum'].data,e))
				W,F,U = item['s.spectrum'].data
			else:
				W,F,U = item['s.spectrum']
				
			W,F,U = scrub([W,F,U])

			temp = [W,F,U]

			template = u.rebin_spec(temp,object[0])
			template = [k.value for k in template]
			if template[1][0]==0.0:
				print 'temp', item['s.source_id']
				continue
	## plot template over template's original flux
			if plot==True:
				plt.plot(temp[0],temp[1])
				plt.plot(template[0],template[1])			
				plt.savefig('/Users/paigegiorla/Publications/'+'{}'.format(short_name)+'/Images/original_over_template/{}'.format(extraction)+'/{}'.format(item['s.source_id'])+'.pdf')
				plt.clf()		

			template = m.normalize_with_c(template,object,norm_to[1])

	## plot template over object
			if plot==True:
				plt.scatter(template[0],template[1], color='red')
				plt.scatter(object[0],object[1], color='blue')
				plt.savefig('/Users/paigegiorla/Publications/'+'{}'.format(short_name)+'/Images/template_over_object/{}'.format(extraction)+'/{}'.format(item['s.source_id'])+'.pdf')
				plt.clf()

			if reduced:
				chi=sum(((template[1]-object[1])**2)/((template[2]**2)+(object[2]**2)))/dof     	#both data sets have uncertainties
			else:
				chi=sum(((template[1]-object[1])**2)/((template[2]**2)+(object[2]**2)))     	#both data sets have uncertainties
			namelist.append(item['s.source_id'])
			sptlist.append(item['t.spectral_type'])
			chilist.append(float(chi))
			specidlist.append(item['s.id'])
			templist.append(template)
			gravlist.append(item['t.gravity'])
			shortnames.append(item['c.shortname'])
		except (ValueError, AttributeError): pass
	output = [chilist,namelist,specidlist,sptlist,templist,gravlist,shortnames]
	output_table = table.Table([chilist,namelist,specidlist,sptlist,templist,gravlist,shortnames],names=('chi_values','source_id','spectra_id','spectral_types','template','gravity','shortnames'), meta={'name':'results from GOF'})
	output_table.sort('chi_values')
# 	dict = {'chi_values','source_id','spectra_id','spectral_types','template','gravity','shortnames'
# 	output_table.write('/Users/paigegiorla/Publications/'+'{}'.format(short_name)+'/Results/{}'.format(extraction)+'.csv',format='csv',fast_writer=False)
	fb = open('/Users/paigegiorla/Publications/'+'{}'.format(short_name)+'/Results/{}'.format(extraction)+'.pkl','wb')
	pickle.dump(output,fb)
	fb.close()
	chisquare=[sptlist,chilist]	
	print "# of objects fit = ", len(output_table)
	print "min chi value = ", output_table[0][0]
	print "spectral type = ", output_table[0][3]
	print "source id = ", output_table[0][1]
	print "spectrum id = ", output_table[0][2]
	print "young?", output_table[0][5]
	
	return chisquare, object, output_table
