import splat

## put all of the standard spectrum objects in a predefined dictionary called SPEX_STDS with keys corresponding to spectral type
splat.initiateStandards()

## flux calibrate all the standards to their absolute magnitudes, which is done with typeToMag to get the abs mag and fluxCalibrate to scale the spectrum
## do this in a for loop to calibrate all of the spectra in the standard list
spectra = []
for d in splat.SPEX_STDS:
# 	mag, mag_e = splat.typeToMag(str(splat.SPEX_STDS),'2MASS J',ref='filippazzo')
	sp = splat.SPEX_STDS[d]
# 	sp.fluxCalibrate('2MASS J',mag,absolute=True)
	spectra.append(sp)
print spectra[0].keys()
print spectra[0].values()	


dict=[]
for i in range(len(spectra)):
	try:
		for x in range(len(spectra)):
			print x, spectra[i].keys()[0], spectra[i+x].values()[0]
			name = spectra[i].keys()[0] + spectra[i+x].values()[0]
			spec = spectra[i].keys()[0] + spectra[i+x].values()[0]
			dict.append({name:spec})
	except: pass	
# print dict[0]
		
## Then you can use compareSpectra to compare your spectrum to the templates (presumably created in a big array of some kind)