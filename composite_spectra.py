import splat

## put all of the standard spectrum objects in a predefined dictionary called SPEX_STDS with keys corresponding to spectral type
## called SPEX_STDS with keys corresponding to the spectral type; e.g., splat.SPEX_STDS[‘M5.0’]  <— this is a spectrum object for the M5 standard
splat.initiateStandards()

## flux calibrate all the standards to their absolute magnitudes, which is done with typeToMag to get the abs mag and fluxCalibrate to scale the spectrum
## do this in a for loop to calibrate all of the spectra in the standard list
mag, mag_e = splat.typeToMag('L5.0','2MASS J',ref='filippazzo’)
sp = splat.SPEX_STDS['L5.0']
sp.fluxCalibrate('2MASS J',mag,absolute=True)

## make templates by adding the spectra together
sp1 = splat.SPEX_STDS['L5.0']
sp2 = splat.SPEX_STDS['T5.0’]
sp = sp1 + sp2

## Then you can use compareSpectra to compare your spectrum to the templates (presumably created in a big array of some kind)