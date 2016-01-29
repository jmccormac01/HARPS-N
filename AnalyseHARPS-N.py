#
# playing with HARPS-N spectra
# simple script to have a fiddle with HARPS data
#

# Notes from Gzregorz Nowak on HARPS-N DRS pipeline output
"""
The data were reduced automatically, directly after observations by
the DRS pipeline and the spectra were cross-correlated using G2V
template mask. Therefore you may find in the scientific_data/
subdirectory not only merged stellar (s1d_A) and sky (s1d_B) spectra
(as well as spectra order by order, i.e. e2ds_A/B FITS files), but
also FITS files with cross-correlation functions (ccf_G2_A/B), and
bisector lines (bis_G2_A/B) and data files with Gaussian parameters
fitted to the CCFs from different orders (_ccf_G2_A/B.tbl).
"""

"""
Questions:
	How exactly is the final combined CCF made? Is it a simple summation?
	Have there been SNR limits applied? 
	There seems to be no Gaussian in the sky CCF, so the fit fails.
	Another way to correct the affected RV points?
"""

from scipy.optimize import curve_fit
import matplotlib.pyplot as pl
from astropy.io import fits
from math import sqrt,ceil
import glob as g
import numpy as np
import seaborn

# gamma_velocity for the planet we are analysing
gamma_vel= -6.68149

# Gaussian = Aexp[ -(x-x0)^2 / (2*width^2) ]
def Gaussian(x, *p): 
	A, x0, width = p
	y = np.zeros(x.shape)
	g=A*np.exp((-(x-x0)**2)/(2*width**2))
	return g

# return estimated A, x0 and fwhm
def estimateGaussianParams(t_data,xn):
	return [min(t_data) - max(t_data),np.median(xn),np.std(xn)] 

# code below is for analysing sky or star CCFs
def analyseCCFs(ccf_G2):
	amp=np.empty(len(ccf_G2))
	rv=np.empty(len(ccf_G2))
	fwhm=np.empty(len(ccf_G2))
	mjd=np.empty(len(ccf_G2))
	combined_ccfs=[]
	for i in range(0,len(ccf_G2)):
		h=fits.open(ccf_G2[i])
		crval1=h[0].header['CRVAL1']
		cdelt1=h[0].header['CDELT1']
		mjd[i]=h[0].header['MJD-OBS']
		data=h[0].data
		
		# prepare the plotting array
		fig = pl.figure(i+1,figsize=(15,15))
		fig.clf()
		fig.suptitle('%s CCFs' % (ccf_G2[i]))
		seaborn.axes_style("darkgrid")
		montage_dim=int(ceil(sqrt(len(data))))
		c=0
		for j in range(0,montage_dim):
			for k in range(0,montage_dim):
				ax = fig.add_subplot(montage_dim, montage_dim, c+1, xticks=[], yticks=[])		
				if c < len(data):
					# try fitting a Gaussian to the ccf
					# tweak y values to help the fit
					# only fit those that are not 0s
					x=np.arange(0,len(data[c]))
					xn=(x+cdelt1)+crval1
					t_data=data[c]-max(data[c])	
					ax.plot(xn,t_data,'b.')
					ax.set_title('order %d' % (c+1))
					if max(data[c]) > 0 and np.std(data[c]) > 0:
						p0=estimateGaussianParams(t_data,xn)
						try:
							coeff,cov=curve_fit(Gaussian,xn,t_data,p0)
							yfit = Gaussian(xn, *coeff)
							print ('[%d/%d] x0=%.4f ' % (c+1,len(data),coeff[1]))
							ax.plot(xn,yfit,'r--')
							if c == len(data)-1:
								amp[i]=coeff[0]
								rv[i]=coeff[1]
								fwhm[i]=coeff[2]
								combined_ccfs.append(data[c])
						except RuntimeError:
							print ('[%d/%d] fit did not converge ' % (c+1,len(data)))
							if c == len(data)-1:
								print ('[%d/%d] FINAL FIT DID NOT CONVERGE! ' % (c+1,len(data)))
								combined_ccfs.append(data[c])
					else:
						print('[%d/%d] Cannot fit zeros' % (c+1,len(data)))
					c+=1
					#print("[%d/%d]" % (c,len(data)))
		# plot the montage of CCFs
		fig.savefig('%s.png' % (ccf_G2[i]),dpi=300)
	return mjd,rv_corr,amp,fwhm,combined_ccfs

# make a plot of the RV results
def makeSummaryPlot(mjd,rv,amp,fwhm):
	# plot the results from the analysis of the RVs/RM
	mjds=mjd-int(mjd[0])
	fig = pl.figure(i+1,figsize=(15,15))
	fig.clf()
	ax = fig.add_subplot(3,1,1)
	ax.plot(mjds,rv,'ro')
	ax.set_ylabel('RV (km/s))')
	ax = fig.add_subplot(3,1,2)
	ax.plot(mjds,amp,'ro')
	ax.set_ylabel('Inverse CCF SNR')
	ax = fig.add_subplot(3,1,3)
	ax.plot(mjds,fwhm,'ro')
	ax.set_xlabel('MJD - int(MJD[0])')
	ax.set_ylabel('FWHM (km/s))')
	pl.savefig('SummaryResults.png',dpi=300)

if __name__ == '__main__':
	# get lists of the different HARPS-N DRS
	stellar_1d=g.glob('*_s1d_A.fits')
	sky_1d=g.glob('*_s1d_B.fits')
	stellar_orders=g.glob('*_e2ds_A.fits')
	sky_orders=g.glob('*_e2ds_B.fits')
	stellar_ccf_G2=g.glob('*ccf_G2_A.fits')
	sky_ccf_G2=g.glob('*ccf_G2_B.fits')
	stellar_bis_G2=g.glob('*bis_G2_A.fits')
	sky_bis_G2=g.glob('*bis_G2_B.fits')
	stellar_gauss_ccf_G2=g.glob('*ccf_G2_A.tbl')
	sky_gauss_ccf_G2=g.glob('*ccf_G2_B.tbl')

	# analyse the stellar spectra
	mjd_star,rv_star,amp_star,fwhm_star,ccfs_star=analyseCCFs(stellar_ccf_G2)
	mjd_sky,rv_sky,amp_sky,fwhm_sky,ccfs_sky=analyseCCFs(sky_ccf_G2)

	
	
	
	
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

	
	
	
	
	
	

