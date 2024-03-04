#============================================================
#	IMAGE TRIM - WCSREMAP - SUBTRACTION
#------------------------------------------------------------
#	2023.04.28	UPDATED BY GREGORY S.H. PAEK
#	2020.02.16	CREATED BY GREGORY S.H. PAEK
#============================================================
import os
import sys
import glob
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
import numpy as np
#============================================================
#	FUNCTION
#============================================================
def trim(inim, position, size, outim='trim.fits'):
	# Load the image and the WCS
	hdu = fits.open(inim)[0]
	wcs = WCS(hdu.header)
	# Make the cutout, including the WCS
	cutout = Cutout2D(
		hdu.data,
		position=position,
		size=size,
		wcs=wcs,
		mode='partial',
		fill_value=np.median(hdu.data),
		)
	# Put the cutout image in the FITS HDU
	hdu.data = cutout.data
	# Update the FITS header with the cutout WCS
	hdu.header.update(cutout.wcs.to_header())
	# Write the cutout to a new FITS file
	hdu.writeto(outim, overwrite=True)
#------------------------------------------------------------
def wcsremap(inim, refim, outim, path_com='/data3/wcsremap/wcsremap-1.0.1/wcsremap'):
	com = f'{path_com} -template {refim} -source {inim} -outim {outim}'
	print(com)
	os.system(com)
	return outim
#------------------------------------------------------------
def hotpants(inim, refim, outim='hd.fits', convim='hc.fits'):
	'''
	inim : Science image
	refim : Reference image
	'''
	com = f'hotpants -c t -n i -iu 6000000000 -il -100000 -tu 6000000000 -tl -100000 -v 0 -inim {inim} -tmplim {refim} -outim {outim} -oci {convim}'
	print(com)
	os.system(com)
#------------------------------------------------------------
def trim_sub_routine(inim, refim, position, size):
	if os.path.dirname(inim) == '':
		inim = f'./{inim}'
	#	Output names
	path = os.path.dirname(inim)
	image = os.path.basename(inim)
	outim = f'{path}/hd{image}'
	convim = f'{path}/hc{image}'
	trinim = f'{path}/tr{image}'
	trrefim = f'{path}/tr{os.path.basename(refim)}'
	wrim = f'{path}/wr{os.path.basename(refim)}'
	#	Trim
	trim(inim, position, size, trinim)
	trim(refim, position, size, trrefim)
	#	WCSremap
	wcsremap(inim=trrefim, refim=trinim, outim=wrim)
	#	Hotpants
	hotpants(trinim, wrim, outim=outim, convim=convim)
	#	ds9 command
	ds9com = f'ds9 {trinim} {wrim} {outim} -lock frame wcs -tile column &'
	return ds9com
#============================================================
#	IMSNG
#------------------------------------------------------------
os.system('ls *.fits')
print('='*60)
print('Input paramters')
print('-'*60)
# imkey = input('Scince image: ')
# if imkey == '':
# 	imkey = 'Calib*m.fits'
# imkey = './Calib*m.fits'
# imkey = './Calib-LOAO-NGC7769-20230209-022233-R-180.com.fits'
# imkey = sys.argv[1]
imkey = './C*m.fits'
imlist = sorted(glob.glob(imkey))
# refim = input('Reference image: ')
# if refim == '':
# 	refim = 'ref.fits'
# refim = './cutout_rings.v3.skycell.1847.000.stk.r.unconv.fits'
# refim = sys.argv[2]
refim = './ref.fits'
#------------------------------------------------------------
#	Exception for files in the current directory
if os.path.dirname(refim) == '':
	refim = f'./{refim}'
#------------------------------------------------------------
print('-'*60)
print('Image trim paramters')
print('-'*60)
#------------------------------------------------------------
# tra = float(input('RA [deg]\t: '))
# tdec = float(input('Dec [deg]\t: '))
# length = float(input('Length [arcmin]\t: '))
# tra, tdec = 60.858, -75.379
# tra, tdec = 357.767,20.150
# length = 6
# tra, tdec = float(sys.argv[3]), float(sys.argv[4])
# length = float(sys.argv[5])
tra, tdec = 210.910674637, +54.3116510708
length = 10

#------------------------------------------------------------
position = SkyCoord(tra, tdec, frame='icrs', unit='deg')
size = u.Quantity((length*2, length*2), u.arcmin)
#============================================================
ds9comlist = []
for inim in imlist:
	ds9com = trim_sub_routine(inim, refim, position, size)
	ds9comlist.append(ds9com)
