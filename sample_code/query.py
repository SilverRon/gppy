#============================================================
#	https://ps1images.stsci.edu/ps1image.html
#	By Sophia Kim 2019.01.22. based on code PS1 suggests on the link above
#	Pan-STARRS DR1 data query
#	from https://michaelmommert.wordpress.com/2017/02/13/accessing-the-gaia-and-pan-starrs-catalogs-using-python/
#	By CS Choi 
#	REVISED AND ORGANIZED BY GREGORY S.H. PAEK
#	UPDATE : 20.01.03
#============================================================
from __future__ import print_function
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.io import ascii
from astropy import units as u
from PIL import Image
from io import BytesIO
import numpy as np
import sys
import requests
import pylab
from astropy.io.votable import parse
import multiprocessing as mp
import os, glob
from imsng import phot
#============================================================
def geturl(ra, dec, size=240, output_size=None, filters="grizy", format="jpg", color=False):
	"""
	Get URL for images in the table
	
	ra, dec = position in degrees
	size = extracted image size in pixels (0.25 arcsec/pixel)
	output_size = output (display) image size in pixels (default = size).
				  output_size has no effect for fits format images.
	filters = string with filters to include
	format = data format (options are "jpg", "png" or "fits")
	color = if True, creates a color image (only for jpg or png format).
			Default is return a list of URLs for single-filter grayscale images.
	Returns a string with the URL
	"""	
	#------------------------------------------------------------
	import requests 
	from astropy.io.votable import parse_single_table 
	#------------------------------------------------------------
	if color and format == "fits":
		raise ValueError("color images are available only for jpg or png formats")
	if format not in ("jpg","png","fits"):
		raise ValueError("format must be one of jpg, png, fits")

	table	= getimages(ra, dec, size=size, filters=filters)
	url		= (	"https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
				"ra={ra}&dec={dec}&size={size}&format={format}").format(**locals())
	if output_size:
		url = url + "&output_size={}".format(output_size)
	# sort filters from red to blue
	flist = ["yzirg".find(x) for x in table['filter']]
	table = table[np.argsort(flist)]
	if color:
		if len(table) > 3:
			# pick 3 filters
			table = table[[0,len(table)//2,len(table)-1]]
		for i, param in enumerate(["red","green","blue"]):
			url = url + "&{}={}".format(param,table['filename'][i])
	else:
		urlbase = url + "&red="
		url = []
		for filename in table['filename']:
			url.append(urlbase+filename)
	return url
#------------------------------------------------------------
def getimages(ra,dec,size=240,filters="grizy"):
	"""
	Query ps1filenames.py service to get a list of images
	ra, dec = position in degrees
	size = image size in pixels (0.25 arcsec/pixel)
	filters = string with filters to include
	Returns a table with the results
	"""	
	service	= "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
	url		= ("{service}?ra={ra}&dec={dec}&size={size}&format=fits"
		   "&filters={filters}").format(**locals())
	table	= Table.read(url, format='ascii')
	return table
#------------------------------------------------------------
def downimage_routine(outim, name, ra, dec, size, output_size, filters, format, save_dir='.'):
	filt = filters
	param_geturl= dict(	ra			= ra,
						dec			= dec,
						size		= size,
						output_size	= None,
						filters		= filt,
						format		= "fits")
	url			= geturl(**param_geturl)
	fh		= fits.open(url[0])
	# newname	= 'Ref'+'-PS1-'+name+'-'+filt+'.fits'
	newname = outim
	fh.writeto(save_dir+'/'+newname, overwrite=True)
	pan, panhd	= fits.getdata(save_dir+'/'+newname, header=True)
	pan0		= np.nan_to_num(pan)
	fits.writeto(save_dir+'/'+newname, pan0, panhd, overwrite=True)
#------------------------------------------------------------
#	SAMPLE LINES FOR DOWNLOADING SINGLE IMAGE FROM PS1
'''
param_down	= dict(  name = 'AT2019dko',
                    ra = 181.4641667,
                    dec = 67.2569444,
                    size = 5000,
                    output_size = None,
                    filters = 'r',
                    format = 'fits',
                    save_dir='.')
try:
	query.downimage_routine(**param_down)
except:
	try:
		query.downimage_routine(**param_down)
	except:
		try:
			query.downimage_routine(**param_down)		
		except:
			query.downimage_routine(**param_down)

#------------------------------------------------------------
#	SAMPLE LINES FOR DOWNLOADING MULTIPLE IMAGES FROM PS1
#	intbl : 'ascii' TABLE
for i in range(len(intbl)):
	print('['+str(i+1)+'/'+str(len(intbl))+']')
	param_down= dict(	name		= intbl['name'][i],
						ra			= intbl['ra'][i],
						dec			= intbl['dec'][i],
						size		= 5000,
						output_size	= None,
						filters		= 'r',
						format		= "fits")
	try:
		downimage_routine(**param_down)
	except:
		try:
			downimage_routine(**param_down)
		except:
			downimage_routine(**param_down)
'''
#------------------------------------------------------------
def querybox(refcatname, obj, racent, decent, path_refcat, radius=0.5, refmagkey=''):
	'''
	reftbl = querybox(**param_query)
	'''
	#------------------------------------------------------------
	#	REF. CATALOG QUERY
	#------------------------------------------------------------
	refcatlist	= glob.glob(path_refcat+'/*.cat')
	#------------------------------------------------------------
	if refcatname	== 'PS1':
		if path_refcat+'/ps1-'+obj+'.cat' not in refcatlist:
			querytbl = phot.ps1_query(obj, racent, decent, path_refcat, radius=radius)
		else:
			querytbl = ascii.read(path_refcat+'/ps1-'+obj+'.cat')
		reftbl, refcat = phot.ps1_Tonry(querytbl, obj)
	#------------------------------------------------------------
	elif refcatname	== 'SDSS':
		if path_refcat+'/sdss-'+obj+'.cat' not in refcatlist:
			querytbl = phot.sdss_query(obj, racent, decent, path_refcat)
		else:
			querytbl = ascii.read(path_refcat+'/sdss-'+obj+'.cat')
		reftbl, refcat = phot.sdss_Blaton(querytbl, obj)
	#------------------------------------------------------------
	elif refcatname == 'APASS':
		if path_refcat+'/apass-'+obj+'.cat' not in refcatlist:
			querytbl = phot.apass_query(obj, racent, decent, path_refcat, radius=radius)
		else:
			querytbl = ascii.read(path_refcat+'/apass-'+obj+'.cat')
		reftbl, refcat = phot.apass_Blaton(querytbl, obj)
	#------------------------------------------------------------
	elif refcatname	== '2MASS':
		if path_refcat+'/2mass-'+obj+'.cat' not in refcatlist:
			querytbl        = phot.twomass_query(obj, racent, decent, path_refcat, band=refmagkey, radius=radius)
		else:
			querytbl = ascii.read(path_refcat+'/2mass-'+obj+'.cat')
		reftbl, refcat = querytbl, '2mass-'+obj+'.cat'
	return reftbl

