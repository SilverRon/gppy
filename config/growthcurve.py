#	GROWTH CURVE IN PHOTOMETRY FOR PYTHON 3.X
#	CREATED	2019.12.26	Gregory S.H. Paek
#============================================================
import os, glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
from astropy.io import ascii
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
from imsng import phot
import time
#============================================================
#	FUNCTION
#============================================================
def getccdinfo(obs, path_obs):
	'''
	GET CCD INFORMATION (GAIN, PIXEL SCALE, FOV)

	gain, pixscale, fov = getccdinfo(obs, path_obs)
	
	INPUT:
	obs = 'SAO'
	path_obs = '/home/sonic/Research/table'

	OUTPUT:
	gain, pixscale, fov
	'''
	obstbl = ascii.read(path_obs+'/obs.txt')
	indx_obs = np.where(obstbl['obs']==obs)
	gain = obstbl[indx_obs]['gain'][0]
	pixscale = obstbl[indx_obs]['pixelscale'][0]
	fov = obstbl[indx_obs]['fov'][0]
	return gain, pixscale, fov
#------------------------------------------------------------
def secom(inim, param_insex, dualmode=False):
	'''
	
	'''
	param_sex = dict(	CONF_NAME = 'default.sex',
						#------------------------------
						#	CATALOG
						#------------------------------
						CATALOG_NAME = 'test.cat',
						CATALOG_TYPE = 'ASCII_HEAD',
						PARAMETERS_NAME = 'default.param',
						#------------------------------
						#	EXTRACTION
						#------------------------------
						DETECT_TYPE = 'CCD',
						DETECT_MINAREA = '5',
						DETECT_MAXAREA = '0',
						DETECT_THRESH = '1.5',
						# ANALYSIS_THRESH = 'RELATIVE',
						ANALYSIS_THRESH = '1.5',						
						FILTER = 'Y',
						FILTER_NAME = 'default.conv',
						DEBLEND_NTHRESH = '64',
						DEBLEND_MINCONT = '0.0001',
						CLEAN = 'Y',
						CLEAN_PARAM = '1.0',
						MASK_TYPE = 'CORRECT',
						#------------------------------
						#	PHOTOMETRY
						#------------------------------
						# PHOT_APERTURES = '3',
						PHOT_AUTOPARAMS = '2.5,3.5',
						PHOT_PETROPARAMS = '2.0,3.5',
						SATUR_LEVEL  = '50000.0',
						SATUR_KEY = 'SQTURATE',
						MAG_ZEROPOINT = '0.0',
						MAG_GAMMA = '4.0',
						GAIN = '1.0',
						GAIN_KEY = 'GAIN',   
						PIXEL_SCALE = '1.0',
						#------------------------------
						#	STAR/GALAXY SEPARATION
						#------------------------------
						SEEING_FWHM = '3.0',
						STARNNW_NAME = 'default.nnw',
						#------------------------------
						#	BACKGROUND
						#------------------------------
						BACK_SIZE = '128',
						BACK_FILTERSIZE = '10',
						BACKPHOTO_TYPE = 'LOCAL',
						#------------------------------
						#	CHECK IMAGE
						#------------------------------
						CHECKIMAGE_TYPE = 'NONE',
						CHECKIMAGE_NAME = 'check.fits',
						#==============================
						#	MEMORY & MISCELLANEOUS
						#==============================
						MEMORY_OBJSTACK = '3000',
						MEMORY_PIXSTACK = '300000',
						MEMORY_BUFSIZE = '1024',
						VERBOSE_TYPE = 'NORMAL',
						HEADER_SUFFIX = '.head',
						WRITE_XML = 'N',
						XML_NAME = 'sex.xml')

	for key in param_insex.keys():
		param_sex[key] = param_insex[key]

	
	secom_normal = 'sex -c {} {} '.format(param_sex['CONF_NAME'], inim)
	secom_dual = 'sex -c {} {} '.format(param_sex['CONF_NAME'], inim)
	for key in param_sex.keys():
		if key != 'CONF_NAME':
			secom_normal += '-{} {} '.format(key, param_sex[key])

	# print(secom_normal)
	os.system(secom_normal)
#------------------------------------------------------------
def querybox(refcatname, obj, racent, decent, path_refcat, radius=0.5):
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
			querytbl = phot.apass_query(obj, racent, decent, path_refcat)
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
#============================================================
#	USER SETTING
#============================================================
path_base = './'
radius = 0.5								#	[DEGREE]
frac = 0.9									#	IMAGE FRACTION [%]
# refcatname = '2MASS'						#	REFERENCE CATALOG
refcatname = 'PS1'
inmagkey = 'MAG_APER'
inmagerkey = 'MAGERR_APER'
refmaglower, refmagupper = 14, 17			#	REF MAG RANGE [MAG]
refmagerupper = 0.1
inmagerupper = 0.1
# aperture = str(input('aperture:'))			#	DIAMETER
# tra, tdec = 44.54382542, -8.9577325		#	GRB190829A (MORE ACCURATE)
tra, tdec = 44.62963541, -8.951577778		#	GRB190829A (MORE ACCURATE)
phottype = 'normal'
apertures = np.arange(2, 50.0, 1.5)
aper_input = ''
for i in apertures:
	aper_input = aper_input+'{},'.format(i)
aper_input = aper_input[:-1]
#------------------------------------------------------------
#	PATH
#------------------------------------------------------------
path_refcat	= '/home/sonic/Research/cat/refcat'
path_obs = '/home/sonic/Research/table'
path_config = '/home/sonic/Research/yourpy/gppy/config'
#============================================================
#	PHOTOMETRY
#============================================================
# inim = '/data1/GRB190829A/ukirt/test2/Calib-UKIRT-GRB190829A-20190830-123028-K-10.fits'
# inim = 'Calib-DOAO-NGC5350-20190503-145856-R-60.fits'
# inim = 'Calib-MAO_SNUCAM-NGC5350-20190712-171053-R-120.fits'
# inim = 'Calib-KMTNET-S190425z_FIELD-20190000-000000-R-0.fits'
# inim = 'Calib-SAO_STX16803-NGC5350-20190521-115820-R-120.fits'
# inim = 'Calib-LOAO-NGC5350-20190503-062427-R-60.fits'
# inim = 'Calib-SAO_STX16803-NGC5350-20190501-190240-R-180.fits'

head = os.path.splitext(inim)[0]

hdr = fits.getheader(inim)
xcent, ycent= hdr['NAXIS1']/2., hdr['NAXIS2']/2.
#------------------------------------------------------------
#	RA, Dec CENTER FOR QUERYING
#------------------------------------------------------------
try:
	w = WCS(inim)
	racent, decent = w.all_pix2world(xcent, ycent, 1)
	racent, decent = racent.item(), decent.item()
except:
	print('BAD WCS INFORMATION?')
	racent, decent = hdr['CRVAL1'], hdr['CRVAL2']
#------------------------------------------------------------
#	DATE-OBS, JD
#------------------------------------------------------------
try:
	date_obs = hdr['date-obs']
	jd = round(Time(date_obs, format='isot', scale='utc').jd, 3)
except:
	date_obs = None
	jd = None
#------------------------------------------------------------
#	NAME INFO
#------------------------------------------------------------
part = inim.split('-')
obs = part[1]
obj = part[2]
refmagkey = part[5]
refmagerkey = refmagkey+'err'
#------------------------------------------------------------
gain, pixscale, fov = getccdinfo(obs, path_obs)
#------------------------------------------------------------
cat = head+'.cat'
seg = head+'.seg.fits'
bkg = head+'.bkg.fits'
sub = head+'.sub.fits'
psf = head+'.psf'
aper = head+'.aper.fits'

param = path_config+'/growthcurve.param'
conv = path_config+'/growthcurve.conv'
nnw = path_config+'/growthcurve.nnw'
conf = path_config+'/growthcurve.sex'
#------------------------------------------------------------
#	SOURCE EXTRACTOR CONFIGURATION
#------------------------------------------------------------
param_insex = dict(	#------------------------------
					#	CATALOG
					#------------------------------
					CATALOG_NAME = cat,
					#------------------------------
					#	CONFIG FILES
					#------------------------------
					CONF_NAME = conf,
					PARAMETERS_NAME = param,
					FILTER_NAME = conv,    
					STARNNW_NAME = nnw,
					#------------------------------
					#	EXTRACTION
					#------------------------------			
					# PSF_NAME = psf,
					DETECT_MINAREA = '5',
					DETECT_THRESH = '3.0',
					DEBLEND_NTHRESH = '64',
					DEBLEND_MINCONT = '0.0001',
					#------------------------------
					#	PHOTOMETRY
					#------------------------------
					PHOT_APERTURES = aper_input,
					SATUR_LEVEL  = '65000.0',
					# MAG_ZEROPOINT = '0.0',
					GAIN = str(gain),
					PIXEL_SCALE = str(pixscale),
					#------------------------------
					#	STAR/GALAXY SEPARATION
					#------------------------------
					SEEING_FWHM = '3.0',
					#------------------------------
					#	BACKGROUND
					#------------------------------
					BACK_SIZE = '128',
					BACK_FILTERSIZE = '10',
					BACKPHOTO_TYPE = 'LOCAL',
					#------------------------------
					#	CHECK IMAGE
					#------------------------------
					CHECKIMAGE_TYPE = 'SEGMENTATION,APERTURES,BACKGROUND,-BACKGROUND',
					CHECKIMAGE_NAME = '{},{},{},{}'.format(seg, aper, bkg, sub),
					)

# peeing, seeing	= phot.psfex(inim, pixscale)
secom(inim, param_insex)
setbl = ascii.read(cat)

reftbl = querybox(refcatname, obj, racent, decent, path_refcat, radius=3.0)
#------------------------------------------------------------
#	CENTER POS. & DIST CUT
#------------------------------------------------------------
deldist = phot.sqsum((xcent-setbl['X_IMAGE']), (ycent-setbl['Y_IMAGE']))
indx_dist = np.where(deldist < np.sqrt(frac)*(xcent+ycent)/2.)
intbl = setbl[indx_dist]
intbl.write(cat, format='ascii', overwrite=True)
#------------------------------------------------------------
#	MATCHING
#------------------------------------------------------------
param_match = dict(	intbl=intbl, reftbl=reftbl,
					inra=intbl['ALPHA_J2000'], indec=intbl['DELTA_J2000'],
					refra=reftbl['ra'], refdec=reftbl['dec'],
					sep=2.0)
mtbl = phot.matching(**param_match)
#------------------------------------------------------------
#	ZEROPOINT CALCULATION
#------------------------------------------------------------
param_st4zp	= dict(	intbl=mtbl,
					inmagerkey=inmagkey,
					refmagkey=refmagkey,
					refmagerkey=refmagerkey,
					refmaglower=refmaglower,
					refmagupper=refmagupper,
					refmagerupper=refmagerupper,
					inmagerupper=inmagerupper)
param_zpcal	= dict(	intbl=phot.star4zp(**param_st4zp),
					inmagkey=inmagkey, inmagerkey=inmagerkey,
					refmagkey=refmagkey, refmagerkey=refmagerkey,
					sigma=2.0)
zp, zper, otbl, xtbl = phot.zpcal(**param_zpcal)
#------------------------------------------------------------
#	
#------------------------------------------------------------
alltbl = vstack([otbl, xtbl])
alltbl = alltbl[alltbl['NUMBER'].argsort()]

def extractonerow(alltbl, number, apertures):
	onerow = alltbl[alltbl['NUMBER']==number]
	# onerow = alltbl[0]
	mag, mager = [], []
	flux, fluxer = [], []
	aper = apertures
	snr = []
	for i in range(len(aper)):
		if i == 0:
			tail = ''
		else:
			tail = '_{}'.format(i)
		mag.append(onerow['MAG_APER'+tail])
		mager.append(onerow['MAGERR_APER'+tail])
		flux.append(onerow['FLUX_APER'+tail])
		fluxer.append(onerow['FLUXERR_APER'+tail])
		snr.append((onerow['FLUX_APER'+tail]/onerow['FLUXERR_APER'+tail]).item())
	magdif = []
	magdifer = []
	fluxdif = []
	fluxdifer = []

	for i in range(len(mag)):
		if i == 0:
			mag0 = mag[0]
			mager0 = mager[0]
			flux0 = flux[0]
			fluxer0 = fluxer[0]
			pass
		else:
			magdif.append(mag[i]-mag0)
			magdifer.append(np.sqrt( mager[i]**2 + mager0**2 ))
			mag0 = mag[i]
			mager0 = mager[i]
			flux0 = flux[i]
			fluxer0 = fluxer[i]
	plt.scatter(aper[1:], snr[1:])
	plt.plot(aper[1:], snr[1:], color='grey', alpha=0.5)
	x, y = aper[1:], snr[1:]
	optaper = x[y == np.max(y)]
	plt.axvline(x=optaper)

	# plt.errorbar(aper[1:], magdif, yerr=magdifer,
				# color='grey', alpha=0.5,
				# marker='o', linestyle='')
	
	return magdif, magdifer, optaper

def meandelmagdif(alltbl, apertures):
	aper = apertures
	y = []
	yer = []
	for i in range(len(aper)):
		if i == 0:
			tail = ''
			mag0 = alltbl['MAG_APER{}'.format(tail)]
		else:
			tail = '_{}'.format(i)
			mag1 = alltbl['MAG_APER{}'.format(tail)]
			y.append(np.mean(mag1-mag0))
			yer.append(np.std(mag1-mag0))
			mag0 = mag1
	plt.errorbar(aper[1:], y, yerr=yer,
				color='gold', alpha=0.5,
				marker='*', markersize=15, linestyle='')
	return y, yer






plt.close('all')
optapers = []
for number in alltbl['NUMBER']:
	magdif, magdifer, optaper = extractonerow(alltbl, number, apertures)
	optapers.append(optaper.item())

# y, yer = meandelmagdif(alltbl)

# plt.plot(aper[1:], yer, color='dodgerblue')
# plt.axvline(aper[1:][yer == np.min(yer)], color='tomato', linestyle='--')

optaper = round(np.mean(optapers), 2)
plt.axvline(x=optaper, color='tomato', label='Opt.Aper {}'.format(optaper))
print(optaper)

down, up = plt.ylim()

plt.xlim([0, 50])
plt.ylim([0, up])
plt.legend(fontsize=20, loc='best', framealpha=1, markerscale=2)
plt.tick_params(which='both', direction='in')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
# plt.gca().invert_yaxis()
plt.title(inim, fontsize=15)
plt.xlabel(r'Diameter [pixel]', fontsize=20)
plt.ylabel(r'SNR', fontsize=20)
plt.grid(color='dimgrey', linestyle='--', linewidth=1, alpha=0.5)
plt.tight_layout()
plt.minorticks_on()



plt.savefig(inim[:-5]+'.gcurve.png', overwrite=True)