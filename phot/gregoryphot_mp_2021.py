#	PHOTOMETRY CODE FOR PYTHON 3.X
#	CREATED	2020.12.10	Gregory S.H. Paek
#============================================================
import os, glob, sys, subprocess
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
from astropy.io import ascii
from astropy.io import fits
from astropy.time import Time
# from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
from datetime import date
import time
import multiprocessing
sys.path.append('..')
sys.path.append('/home/gecko/gppy')
from phot import gpphot
from util import query
from util import tool
from phot import gcurve
from preprocess import calib
starttime = time.time()
#============================================================
#	FUNCTION
#============================================================
def file2dict(path_infile):
	out_dict = dict()
	f = open(path_infile)
	for line in f:
		key, val = line.split()
		out_dict[key] = val
	return out_dict
#------------------------------------------------------------
def phot_routine(inim):
	#------------------------------------------------------------
	#	INFO. from file name
	#------------------------------------------------------------
	hdr = fits.getheader(inim)
	part = os.path.splitext(os.path.basename(inim))[0].split('-')
	head = os.path.splitext(inim)[0]

	obs = part[1]
	# obj = part[2]
	# refmagkey = part[5]
	obj = hdr['OBJECT']
	refmagkey = hdr['FILTER']
	
	refmagerkey = refmagkey+'err'
	#------------------------------------------------------------
	print(inim, obs, obj, refmagkey, refmagerkey)
	if obs == 'GALEX':
		obsdict = dict(
		gain=float(hdr[f"{refmagkey[0:1]}SXGAIN"])*u.electron/u.second,
		pixelscale=1.5*u.arcsec/u.pixel,
		fov=72,
		)
	else:
		obsdict = tool.getccdinfo(obs, path_obs)
	gain = obsdict['gain']
	pixscale = obsdict['pixelscale']
	fov = obsdict['fov']
	# rdnoise = obsdict['readoutnoise']
	if obs == 'SAO_C361K':
		pixscale = pixscale*hdr['XBINNING']
	#------------------------------------------------------------
	#	OUTPUT NAMES
	cat = head+'.cat'
	cat_gc = head+'.gcurve.cat'
	seg = head+'.seg.fits'
	bkg = head+'.bkg.fits'
	sub = head+'.sub.fits'
	psf = head+'.psf'
	aper = head+'.aper.fits'
	#	GROWTH CURVE NAMES
	param_gc = path_config+'/growthcurve.param'
	conv_gc = path_config+'/growthcurve.conv'
	nnw_gc = path_config+'/growthcurve.nnw'
	conf_gc = path_config+'/growthcurve.sex'
	#	PHOTOMETRY NAMES
	param = path_config+'/gregoryphot.param'
	conv = path_config+'/gregoryphot.conv'
	nnw = path_config+'/gregoryphot.nnw'
	conf = path_config+'/gregoryphot.sex'
	#------------------------------------------------------------
	#	Aperture determine (diameter for SE input)
	# aper_lower = 2*0.5*seeing_assume/pixscale
	# aper_upper = 2*1.5*seeing_assume/pixscale
	aper_lower = 1.0 * u.arcsecond/pixscale
	aper_upper = 10.0 * u.arcsecond/pixscale
	# apertures = np.linspace(aper_lower.value, aper_upper.value)
	apertures = np.linspace(aper_lower, aper_upper, 32)
	aper_input = ''
	for i in apertures.value: aper_input = aper_input+'{},'.format(i)
	aper_input = aper_input[:-1]

	print('-'*60)
	# print('[{}/{}] {}'.format(j+1, len(imlist), inim))
	print(inim)
	print('{}\t{} in {}-band'.format(obs, obj, refmagkey))
	print('-'*60)
	#------------------------------------------------------------
	#	INFO. from file header
	#------------------------------------------------------------
	hdul = fits.open(inim)
	hdr = hdul[0].header
	# hdr = fits.getheader(inim)
	#	RA, Dec center for reference catalog query
	xcent, ycent= hdr['NAXIS1']/2., hdr['NAXIS2']/2.
	w = WCS(inim)
	racent, decent = w.all_pix2world(xcent, ycent, 1)
	racent, decent = racent.item(), decent.item()
	# print('BAD WCS INFORMATION?')
	# racent, decent = hdr['CRVAL1'], hdr['CRVAL2']
	#------------------------------------------------------------
	#	DATE-OBS, JD
	#------------------------------------------------------------
	try:
		date_obs = hdr['date-obs']
		t = Time(date_obs, format='isot', scale='utc')
		jd = round(t.jd, 3)
		mjd = round(t.mjd, 3)
	except:
		print('SKIP : No DATE-OBS information on header')
		date_obs = None
		jd = None
		mjd = None
	#------------------------------------------------------------
	#	SOURCE EXTRACTOR CONFIGURATION FOR GROTH CURVE
	#------------------------------------------------------------
	param_gcurve = dict(#------------------------------
						#	CATALOG
						#------------------------------
						CATALOG_NAME = cat_gc,
						#------------------------------
						#	CONFIG FILES
						#------------------------------
						CONF_NAME = conf_gc,
						PARAMETERS_NAME = param_gc,
						FILTER_NAME = conv_gc,    
						STARNNW_NAME = nnw_gc,
						#------------------------------
						#	EXTRACTION
						#------------------------------			
						# PSF_NAME = psf,
						DETECT_MINAREA = DETECT_MINAREA,
						DETECT_THRESH = DETECT_THRESH,
						DEBLEND_NTHRESH = DEBLEND_NTHRESH,
						DEBLEND_MINCONT = DEBLEND_MINCONT,
						#------------------------------
						#	PHOTOMETRY
						#------------------------------
						PHOT_APERTURES = aper_input,
						SATUR_LEVEL  = '65000.0',
						# MAG_ZEROPOINT = '0.0',
						GAIN = str(gain.value),
						PIXEL_SCALE = str(pixscale.value),
						#------------------------------
						#	STAR/GALAXY SEPARATION
						#------------------------------
						SEEING_FWHM = str(seeing_assume.value),
						#------------------------------
						#	BACKGROUND
						#------------------------------
						BACK_SIZE = BACK_SIZE,
						BACK_FILTERSIZE = BACK_FILTERSIZE,
						BACKPHOTO_TYPE = BACKPHOTO_TYPE,
						#------------------------------
						#	CHECK IMAGE
						#------------------------------
						# CHECKIMAGE_TYPE = 'SEGMENTATION,APERTURES,BACKGROUND,-BACKGROUND',
						# CHECKIMAGE_NAME = '{},{},{},{}'.format(seg, aper, bkg, sub),
						)
	print('1. GROWTH CURVE')
	os.system(gpphot.sexcom(inim, param_gcurve))
	setbl = ascii.read(cat_gc)
	# setbl.write(cat, format='ascii', overwrite=True)
	reftbl = query.querybox(refcatname, obj, racent, decent, path_refcat, radius=refqueryradius.value, refmagkey=refmagkey)
	reftbl['NUMBER'] = np.arange(0, len(reftbl), 1)
	#	CENTER POS. & DIST CUT
	deldist = tool.sqsum((xcent-setbl['X_IMAGE']), (ycent-setbl['Y_IMAGE']))
	indx_dist = np.where(deldist < frac*(xcent+ycent)/2.)
	intbl = setbl[indx_dist]
	#	MATCHING
	print(reftbl)
	param_match = dict(	
						intbl=intbl, reftbl=reftbl,
						inra=intbl['ALPHA_J2000'], indec=intbl['DELTA_J2000'],
						refra=reftbl['ra'], refdec=reftbl['dec'],
						# sep=float(seeing_input),
						sep=5.0,
						)
	mtbl = gpphot.matching(**param_match)
	#	Exception for bad astrometry --> insufficient matced stars
	if len(mtbl) < 3:
		print('Suspect bad astrometry')
		if os.path.dirname(inim) == '':
			path_im = '.'
		else:
			path_im = os.path.dirname(inim)
		calib.astrometry('{}/{}'.format(path_im, inim), pixscale.value, ra=None, dec=None, fov=1, cpulimit=60)
		outim = '{}/a{}'.format(path_im, os.path.basename(inim))
		mvcom = 'mv {} {}'.format(outim, inim)
		print(mvcom)
		os.system(mvcom)
	else:
		pass
	#	ZEROPOINT CALCULATION
	param_st4zp = dict(	intbl=mtbl,
						inmagerkey='MAG_AUTO',
						refmagkey=refmagkey,
						refmagerkey=refmagerkey,
						refmaglower=refmaglower,
						refmagupper=refmagupper,
						refmagerupper=refmagerupper,
						inmagerupper=inmagerupper,
						flagcut=flagcut,
						)
	st4tbl = gpphot.star4zp(**param_st4zp)
	param_zpcal = dict(	intbl=st4tbl,
						inmagkey='MAG_AUTO', inmagerkey='MAGERR_AUTO',
						refmagkey=refmagkey, refmagerkey=refmagerkey,
						sigma=2.0,
						# method='weightedmean',
						)
	'''
	try:
		_, _, otbl, xtbl = gpphot.zpcal(**param_zpcal)
	except Exception as e:
		print(inim, e)
		failist.append(inim)
	'''
	_, _, otbl, xtbl = gpphot.zpcal(**param_zpcal)

	alltbl = vstack([otbl, xtbl])
	alltbl = alltbl[alltbl['NUMBER'].argsort()]
	seeing = np.median(alltbl['FWHM_WORLD'])*3600*u.arcsecond
	peeing = np.median(alltbl['FWHM_IMAGE'])*u.arcsecond
	seeing_input = str(seeing.value)

	optaper, optapers = gcurve.gcurveplot(inim, alltbl, apertures=apertures, pixelscale=pixscale)
	#------------------------------------------------------------
	#	APERTURE SETTING
	#------------------------------------------------------------
	# for inmagkey, inmagerkey in zip():
	inmagkeys = [
				'MAG_AUTO',
				'MAG_APER',
				'MAG_APER_1',
				'MAG_APER_2',
				'MAG_APER_3',
				'MAG_APER_4',
				'MAG_APER_5',
				]
	inmagerkeys = []
	for key in inmagkeys: inmagerkeys.append(key.replace('MAG_', 'MAGERR_'))
	aperkeys = []
	for key in inmagkeys: aperkeys.append(key.replace('MAG_', ''))
	aperlist = [0, optaper, 2*0.6731*peeing.value, peeing.value*2, peeing.value*3, (3*u.arcsecond/pixscale).value, (5*u.arcsecond/pixscale).value]	
	aperdiscription = ['MAG_AUTO DIAMETER [pix]', 'BEST APERTURE DIAMETER in SNR curve [pix]', 'BEST GAUSSIAN APERTURE DIAMETER [pix]', '2*SEEING APERTURE DIAMETER [pix]', '3*SEEING APERTURE DIAMETER [pix]', """FIXED 3" APERTURE DIAMETER [pix]""", """FIXED 5" APERTURE DIAMETER [pix]""",]
	# aperlist = [0, optaper, 2*0.6731*peeing.value, peeing.value*2, peeing.value*3, peeing.value*4, peeing.value*5]
	# aperdiscription = ['MAG_AUTO DIAMETER [pix]', 'BEST APERTURE DIAMETER in SNR curve [pix]', 'BEST GAUSSIAN APERTURE DIAMETER [pix]', '2*SEEING APERTURE DIAMETER [pix]', '3*SEEING APERTURE DIAMETER [pix]', '4*SEEING APERTURE DIAMETER [pix]', '5*SEEING APERTURE DIAMETER [pix]',]
	#------------------------------------------------------------
	#	SOURCE EXTRACTOR CONFIGURATION FOR PHOTOMETRY
	#------------------------------------------------------------
	# optaper = 50
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
						DETECT_MINAREA = DETECT_MINAREA,
						DETECT_THRESH = DETECT_THRESH,
						DEBLEND_NTHRESH = DEBLEND_NTHRESH,
						DEBLEND_MINCONT = DEBLEND_MINCONT,
						#------------------------------
						#	PHOTOMETRY
						#------------------------------
						#	DIAMETER
						#	OPT.APER, (SEEING x2), x3, x4, x5
						#	MAG_APER	OPT.APER
						#	MAG_APER_1	OPT.GAUSSIAN.APER
						#	MAG_APER_2	SEEINGx2
						#	...
						PHOT_APERTURES = '{},{},{},{},{},{}'.format(aperlist[1], aperlist[2], aperlist[3], aperlist[4], aperlist[5], aperlist[6], ),
						# PHOT_APERTURES = '{},{},{},{},{},{}'.format(optaper, 2*0.6731*peeing.value, peeing.value*2, peeing.value*3, peeing.value*4, peeing.value*5),
						# PHOT_APERTURES = ','.join(aperlist[1:]),
						SATUR_LEVEL  = '65000.0',
						# MAG_ZEROPOINT = '0.0',
						GAIN = str(gain.value),
						PIXEL_SCALE = str(pixscale.value),
						#------------------------------
						#	STAR/GALAXY SEPARATION
						#------------------------------
						SEEING_FWHM = str(seeing_assume.value),
						#------------------------------
						#	BACKGROUND
						#------------------------------
						BACK_SIZE = BACK_SIZE,
						BACK_FILTERSIZE = BACK_FILTERSIZE,
						BACKPHOTO_TYPE = BACKPHOTO_TYPE,
						#------------------------------
						#	CHECK IMAGE
						#------------------------------
						# CHECKIMAGE_TYPE = 'SEGMENTATION,APERTURES,BACKGROUND,-BACKGROUND',
						# CHECKIMAGE_NAME = '{},{},{},{}'.format(seg, aper, bkg, sub),
						)
	if check == True:
		param_insex['CHECKIMAGE_TYPE'] = 'SEGMENTATION,APERTURES,BACKGROUND,-BACKGROUND'
		param_insex['CHECKIMAGE_NAME'] = '{},{},{},{}'.format(seg, aper, bkg, sub)
	else:
		pass
	print('2. SOURCE EXTRACTOR')
	com = gpphot.sexcom(inim, param_insex)
	sexout = subprocess.getoutput(com)
	line = [s for s in sexout.split('\n') if 'RMS' in s]
	skymed, skysig = float(line[0].split('Background:')[1].split('RMS:')[0]), float(line[0].split('RMS:')[1].split('/')[0])
	os.system('rm {} {} {} {}'.format(seg, aper, bkg, sub))
	setbl = ascii.read(cat)
	setbl['obs'] = obs
	setbl['obj'] = obj
	setbl['filter'] = refmagkey
	setbl['date-obs'] = date_obs
	setbl['jd'] = jd
	setbl['mjd'] = mjd
	#------------------------------------------------------------
	#	CENTER POS. & DIST CUT
	#------------------------------------------------------------
	deldist = tool.sqsum((xcent-setbl['X_IMAGE']), (ycent-setbl['Y_IMAGE']))
	# indx_dist = np.where(deldist < np.sqrt(frac)*(xcent+ycent)/2.)
	indx_dist = np.where(deldist < frac*(xcent+ycent)/2.)
	# intbl = setbl
	# intbl.write(cat, format='ascii', overwrite=True)
	frctbl = setbl[indx_dist]
	#------------------------------------------------------------
	#	MATCHING
	#------------------------------------------------------------
	param_match = dict(
						intbl=frctbl, reftbl=reftbl,
						inra=frctbl['ALPHA_J2000'], indec=frctbl['DELTA_J2000'],
						refra=reftbl['ra'], refdec=reftbl['dec'],
						sep=float(seeing_input),
					)
						# sep=3.0)
	print('3. MATCHING')
	mtbl = gpphot.matching(**param_match)
	# mtbl.write(cat, format='ascii', overwrite=True)
	intbl = setbl
	intbl['FWHM_IMAGE'] = intbl['FWHM_IMAGE'] 
	intbl['FWHM_WORLD'] = intbl['FWHM_WORLD'].to(u.arcsecond)
	ellipticity = np.median(intbl['ELLIPTICITY'])
	elongation = np.median(intbl['ELONGATION'])
	#------------------------------------------------------------
	#	Header
	#------------------------------------------------------------
	tool.puthdr(inim, 'AUTHOR', 'Gregory S.H. Paek',		hdrcomment='PHOTOMETRY AUTHOR')
	tool.puthdr(inim, 'JD', jd, hdrcomment='Julian date')
	tool.puthdr(inim, 'MJD', mjd, hdrcomment='Modified Julian date')
	tool.puthdr(inim, 'PHOTIME',date.today().isoformat(),	hdrcomment='PHTOMETRY TIME [KR]')
	tool.puthdr(inim, 'SEEING',	round(seeing.value, 3),	hdrcomment='SEEING [arcsec]')
	tool.puthdr(inim, 'PEEING',	round(peeing.value, 3),	hdrcomment='SEEING [pixel]')
	tool.puthdr(inim, 'ELLIP', round(ellipticity, 3), hdrcomment='ELLIPTICITY 1-B/A [0-1]')
	tool.puthdr(inim, 'ELONG', round(elongation, 3), hdrcomment='ELONGATION A/B [1-]')
	tool.puthdr(inim, 'APERPIX', optaper,					hdrcomment='APERTURE DIAMETER [pixel]')
	tool.puthdr(inim, 'APER', round(optaper*pixscale.value, 3), 	hdrcomment='APERTURE DIAMETER [arcsec]')
	tool.puthdr(inim, 'NAPER', round(optaper/peeing.value, 3), 	hdrcomment='N = APERTURE/PEEING')
	tool.puthdr(inim, 'SKYSIG',	round(skysig, 3),			hdrcomment='SKY SIGMA VALUE')
	tool.puthdr(inim, 'SKYVAL',	round(skymed, 3),			hdrcomment='SKY MEDIAN VALUE')
	tool.puthdr(inim, 'REFCAT', refcatname,					hdrcomment='REFERENCE CATALOG NAME')
	tool.puthdr(inim, 'MAGLOW', refmaglower,				hdrcomment='REF MAG RANGE, BRIGHT LIMIT')
	tool.puthdr(inim, 'MAGUP', refmagupper,					hdrcomment='REF MAG RANGE, DIM LIMIT')
	#------------------------------------------------------------
	#	ZEROPOINT CALCULATION
	#------------------------------------------------------------
	onetbl = Table()
	onetbl['obs'] = [obs]
	onetbl['field'] = obj
	onetbl['filter'] = refmagkey
	onetbl['racent'] = racent *u.degree
	onetbl['decent'] = decent *u.degree
	onetbl['date-obs'] = date_obs
	onetbl['jd'] = jd
	onetbl['mjd'] = mjd
	onetbl['skyval'] = skymed
	onetbl['skysig'] = skysig
	onetbl['seeing'] = seeing
	onetbl['peeing'] = peeing
	onetbl['ellipticity'] = ellipticity
	onetbl['elongation'] = elongation
	onetbl['image'] = inim
	# k = 0
	for k in range(len(inmagkeys)):
		aper = aperlist[k]
		aperkey = aperkeys[k]
		aperdisc = aperdiscription[k]
		inmagkey = inmagkeys[k]
		inmagerkey = inmagerkeys[k]
		print('4. ZERO POINT CALCULATION for {} [{}/{}]'.format(inmagkey, k+1, len(inmagkeys)))
		try:
			param_st4zp	= dict(	
								intbl=mtbl,
								inmagerkey=inmagerkey,
								refmagkey=refmagkey,
								refmagerkey=refmagerkey,
								refmaglower=refmaglower,
								refmagupper=refmagupper,
								refmagerupper=refmagerupper,
								inmagerupper=inmagerupper,
								verbose=True,
								plot=True,
								flagcut=flagcut,
								plotout='{}.{}.star.png'.format(inim[:-5], inmagkey)
								# plotout=inim[:-5]+'.star.png'
								)
			st4tbl = gpphot.star4zp(**param_st4zp)
			param_zpcal	= dict(	
								intbl=st4tbl,
								inmagkey=inmagkey, inmagerkey=inmagerkey,
								refmagkey=refmagkey, refmagerkey=refmagerkey,
								sigma=2.0,
								# method='weightedmean',
								)
			zp, zper, otbl, xtbl = gpphot.zpcal(**param_zpcal)
			#------------------------------------------------------------
			#	ZEROPOINT PLOT
			#------------------------------------------------------------
			param_zpplot = dict(
								outname='{}.zpcal.png'.format(head, inmagkey),
								# alltbl=mtbl,
								otbl=otbl, xtbl=xtbl,
								inmagkey=inmagkey, inmagerkey=inmagerkey,
								refmagkey=refmagkey, refmagerkey=refmagerkey,
								refmaglower=refmaglower, refmagupper=refmagupper,
								zp=zp, zper=zper
								)
			param_plot	= dict(
								inim		= inim,
								numb_list	= otbl['NUMBER'],
								xim_list	= otbl['X_IMAGE'],
								yim_list	= otbl['Y_IMAGE'],
								add			= True,
								numb_addlist= xtbl['NUMBER'],
								xim_addlist	= xtbl['X_IMAGE'],
								yim_addlist	= xtbl['Y_IMAGE']
								)
			gpphot.zpplot(**param_zpplot)
			if inmagkey == 'MAG_APER':
				print('4-1. PLOT')
				try:
					gpphot.plotshow(**param_plot)
					plt.close('all')
				except:
					print('FAIL TO DRAW ZEROPOINT GRAPH')
					pass
			else:
				pass
			#------------------------------------------------------------
			#	
			#------------------------------------------------------------
			# aper = 2*peeing
			if inmagkey == 'MAG_AUTO':
				ul_3sig = 0
				ul_5sig = 0
			else:
				ul_3sig = gpphot.limitmag(3, zp, aper, skysig)
				ul_5sig = gpphot.limitmag(5, zp, aper, skysig)
			#------------------------------------------------------------
			if inmagkey == 'MAG_AUTO':
				pass
			else:
				onetbl[inmagkey.replace('MAG_', '').replace('APER', 'aperture')] = aper
			onetbl[inmagkey.replace('MAG', 'STDNUMB').lower()] = len(otbl)
			onetbl[inmagkey.replace('MAG', 'ZP').lower()] = round(zp, 3) *u.ABmag
			onetbl[inmagkey.replace('MAG', 'ZPER').lower()] = round(zper, 3) *u.ABmag
			onetbl[inmagkey.replace('MAG', 'UL3SIG').lower()] = round(ul_3sig, 3) *u.ABmag
			onetbl[inmagkey.replace('MAG', 'UL5SIG').lower()] = round(ul_5sig, 3) *u.ABmag

			#------------------------------------------------------------
			#	ADD HEADER INFO
			#------------------------------------------------------------
			print('4-2. CHANGE HEADER')
			tool.puthdr(inim, inmagkey.replace('MAG_', ''), round(aper, 3), hdrcomment='{}'.format(aperdisc))
			tool.puthdr(inim, 'STDNUM_{}'.format(k), len(otbl), hdrcomment='# OF STD STARS for {}'.format(inmagkey))
			tool.puthdr(inim, 'ZP_{}'.format(k), round(zp, 3), hdrcomment='ZERO POINT for {}'.format(inmagkey))
			tool.puthdr(inim, 'ZPER_{}'.format(k), round(zper, 3), hdrcomment='ZERO POINT ERROR for {}'.format(inmagkey))
			tool.puthdr(inim, 'UL3_{}'.format(k), round(ul_3sig, 3), hdrcomment='3 sigma limit mag for {}'.format(inmagkey))
			tool.puthdr(inim, 'UL5_{}'.format(k), round(ul_5sig, 5), hdrcomment='5 sigma limit mag for {}'.format(inmagkey))
			#------------------------------------------------------------
			#	NORMAL PHOTOMETRY
			#------------------------------------------------------------
			print('4-3. PHOTOMETRY TABLE')
			intbl[inmagkey.lower()] = zp + intbl[inmagkey]
			intbl[inmagerkey.lower()] = tool.sqsum(zper, intbl[inmagerkey])
		except Exception as e:
			# print('Maybe the number of standard stars is insufficient.')
			print(inim, e)
			failist.append(inim)

	if transient == True:
		for tname, tra, tdec in zip(targtbl['name'], targtbl['ra'], targtbl['dec']):
			indx_targ, sep, _ = gpphot.targetfind(tra, tdec, intbl['ALPHA_J2000'], intbl['DELTA_J2000'], sep=seeing.value)
			if sep.to(u.arcsecond) < seeing:
				radeg, dedeg = intbl['ALPHA_J2000'][indx_targ].item() *u.degree, intbl['DELTA_J2000'][indx_targ].item() *u.degree
				# trtbl = intbl[indx_targ]
				trtbl = Table()
				for key in onetbl.keys(): trtbl[key] = onetbl[key]
				# trtbl = onetbl
				trtbl['name'] = tname
				trtbl['ra'] = radeg
				trtbl['dec'] = dedeg
				for key in ['NUMBER','X_IMAGE','Y_IMAGE','BACKGROUND','THRESHOLD','FLAGS','ELONGATION','ELLIPTICITY','CLASS_STAR','FWHM_IMAGE','FWHM_WORLD',]:
					trtbl[key] = intbl[key][indx_targ]
			else:
				trtbl = onetbl
			trtblist.append(trtbl)
	intbl.write('{}.phot.cat'.format(head), format='ascii.tab', overwrite=True)
	onetbl.write('{}.phot.summary.dat'.format(head), format='ascii.tab', overwrite=True)
	tblist.append(onetbl)
	plt.close('all')
#============================================================
#	USER SETTING
#============================================================
#	PATH
#------------------------------------------------------------
try:
	# obs = (sys.argv[1]).upper()
	path_base = sys.argv[1]
except:
	path_base = '.'
# print(path_base)
path_refcat	= '/data3/paek/factory/refcat'
# path_obs = '/home/paek/table/obs.txt'
path_obs = '/home/paek/table/obs.dat'
path_config = '/home/paek/config'
# path_target = '/home/sonic/Research/gppy/table/transient.dat'
path_target = './transient.dat'
path_gphot = f'{path_base}/gphot.config'
path_default_gphot = '/home/paek/config/gphot.config'
#------------------------------------------------------------
print(path_gphot)
if os.path.exists(path_gphot) == True:
	gphot_dict = file2dict(path_gphot)
else:
	gphot_dict = file2dict(path_default_gphot)
	print('There is no gregoryphot configuration. Use default.')
# gphot_dict = file2dict(path_gphot)
if len(glob.glob(path_target)) == 0:
	transient = False
	print('No transient catalog for photometry')
else:
	targtbl = ascii.read(path_target)
	transient = True
#------------------------------------------------------------
imkey = gphot_dict['imkey']
refqueryradius = float(gphot_dict['refqueryradius'])*u.degree
frac = float(gphot_dict['photfraction'])
refcatname = gphot_dict['refcatname']
refmaglower = float(gphot_dict['refmaglower'])
refmagupper = float(gphot_dict['refmagupper'])
refmagerupper = float(gphot_dict['refmagerupper'])
inmagerupper = float(gphot_dict['inmagerupper'])
flagcut = int(gphot_dict['flagcut'])
check = bool(gphot_dict['check'])

DETECT_MINAREA = gphot_dict['DETECT_MINAREA']
DETECT_THRESH = gphot_dict['DETECT_THRESH']
DEBLEND_NTHRESH = gphot_dict['DEBLEND_NTHRESH']
DEBLEND_MINCONT = gphot_dict['DEBLEND_MINCONT']
BACK_SIZE = gphot_dict['BACK_SIZE']
BACK_FILTERSIZE = gphot_dict['BACK_FILTERSIZE']
BACKPHOTO_TYPE = gphot_dict['BACKPHOTO_TYPE']
#------------------------------------------------------------
seeing_assume = 3.0 * u.arcsecond
#------------------------------------------------------------
# inmagkey = 'MAG_APER'
# inmagerkey = 'MAGERR_APER'
# inmagkey = 'MAG_AUTO'
# inmagerkey = 'MAGERR_AUTO'
#------------------------------------------------------------
imlist = sorted(glob.glob(imkey))
# ncore = 8
# ncore = 4
ncore = int(sys.argv[2])
print('#\t{} images to do photometry'.format(len(imlist)))
print('='*60)
# for i, img in enumerate(imlist): print('{}\t{}'.format(i, os.path.basename(img)))
for i, img in enumerate(imlist): print('{}\t{}'.format(i, img))
print('='*60)
#------------------------------------------------------------
# inim = 'Calib-LOAO-NGC2207-20201207-074859-R-120_gregister-com.fits'
tblist = []
trtblist = []
failist = []
# j = 0
# inim = imlist[j]

if __name__ == '__main__':
	pool = multiprocessing.Pool(processes=ncore)
	pool.map(phot_routine, imlist)
	pool.close()
	pool.join()

print('='*60)
print('5. RESULT TO FILES')
# comtbl = vstack(tblist)
# comtbl = vstack(sorted(glob.glob('{}/*.phot.summary.dat'.format(path_base))))
comtbl = vstack([ascii.read(tbl) for tbl in sorted(glob.glob('{}/*.phot.summary.dat'.format(path_base)))])
if 'com' in imkey:
	tbl_outname = 'phot.com'
elif '0.fits' in imkey:
	tbl_outname = 'phot.single'
else:
	tbl_outname = 'phot'

comtbl.write('{}/{}.summary.dat'.format(path_base, tbl_outname), format='ascii.tab', overwrite=True)
if transient == True:
	targetbl = vstack(trtblist)
	targetbl.write('{}/{}.dat'.format(path_base, tbl_outname), format='ascii.tab', overwrite=True)

print('{} failed images'.format(len(failist)))
f = open('{}/{}.fail.list'.format(path_base, tbl_outname), 'w')
f.write('image\n')
for failim in failist: f.write('{}\n'.format(failim))
f.close()
#------------------------------------------------------------
#	Time
#------------------------------------------------------------
delt = time.time() - starttime
dimen = 'seconds'
if delt > 60.:
	delt = delt/60.
	dimen = 'mins'
if delt > 60.:
	delt = delt/60.
	dimen = 'hours'
print('PHOTOMETRY IS DONE.\t({} {})'.format(round(delt, 3), dimen))
