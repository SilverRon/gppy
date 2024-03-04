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
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
sys.path.append('..')
sys.path.append('/home/gecko/gppy')
from phot import gpphot
from util import query
from util import tool
from phot import gcurve
from datetime import date
import time
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

#============================================================
#	USER SETTING
#============================================================
#	PATH
#------------------------------------------------------------
try:
	path_base = sys.argv[1]
except:
	path_base = '.'
# print(path_base)
path_refcat	= '/data3/paek/factory/refcat'
path_obs = '/home/paek/table/obs.dat'
path_config = '/home/paek/config'
# path_target = '/home/sonic/Research/gppy/table/transient.dat'
path_target = './transient.dat'
path_gphot = '{}/gphot.config'.format(path_base)
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
#------------------------------------------------------------
# inmagkey = 'MAG_APER'
# inmagerkey = 'MAGERR_APER'
# inmagkey = 'MAG_AUTO'
# inmagerkey = 'MAGERR_AUTO'
#------------------------------------------------------------
imlist = sorted(glob.glob(imkey))
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
for j, inim in enumerate(imlist):
	#------------------------------------------------------------
	#	INFO. from file name
	#------------------------------------------------------------
	part = os.path.splitext(os.path.basename(inim))[0].split('-')
	head = os.path.splitext(inim)[0]
	obs = part[1]
	obj = part[2]
	refmagkey = part[5]
	refmagerkey = refmagkey+'err'
	#------------------------------------------------------------
	obsdict = tool.getccdinfo(obs, path_obs)
	gain = obsdict['gain']
	pixscale = obsdict['pixelscale']
	if obs == 'SAO_C361K':
		pixscale = pixscale*fits.getheader(inim)['XBINNING']

	fov = obsdict['fov']
	# rdnoise = obsdict['readoutnoise']
	#------------------------------------------------------------
	#	OUTPUT NAMES
	cat = head+'.cat'
	cat_gc = head+'.gcurve.cat'
	seg = head+'.seg.fits'
	bkg = head+'.bkg.fits'
	sub = head+'.sub.fits'
	psf = head+'.psf'
	aper = head+'.aper.fits'
	#	PHOTOMETRY NAMES
	param = path_config+'/gregoryphot.param'
	conv = path_config+'/gregoryphot.conv'
	nnw = path_config+'/gregoryphot.nnw'
	conf = path_config+'/gregoryphot.sex'
	#------------------------------------------------------------
	print('-'*60)
	print('[{}/{}] {}'.format(j+1, len(imlist), inim))
	print('{}\t{} in {}-band'.format(obs, obj, refmagkey))
	print('-'*60)
	#------------------------------------------------------------
	#	INFO. from file header
	#------------------------------------------------------------
	hdr = fits.getheader(inim)
	#	RA, Dec center for reference catalog query
	xcent, ycent= hdr['NAXIS1']/2., hdr['NAXIS2']/2.
	w = WCS(inim)
	racent, decent = w.all_pix2world(xcent, ycent, 1)
	racent, decent = racent.item(), decent.item()
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
	seeing = hdr['seeing']*u.arcsecond
	peeing = hdr['peeing']*u.pix
	seeing_input = str(seeing.value)

	optaper = hdr['APERPIX']
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
						SEEING_FWHM = str(seeing.value),
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
	print(com)
	os.system(com)
	# sexout = subprocess.getoutput(com)
	# line = [s for s in sexout.split('\n') if 'RMS' in s]
	# skymed, skysig = float(line[0].split('Background:')[1].split('RMS:')[0]), float(line[0].split('RMS:')[1].split('/')[0])
	skymed = hdr['SKYVAL']
	skysig = hdr['SKYSIG']
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
	intbl = setbl
	intbl['FWHM_IMAGE'] = intbl['FWHM_IMAGE'] 
	intbl['FWHM_WORLD'] = intbl['FWHM_WORLD'].to(u.arcsecond)
	ellipticity = hdr['ELLIP']
	elongation = hdr['ELONG']
	#------------------------------------------------------------
	#	Header
	#------------------------------------------------------------
	tool.puthdr(inim, 'SUBTIME', date.today().isoformat(), hdrcomment='PHTOMETRY TIME [KR]')
	#------------------------------------------------------------
	#	ZEROPOINT CALCULATION
	#------------------------------------------------------------
	onetbl = Table()
	onetbl['obs'] = [obs]
	onetbl['field'] = hdr['OBJECT']
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
		print('ZERO POINT CALCULATION for {}'.format(inmagkey))
		try:
			if inmagkey == 'MAG_AUTO':
				pass
			else:
				onetbl[inmagkey.replace('MAG_', '').replace('APER', 'aperture')] = aper
			
			zp = hdr['ZP_{}'.format(k)]
			zper = hdr['ZPER_{}'.format(k)]
			onetbl[inmagkey.replace('MAG', 'STDNUMB').lower()] = 0
			onetbl[inmagkey.replace('MAG', 'ZP').lower()] = zp
			onetbl[inmagkey.replace('MAG', 'ZPER').lower()] = zper
			onetbl[inmagkey.replace('MAG', 'UL3SIG').lower()] = hdr['UL3_{}'.format(k)]
			onetbl[inmagkey.replace('MAG', 'UL5SIG').lower()] = hdr['UL5_{}'.format(k)]
			#------------------------------------------------------------
			#	NORMAL PHOTOMETRY
			#------------------------------------------------------------
			# print('4-3. PHOTOMETRY TABLE')
			intbl[inmagkey.lower()] = zp + intbl[inmagkey]
			intbl[inmagerkey.lower()] = tool.sqsum(zper, intbl[inmagerkey])
		except:
			print('Maybe ZP for {} is insufficient.'.format(inmagkey))
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
	intbl.write('{}.phot_sub.cat'.format(head), format='ascii.tab', overwrite=True)
	onetbl.write('{}.phot_sub.summary.dat'.format(head), format='ascii.tab', overwrite=True)
	tblist.append(onetbl)
	plt.close('all')

print('='*60)
print('5. RESULT TO FILES')
comtbl = vstack(tblist)
if 'com' in imkey:
	tbl_outname = 'phot_sub.com'
elif '0.fits' in imkey:
	tbl_outname = 'phot_sub.single'
else:
	tbl_outname = 'phot_sub'

comtbl.write('{}/{}.summary.dat'.format(path_base, tbl_outname), format='ascii.tab', overwrite=True)
if transient == True:
	targetbl = vstack(trtblist)
	targetbl.write('{}/{}.dat'.format(path_base, tbl_outname), format='ascii.tab', overwrite=True)

print('{} failed images'.format(len(failist)))
f = open('{}/{}.fail_sub.list'.format(path_base, tbl_outname), 'w')
f.write('image\n')
for failim in failist: f.write('{}\n'.format(failim))
f.close()

print('PHOTOMETRY IS DONE.\t({} seconds)'.format(round((time.time() - starttime), 3)))
