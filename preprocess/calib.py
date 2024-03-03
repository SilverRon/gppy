#============================================================
#   CALIBTRATION ROUTINE FUNCTION FOR IMSNG TELESCOPES
#	LOAO, DOAO, SOAO, ...
#   18.12.04	UPDATED BY Gregory S.H. Paek
#	19.08.08	CREATED BY Gregory S.H. Paek	(FUNCTIONALIZATION)
#============================================================
#	FUNCTION
#============================================================
import os, glob
import numpy as np
from astropy import units as u
from astropy.io import fits
import astropy.io.ascii as ascii
import sys
from astropy.nddata import CCDData
import ccdproc
import matplotlib.pyplot as plt
sys.path.append('..')
from util import tool

def bn_median(masked_array, axis=None):
    """
    Perform fast median on masked array

    Parameters

    masked_array : `numpy.ma.masked_array`
        Array of which to find the median.

    axis : int, optional
        Axis along which to perform the median. Default is to find the median of
        the flattened array.
    """
    import numpy as np
    import bottleneck as bn
    data = masked_array.filled(fill_value=np.NaN)
    med = bn.nanmedian(data, axis=axis)
    # construct a masked array result, setting the mask from any NaN entries
    return np.ma.array(med, mask=np.isnan(med))

def folder_check(datalist, path_raw, path_factory):
	'''
	CHECK NEW OBSERVATION DATA FOLDER
	'''
	newlist = []
	#	IS THERE NEW FILE?
	for data in glob.glob(path_raw+'/*'):
		if data not in datalist:
			print('NEW DATA\t: {}'.format(data))
			newlist.append(data)
	#	MOVE NEW FILE TO FACTORY FOLDER
	for data in newlist:
		#	PREVENT OVERWRITE
		if path_factory+'/'+os.path.basename(data) in glob.glob(path_factory+'/*'):
			rmcom = 'rm -rf '+path_factory+'/'+os.path.basename(data)
			print(rmcom) ; os.system(rmcom)
		cpcom = 'cp -r '+data+' '+path_factory
		print(cpcom) ; os.system(cpcom)
	#	CHANGE TO WORKING PATH
	for i in range(len(newlist)):
		newlist[i] = path_factory+'/'+os.path.basename(newlist[i])
	return newlist
#------------------------------------------------------------
def changehdr(inim, where, what):
	from astropy.io import fits
	data, hdr = fits.getdata(inim, header = True)
	hdr[where] = what
	fits.writeto(inim, data, hdr, clobber=True)
#------------------------------------------------------------
def correcthdr_routine(path_data, hdrtbl, obs):
	'''
	1.	NGC337		-> NGC0337
		ngc0337		-> NGC0337
	'''
	from astropy.io import fits
	from astropy.time import Time
	comment = '-'*60+'\n' \
			+ 'CHANGE TO UNIVERSAL HEADER SERIES ...\n' \
			+ '-'*60+'\n'
	print(comment)
	objfilterlist = []
	objexptimelist = []
	flatfilterlist = []
	darkexptimelist = []
	imlist = []
	for ims in ('{}/*.fits'.format(path_data), '{}/*.fit'.format(path_data), '{}/*.fts'.format(path_data)):
		imlist.extend(sorted(glob.glob(ims)))
	# for inim in glob.glob(path_data+'/*.fit*'):
	for inim in imlist:
		print(inim)
		# tool_tbd.puthdr(inim, hdrkey='OBS', hdrval=obs, hdrcomment='observation location')
		fits.setval(inim, 'OBS', value=obs, comment='observation location', ignore_missing_end=True)
		data, hdr   = fits.getdata(inim, header=True)
		#	CHECK IMAGETYP HDR
		if hdr['IMAGETYP'] == 'Light': hdr['IMAGETYP'] = 'object'
		if hdr['IMAGETYP'] == 'zero' : hdr['IMAGETYP'] = 'Bias'
		if hdr['IMAGETYP'] == 'Dark' : hdr['IMAGETYP'] = 'dark'
		if hdr['IMAGETYP'] == 'Flat Field' : hdr['IMAGETYP'] = 'FLAT'
		#	CHECK FILTER HDR
		if hdr['IMAGETYP'] in ['object', 'Bias', 'dark',]:
			if 'FILTER' in hdr.keys():
				if ((hdr['FILTER'] == 'U101') | (hdr['FILTER'] == 1)): hdr['FILTER'] = 'U'
				if ((hdr['FILTER'] == 'B102') | (hdr['FILTER'] == 2)): hdr['FILTER'] = 'B'
				if ((hdr['FILTER'] == 'V103') | (hdr['FILTER'] == 3)): hdr['FILTER'] = 'V'
				if ((hdr['FILTER'] == 'R104') | (hdr['FILTER'] == 4)): hdr['FILTER'] = 'R'
				if ((hdr['FILTER'] == 'I105') | (hdr['FILTER'] == 5)): hdr['FILTER'] = 'I'
				if (hdr['FILTER'] == '0') | (hdr['FILTER'] == 0):
					hdr['FILTER'] = 'R'
					print('PLEASE CHECK SOAO FILTER HDR (FILTER:0?)')
		#	CHECK DATE-OBS
		if 'T' not in hdr['DATE-OBS']: hdr['DATE-OBS'] = hdr['DATE-OBS']+'T'+hdr['TIME-OBS']
		t = Time(hdr['DATE-OBS'], format='isot')
		hdr['JD'] = t.jd
		hdr['MJD'] = t.mjd
		#	CHECK OBJECT HDR
		if 'OBJECT' in hdr.keys():
			if 'ngc' in hdr['OBJECT']:
				while len(hdr['OBJECT'])<7:
					head = hdr['OBJECT'][0:3]
					tail = hdr['OBJECT'][3: ]
					tail = '0'+tail
					hdr['OBJECT'] = head+tail
			_obj = fits.getheader(inim)['object']
			#	GW event
			if (_obj[:2] == 'S2') | (_obj[:3] == 'MS2'):
				obj = _obj
				pass
			#	Other target
			else:
				obj = _obj.upper()
			hdr['OBJECT'] = obj
		else:
			# No OBJECT Header
			pass
		#	CHECK USER SETTING HDR
		for i in range(len(hdrtbl)):
			key = hdrtbl['key'][i]
			val = hdrtbl['val'][i]
			newval = hdrtbl['newval'][i]
			if key in hdr.keys():
				if hdr[key] == val:
					# print(hdr[key], key, newval)
					hdr[key] = newval
		fits.writeto(inim, data, hdr, overwrite=True)
		if hdr['IMAGETYP'] == 'object':
			objfilterlist.append(hdr['FILTER'])
			objexptimelist.append(hdr['EXPTIME'])
		if hdr['IMAGETYP'].upper() == 'FLAT':
			flatfilterlist.append(hdr['FILTER'])
		if hdr['IMAGETYP'].upper() == 'DARK':
			darkexptimelist.append(hdr['EXPTIME'])
	#	OBJECT FILTER, FLAT FILTER
	objfilterlist = list(set(objfilterlist))#;		objfilterlist.sort()
	objexptimelist = list(set(objexptimelist))#;		objexptimelist.sort()
	flatfilterlist = list(set(flatfilterlist))#;		flatfilterlist.sort()
	darkexptimelist = list(set(darkexptimelist))#;	darkexptimelist.sort()
	return objfilterlist, objexptimelist, flatfilterlist, darkexptimelist, t
#------------------------------------------------------------
def isot_to_mjd(time):      # 20181026 to 2018-10-26T00:00:00:000 to MJD form
	from astropy.time import Time
	yr  = time[0:4]         # year
	mo  = time[4:6]         # month
	da  = time[6:8]         # day
	isot = yr+'-'+mo+'-'+da+'T00:00:00.000'         #	ignore hour:min:sec
	t = Time(isot, format='isot', scale='utc')		#	transform to MJD
	return t.mjd
#------------------------------------------------------------
def call_images(path_data):
	import ccdproc
	from astropy.io import ascii
	images = ccdproc.ImageFileCollection(path_data, keywords='*')
	#	DATA SUMMARY
	ascii.write(images.summary, path_data+'/summary.txt', format='fixed_width_two_line', overwrite=True)
	return images
#------------------------------------------------------------
def getobsinfo(obs, obstbl):
	import numpy as np
	from astropy import units as u
	indx_obs = np.where(obs == obstbl['obs'])
	#	CCD INFORMATION
	gain = obstbl['gain'][indx_obs].item() * u.electron / u.adu
	rdnoise = obstbl['RDnoise'][indx_obs].item() * u.electron
	pixscale = obstbl['pixelscale'][indx_obs].item()
	fov = obstbl['fov'][indx_obs].item()
	# dark        = float(obstbl['dark'][indx_obs][0])
	obsinfo = dict(	
					gain=gain,
					rdnoise=rdnoise,
					pixscale=pixscale,
					fov=fov,
					)
	return obsinfo
#------------------------------------------------------------
#	CCD PROCESSING
#------------------------------------------------------------
def master_zero(images, fig=False):
	import ccdproc
	import os
	from astropy.nddata import CCDData
	import matplotlib.pyplot as plt
	comment     = '-'*60+'\n' \
				+ 'MAKING MASTER ZERO\n' \
				+ '-'*60+'\n'
	print(comment)
	path_data = images.location
	#	HEADER FOR MASTER FRAME
	zerolist   = []
	for hdu, fname in images.hdus(imagetyp='Bias', return_fname=True):
		meta = hdu.header
		zerolist.append(ccdproc.CCDData(data=hdu.data, meta=meta, unit="adu"))
	#	HEADER FOR MASTER FRAME
	n = 0
	for hdu, fname in images.hdus(imagetyp='Bias', return_fname=True):	
		n += 1
		newmeta = meta
		newmeta['FILENAME'] = 'zero.fits'
		newmeta['COMB{}'.format(n)] = fname
	print('{} ZERO IMAGES WILL BE COMBINED.'.format(len(zerolist)))
	zeros = ccdproc.Combiner(zerolist)
	# mzero = zeros.median_combine()
	mzero = zeros.median_combine()
	mzero.header  = newmeta
	#	SAVE MASTER FRAME
	if '{}/zero.fits'.format(path_data) in glob.glob(path_data+'/*'):
		os.system('rm {}/zero.fits'.format(path_data))
	mzero.write('{}/zero.fits'.format(path_data))
	if fig != False:
		imstats = lambda dat: (dat.min(), dat.max(), dat.mean(), dat.std())
		zero_min, zero_max, zero_mean, zero_std = imstats(np.asarray(mzero))
		plt.figure(figsize=(15, 15))
		plt.imshow(mzero, vmax=zero_mean + 4*zero_std, vmin=zero_mean - 4*zero_std)
		plt.savefig(path_data+'/zero.png')
	return mzero
#------------------------------------------------------------
def master_zero_bad(images, fig=False):
	comment     = '-'*60+'\n' \
				+ 'MAKING MASTER ZERO\n' \
				+ '-'*60+'\n'
	print(comment)
	path_data = images.location
	#	HEADER FOR MASTER FRAME
	# fname_sample = images.filter(imagetyp='Bias Frame',).files[0]
	# with fits.open(fname_sample) as hdul:
	# 	meta = hdul[0]
	# hdulist = []
	# for hdu, fname in images.hdus(imagetyp='Bias Frame', return_fname=True):
	for hdu, fname in images.hdus(imagetyp='bias', return_fname=True):
		meta = hdu.header
		# hdulist.append(meta)
		break
	#	New Meta data based on the first image
	newmeta = meta.copy()
	imlist = images.filter(imagetyp='Bias Frame',).files
	for n, inim in enumerate(imlist):
		n += 1
		newmeta[f'COMB{n}'] = os.path.basename(inim)
	newmeta['FILENAME'] = 'zero.fits'
	print(f'{len(imlist)} ZERO IMAGES WILL BE COMBINED.')

	combined_data = bn_median(np.ma.array([fits.getdata(inim) for inim in imlist]), axis=0)
	mzero = CCDData(data=combined_data.data, meta=newmeta, unit='adu')

	#	SAVE MASTER FRAME
	if os.path.exists(f"{path_data}/zero.fits"):
		os.system(f'rm {path_data}/zero.fits')
	mzero.write(f'{path_data}/zero.fits')
	# if fig != False:
	# 	imstats = lambda dat: (dat.min(), dat.max(), dat.mean(), dat.std())
	# 	zero_min, zero_max, zero_mean, zero_std = imstats(np.asarray(mzero))
	# 	plt.figure(figsize=(15, 15))
	# 	plt.imshow(mzero, vmax=zero_mean + 4*zero_std, vmin=zero_mean - 4*zero_std)
	# 	plt.savefig(path_data+'/zero.png')
	return mzero
#------------------------------------------------------------
def master_dark(images, mzero, exptime, fig=False):
	# import ccdproc
	# import os
	# from astropy.nddata import CCDData
	# import matplotlib.pyplot as plt
	comment     = '-'*60+'\n' \
				+ 'MAKING MASTER DARK\n' \
				+ '-'*60+'\n'
	print(comment)
	path_data = images.location
	mdark_name = 'dark-{}.fits'.format(int(exptime))
	#	HEADER FOR MASTER FRAME
	zdarklist   = []
	for hdu, fname in images.hdus(imagetyp='dark', exptime=exptime, return_fname=True):
		meta = hdu.header
		dark = ccdproc.CCDData(data=hdu.data, meta=meta, unit="adu")
		zdark = ccdproc.subtract_bias(dark, mzero)
		zdark.meta['SUBBIAS'] = mzero.meta['FILENAME']
		zdarklist.append(ccdproc.CCDData(data=zdark.data, meta=meta, unit="adu"))
	#	HEADER FOR MASTER FRAME
	n = 0
	for hdu, fname in images.hdus(imagetyp='dark', exptime=exptime, return_fname=True):	
		n += 1
		newmeta = meta
		newmeta['FILENAME'] = mdark_name
		newmeta['COMB{}'.format(n)] = fname
	# print('{} DARK({} sec) IMAGES WILL BE COMBINED.'.format(len(zdarklist)), exptime)
	darks = ccdproc.Combiner(zdarklist)
	mdark = darks.median_combine()
	mdark.header  = newmeta
	#	SAVE MASTER FRAME
	if '{}/{}'.format(path_data, mdark_name) in glob.glob(path_data+'/*'):
		os.system('rm {}/{}'.format(path_data, mdark_name))
	mdark.write('{}/{}'.format(path_data, mdark_name))
	if fig != False:
		imstats = lambda dat: (dat.min(), dat.max(), dat.mean(), dat.std())
		zero_min, zero_max, zero_mean, zero_std = imstats(np.asarray(mzero))
		plt.figure(figsize=(15, 15))
		plt.imshow(mzero, vmax=zero_mean + 4*zero_std, vmin=zero_mean - 4*zero_std)
		plt.savefig(path_data+'/zero.png')
	return mdark
#------------------------------------------------------------
def master_flat(images, mzero, filte, mdark=False, fig=False):
	import os
	import ccdproc
	from astropy.nddata import CCDData
	import matplotlib.pyplot as plt
	comment     = '-'*60+'\n' \
				+ 'MAKING MASTER FLAT\n' \
				+ '-'*60+'\n'
	print(comment)
	path_data = images.location
	if mdark == False:
		#	FLAT - MZERO = zFLAT
		print('SUBTRACT MASTER ZERO FROM {} FLAT'.format(filte))
		outname = 'n'+filte+'.fits'
		zflatlist = []
		for hdu, fname in images.hdus(imagetyp='flat', filter=filte, return_fname=True):
			meta = hdu.header
			flat = ccdproc.CCDData(data=hdu.data, meta=meta, unit="adu")
			zflat = ccdproc.subtract_bias(flat, mzero)
			zflat.meta['SUBBIAS'] = mzero.meta['FILENAME']
			zflatlist.append(ccdproc.CCDData(data=zflat.data, meta=meta, unit="adu"))
	else:
		#	FLAT - MZERO = zFLAT
		print('SUBTRACT MASTER ZERO & DARK FROM {} FLAT'.format(filte))
		outname = 'n'+filte+'.fits'
		zflatlist = []
		for hdu, fname in images.hdus(imagetyp='flat', filter=filte, return_fname=True):
			meta = hdu.header
			flat = ccdproc.CCDData(data=hdu.data, meta=meta, unit="adu")
			zflat = ccdproc.subtract_bias(flat, mzero)
			dzflat = ccdproc.subtract_dark(	ccd=zflat, master=mdark,
											# data_exposure=zflat.meta['EXPTIME'],
											exposure_time='EXPTIME',
											exposure_unit=u.second,
											# dark_exposure=mdark.meta['EXPTIME'],
											scale=True)
			meta['SUBBIAS'] = mzero.meta['FILENAME']
			meta['SUBDARK'] = mdark.meta['FILENAME']
			zflatlist.append(ccdproc.CCDData(data=dzflat.data, meta=meta, unit="adu"))
	#	HEADER FOR MASTER FRAME
	newmeta = meta
	newmeta['FILENAME'] = outname
	n = 0
	for hdu, fname in images.hdus(imagetyp='flat', filter=filte, return_fname=True):
		n += 1
		newmeta['COMB{}'.format(n)] = fname
	#	zFLATs -> MASTER FRAME
	print('{} {} FLAT IMAGES WILL BE COMBINED.'.format(len(zflatlist), filte))
	flat_combiner = ccdproc.Combiner(zflatlist)
	flat_combiner.minmax_clipping()
	def scaling_func(arr): return 1/np.ma.median(arr)
	flat_combiner.scaling = scaling_func
	mflat = flat_combiner.median_combine()
	#	SAVE MASTER FRAME
	mflat.header = newmeta
	if '{}/{}'.format(path_data, outname) in glob.glob(path_data+'/*'):
		os.system('rm {}/{}'.format(path_data, outname))
	mflat.write(path_data+'/'+outname)
	# mflat_electron = ccdproc.gain_correct(mflat, gain=gain)
	if fig != False:
		imstats = lambda dat: (dat.min(), dat.max(), dat.mean(), dat.std())
		f_min, f_max, f_mean, f_std = imstats(np.asarray(mflat))
		plt.figure(figsize=(15, 15))
		plt.imshow(mflat, vmin=f_mean-5*f_std, vmax=f_mean+5*f_std)
		plt.savefig(path_data+'/'+outname[:-5]+'.png')

	return mflat
#------------------------------------------------------------
def calibration(images, mzero, mflat, filte, mdark=False):
	import ccdproc
	comment     = '-'*60+'\n' \
				+ 'OBJECT CORRECTION\n' \
				+ '-'*60+'\n'
	print(comment)
	path_data = images.location
	objlist = []
	for hdu, fname in images.hdus(imagetyp='object', filter=filte, return_fname=True):
		meta = hdu.header
		objlist.append(meta['object'])
		obj = ccdproc.CCDData(data=hdu.data, meta=meta, unit="adu")
		if mdark == False:
			#	ZERO CORRECTION
			zobj = ccdproc.subtract_bias(obj, mzero)
			meta['SUBBIAS'] = mzero.meta['FILENAME']
		else:
			#	ZERO & DARK CORRECTION
			zobj = ccdproc.subtract_bias(obj, mzero)
			meta['SUBBIAS'] = mzero.meta['FILENAME']
			zobj = ccdproc.subtract_dark(	ccd=zobj, master=mdark,
											# dark_exposure=mdark.meta['EXPTIME'],
											# data_exposure=dzobj.meta['EXPTIME'],
											exposure_time='EXPTIME',
											exposure_unit=u.second,
											scale=True)
			meta['SUBDARK'] = mdark.meta['FILENAME']
		#	FLAT CORRECTION
		fzobj = ccdproc.flat_correct(zobj, mflat)
		meta['DIVFLAT'] = mflat.meta['FILENAME']
		fzobj.header = meta
		# fzobj.write(path_data+'/fz'+fname)
		fzobj.write('{}/fz{}'.format(path_data, fname), overwrite=True)
#------------------------------------------------------------
def astrometry(inim, pixscale, ra=None, dec=None, fov=1, cpulimit=60):
	'''
	ra : hh:mm:ss
	dec : dd:mm:ss
	fov [deg]
	'''
	import os
	upscl = str(pixscale + pixscale*0.10)
	loscl = str(pixscale - pixscale*0.10)
	outname = os.path.dirname(inim)+'/a'+os.path.basename(inim)
	# com     ='/usr/local/astrometry/bin/solve-field/solve-field '+inim \

	# com     ='solve-field '+inim \
	com     ='/usr/local/astrometry/bin/solve-field '+inim \
			+' --scale-unit arcsecperpix --scale-low '+loscl+' --scale-high '+upscl \
			+' --no-plots --new-fits '+outname+' --overwrite --use-sextractor --cpulimit {}'.format(cpulimit)
	if ((ra != None) & (dec != None) & (fov != None)):
		com	= com + ' --ra {} --dec {} --radius {}'.format(ra, dec, fov)
	else:
		pass
	print(com); os.system(com)
#------------------------------------------------------------
def fnamechange(inim, obs, com='mv'):
	import os
	from astropy.io import fits
	hdr = fits.getheader(inim)
	dateobs = hdr['DATE-OBS']
	datestr = dateobs[0:4]+dateobs[5:7]+dateobs[8:10]
	timestr = dateobs[11:13]+dateobs[14:16]+dateobs[17:19]
	_objname = hdr['OBJECT']
	#	GW event
	if (_objname[:2] == 'S2') | (_objname[:3] == 'MS2'):
		objname = _objname
		pass
	else:
		objname = _objname.upper()
	filname = hdr['FILTER']  # str(hdr['FILTER'])
	exptime = str(int(hdr['EXPTIME']))
	newname = 'Calib-{}-{}-{}-{}-{}-{}.fits'.format(obs, objname, datestr, timestr, filname, exptime)
	if com == 'cp':
		cpcom = 'cp {} {}/{}'.format(inim, os.path.dirname(inim), newname)
		print(cpcom)
		os.system(cpcom)
	if com == 'mv':
		mvcom = 'mv {} {}/{}'.format(inim, os.path.dirname(inim), newname)
		print(mvcom)
		os.system(mvcom)
#------------------------------------------------------------
# def movecalib(inim, path_gal):
# 	import os, glob
# 	img = os.path.basename(inim)
# 	part = img.split('-')
# 	obs = part[1]
# 	obj = part[2]
# 	filte = part[5]
# 	#	MAKE SAVE PATH
# 	path_goal = path_gal
# 	for folder in [obj, obs, filte]:
# 		path_goal = path_goal+'/'+folder
# 		mkcom = 'mkdir {}'.format(path_goal)
# 		print(mkcom); os.system(mkcom)
# 	mvcom = 'mv {} {}'.format(inim, path_goal)
# 	print(mvcom); os.system(mvcom)
#------------------------------------------------------------
def movecalib(inim, path_gal='/data3/IMSNG/IMSNGgalaxies', additional=False):
	import os
	'''
	path_gal = '/data3/IMSNG/IMSNGgalaxies'
	inim = '/data3/paek/factory/doao/20201209-1m-IMSNG/Calib-DOAO-NGC0000-20201209-113240-B-60.fits'
	inim = '/data3/paek/factory/doao/20201209-1m-IMSNG/Calib-DOAO-ESO000-123-20201209-113240-B-60.fits'
	'''
	part = os.path.basename(inim).split('-')
	#	(obj name) = 'NGC0000'
	if (len(part)==7) & (os.path.basename(inim).count('-')==6):
		obs = part[1]
		obj = part[2]
		filte = part[5]
	#	(obj name) = 'ESO000-123' containing additional '-'
	elif (len(part)==8) & (os.path.basename(inim).count('-')==7):
		obs = part[1]
		obj = part[2]
		filte = obj = '-'.join([part[2], part[3]])
	else:
		pass

	#	Path
	path_obj = '{}/{}'.format(path_gal, obj)
	path_obs = '{}/{}/{}'.format(path_gal, obj, obs)
	path_filte = '{}/{}/{}/{}'.format(path_gal, obj, obs, filte)
	if additional != False:
		path_filte = '{}/{}'.format(path_filte, additional)

	#	Command
	mkcom0 = 'mkdir {}'.format(path_obj)
	mkcom1 = 'mkdir {}'.format(path_obs)
	mkcom2 = 'mkdir {}'.format(path_filte)
	mvcom = 'mv {} {}'.format(inim, path_filte)

	#	Create folder
	if os.path.exists(path_obj):
		if os.path.exists(path_obs):
			if os.path.exists(path_filte):
				pass
			else:
				print(mkcom2)
				os.system(mkcom2)
		else:
			print(mkcom1)
			os.system(mkcom1)
	else:
		print(mkcom0)
		os.system(mkcom0)
		print(mkcom1)
		os.system(mkcom1)
		print(mkcom2)
		os.system(mkcom2)
	#	Move image
	# print(mvcom)
	os.system(mvcom)
#------------------------------------------------------------
def defringe(inim, dfim, dfdat):
	'''
	inim : image to remove fringe
	dfim : master fringe image
	dfdat : master fringe data
	'''
	data, hdr = fits.getdata(inim, header=True)
	#	Master fringe
	dataf, hdrf	= fits.getdata(dfim, header=True)

	master_fri = fringe_cal(dfim, dfdat)
	infr_list = []
	fscale = np.median(fringe_cal(inim, dfdat)/master_fri)
	# print(fscale)

	fri_scaled = dataf*fscale
	dfdata = data-fri_scaled
	outim = '{}/df{}'.format(os.path.dirname(inim), os.path.basename(inim))
	fits.writeto(outim, dfdata, hdr, overwrite=True)
	# print(np.median(data), np.median(fri_scaled), np.median(dfdata))

	return outim
#------------------------------------------------------------
def fringe_cal(inim, dfdat):
	data, hdr	= fits.getdata(inim, header=True)
	#subframe	= data[512:1537, 512:1537]
	#skymean, skymed, skysig		= bkgest_mask(inim)
	#subdata		= data-skymed

	dfr_list	= []
	dftbl = ascii.read(dfdat)

	xb1, xb2 = dftbl['xb']-5, dftbl['xb']+5
	yb1, yb2 = dftbl['yb']-5, dftbl['yb']+5

	xf1, xf2 = dftbl['xf']-5, dftbl['xf']+5
	yf1, yf2 = dftbl['yf']-5, dftbl['yf']+5

	for n in range(len(dftbl)):
		'''
		fringe_b= np.median(subdata[yb1[n]:yb2[n], xb1[n]:xb2[n]])
		fringe_f= np.median(subdata[yf1[n]:yf2[n], xf1[n]:xf2[n]])
		'''
		fringe_b= np.median(data[yb1[n]:yb2[n], xb1[n]:xb2[n]])
		fringe_f= np.median(data[yf1[n]:yf2[n], xf1[n]:xf2[n]])		
		dfringe	= fringe_b-fringe_f
		dfr_list.append(dfringe)
	return np.array(dfr_list)
#------------------------------------------------------------
def cr_removal(inim, outim, gain, rdnoise):
	'''
	inim 
	obs = 'LOAO'
	gain = 2.68
	rdnoise = 4.84
	'''
	import os
	from astroscrappy import detect_cosmics
	from astropy.io import fits
	import time

	data, hdr = fits.getdata(inim, header=True)

	param_cr = dict(
					indat=data,
					sigclip=4.5,
					sigfrac=0.3,
					objlim=5.0,
					gain=gain, readnoise=rdnoise, 
					# pssl=0.0,
					niter=4,
					sepmed=True,
					cleantype='meanmask',
					fsmode='median',
					psfmodel='gauss',
					psffwhm=hdr['seeing'],
					#	Notouch
					inmask=None,
					satlevel=65536.0,
					psfsize=7, psfk=None, psfbeta=4.765,
					verbose=False
					)

	time_st = time.time()
	mcrdata, crdata = detect_cosmics(**param_cr)
	fits.writeto(outim, crdata, hdr, overwrite=True)
	# tool_tbd.puthdr(outim,  )
	time_delta = time.time() - time_st
	# fits.writeto('{}/mcr{}'.format(os.path.dirname(inim), os.path.basename(inim)), mcrdata, hdr, overwrite=True)
	#------------------------------------------------------------
	# outbl = Table()
	indx_cr = np.where(mcrdata == True)
	# outbl['number'] = np.arange(1, len(indx_cr[0])+1, 1)
	# outbl['x'] = indx_cr[1]+1
	# outbl['y'] = indx_cr[0]+1
	# outbl['val_ori'] = data[indx_cr]
	# outbl['val_cor'] = crdata[indx_cr]
	# print('Remove {} cosmic-ray for {} [{} sec]'.format(len(outbl), inim, round(time_delta, 3)))
	print('Remove {} cosmic-ray for {} [{} sec]'.format(len(indx_cr[0]), inim, round(time_delta, 3)))
	# tool_tbd.puthdr(outim, 'CR', len(indx_cr[0]), '# of removed cosmic ray by Astroscrappy')
	fits.setval(outim, 'CR', value=len(indx_cr[0]), comment='# of removed cosmic ray by Astroscrappy', ignore_missing_end=True)
	# return outbl

