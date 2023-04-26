#   ZEROPOINT CALCULATION
#   19.03.09    MERRY CHRISTMAS!
#   MADE BY GREGORY S.H. PAEK
#	UPDATE : 20.01.16
#=========================================================================#
# import os, glob, subprocess
import os
# print(os.getcwd())
import sys
sys.path.append('..')
from util import tool
import numpy as np
# from astropy.table import Table
from astropy.table import vstack
import matplotlib.pyplot as plt
# from numba import jit
#-------------------------------------------------------------------------#
#	Bottleneck for faster processing
#-------------------------------------------------------------------------#
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
    import bottleneck as bn
    data = masked_array.filled(fill_value=np.NaN)
    med = bn.nanmedian(data, axis=axis)
    # construct a masked array result, setting the mask from any NaN entries
    return np.ma.array(med, mask=np.isnan(med))
#-------------------------------------------------------------------------#
def which_obs(obs, path_and_file):
	"""
	=====================================================================
	GIVE CCD INFORMATION
	---------------------------------------------------------------------
	INPUT   :   observatory name, and path+table_name(ascii format)
				i.e.) path_and_file='/home/sonic/Research/table/obs.txt'
	OUTPUT  :   GAIN        []
				Pixel Scale ["/pix]
	---------------------------------------------------------------------
	obs         ccd         gain    RDnoise     dark    pixelscale
	SOAO        FLI         1.43    13.60       0.0     0.44540
	LOAO        E2V         2.68    4.84        0.0     0.794
	LOAO_FLI    FLI         3       15.0        0.0     1.28
	DOAO_sophia PIsophia    0.928   7.63        0.0     0.3855
	DOAO        FLI         1.27    14.0        0.005   0.464
	oldBOAO     UNKNOWN     8.9     20.0        0.0     1.7
	BOAO        UNKNOWN     5.0     2.0         0.0     0.2145
	BOAO        UNKNOWN     5.0     5.0         0.0     0.2145    
	SAO         SBIG        1.35    9.0         0.0     0.31
	30inch      UNKNOWN     1.60    5.87        0.0     1.3553
	CCA250      MLI16803    1.46    10.35       0.0     2.06
	MAIDANAK    SNU4kCAM    1.45    4.7         0.0     0.266
	MAIDANAK    UNKNOWN     5.40    9.0         0.0     0.266
	LSGT        SNUCAMII    1.15    6.0         0.2     0.92
	UKIRT       WFCAM       99.0    0.0         0.0     0.4
	=====================================================================
	"""
	import numpy as np
	from astropy.io import ascii
	obsinfo     = ascii.read(path_and_file)
	indx_obs    = np.where( obs == obsinfo['obs'] )
	gain        = obsinfo['gain'][indx_obs]
	pixscale    = obsinfo['pixelscale'][indx_obs]
	return gain, pixscale

#-------------------------------------------------------------------------#
def image_list(imlist):
	"""
	INPUT   :   imlist = glob.glob('Calib-*.fits')
	OUTPUT  :   observatory list
				object list
				fillter list
	"""
	obslist = []
	objlist = []
	fillist = []
	for img in imlist:
		sp  = img.split('-')
		obslist.append(sp[1])
		objlist.append(sp[2])
		fillist.append(sp[5])
	obslist = list(set(obslist))
	objlist = list(set(objlist))
	fillist = list(set(fillist))
	return obslist, objlist, fillist


#-------------------------------------------------------------------------#
def secom(inim, gain, pixscale, zp=0, seeing=3, det_sigma=3, backsize=str(64), backfiltersize=str(3), psf=False, dual=False, detect='detection.fits', check=False):
	"""
	SourceEXtractor
	APERTURE    3", 5", 7",
				1.0seeing, 1.2seeing ,1.5seeing ,1.7seeing ,2.0seeeing
	INPUT   :   (image).fits
				aperture    []
				seeing_fwhm [pixel]
	OUTPUT  :   no return
				.cat
	"""
	import numpy as np
	import os
	from astropy.io import ascii
	#   FILE CHECK
	#	CONFIG FILES (USER BASE PATH)
	configfile      = '/home/sonic/Research/yourpy/config/targetphot.sex'
	paramfile       = '/home/sonic/Research/yourpy/config/targetphot.param'
	nnwfile		    = '/home/sonic/Research/yourpy/config/targetphot.nnw'
	convfile	    = '/home/sonic/Research/yourpy/config/targetphot.conv'
	try:
		comment = 'SourceEXtractor START\n' \
				+ 'IMAGE\t\t: '+inim+'\n' \
				+ 'GAIN\t\t: '+str(gain)+'\n' \
				+ 'PIXSCALE\t: '+str(pixscale)+'\n' \
				+ 'DETECTION SIGMA\t: '+str(det_sigma)+'\n' \
				+ 'PARAM\t\t: '+paramfile+'\n' \
				+ 'BACKSIZE\t: '+backsize+'\n' \
				+ 'BACKFILTER\t: '+backfiltersize+'\n' \
				+ 'CONFIG\t\t: '+configfile+'\n' \
				+ 'NNW\t\t: '+nnwfile+'\n' \
				+ 'CONVOLVE\t: '+convfile
		print(comment)
	except:
		comment = 'CHECK configfile/paramfile/nnewfile/convfile or others.'
		print(comment)
	#   FILE NAME
	cat     = inim[:-5]+'.cat'
	seg     = inim[:-5]+'.seg.fits'
	bkg     = inim[:-5]+'.bkg.fits'
	sub     = inim[:-5]+'.sub.fits'
	impsf   = inim[:-5]+'.psf'
	aper    = inim[:-5]+'.aper.fits'

	#   BASIC INFO.
	det_area        = 5
	det_thresh      = det_sigma/np.sqrt(det_area)
	detecminarea    = str(det_area)
	detectthresh    = str(det_thresh)
	#   OPTION
	aperture        = '%.2f'%(3./pixscale)+','+'%.2f'%(5./pixscale)+','+'%.2f'%(7./pixscale)+','+'%.2f'%(seeing)+','+'%.2f'%(1.2*seeing)+','+'%.2f'%(1.5*seeing)+','+'%.2f'%(1.7*seeing)+','+'%.2f'%(2.0*seeing)
	option0	= ' -MAG_ZEROPOINT {} '.format(zp)
	option1 = ' -CATALOG_NAME '+cat+' -PARAMETERS_NAME '+paramfile
	option2 = ' -DETECT_MINAREA '+detecminarea+' -DETECT_THRESH '+detectthresh \
			+' -FILTER Y '+'-FILTER_NAME '+convfile
	option3 = ' -PHOT_APERTURES '+aperture \
			+' -GAIN '+'%.1f'%(gain)+' -PIXEL_SCALE '+'%.1f'%(pixscale)
	#option4 = ' -SEEING_FWHM '+'%.3f'%(seeing)+' -STARNNW_NAME '+nnwfile
	option4 = ' -SEEING_FWHM '+'%.3f'%(seeing)+' -STARNNW_NAME '+nnwfile
	option5 = ' -BACK_SIZE '+ backsize \
			+ ' -BACK_FILTERSIZE '+ backfiltersize+' -BACKPHOTO_TYPE LOCAL'
	option6 = ' -CHECKIMAGE_TYPE SEGMENTATION,APERTURES,BACKGROUND,-BACKGROUND'
	option7 = ' -CHECKIMAGE_NAME '+seg+','+aper+','+','+bkg+','+sub
	if psf	!= False:
		option8 = ' -PSF_NAME '+impsf
	else:
		option8 = ''
	#   COMMAND
	#   detect = detection.fits is fine image show target significantly and have good seeing and many stars i
	dualcom ='sex -c '+configfile+' '+detect+' , '+inim+' -CATALOG_NAME dual'+cat+' -PARAMETERS_NAME '+paramfile+ ' '+option0+' '+option2+' '+option3+' '+option4+' '+option5+' '+option6+' '+option7+' '+option8
	sglcom  ='sex -c '+configfile+' '+inim+' '+option0+' '+option1+' '+option2+' '+option3+' '+option4+' '+option5+' '+option6+' '+option7+' '+option8
	clearcom='sex -c '+configfile+' '+inim+' '+option0+' '+option1+' '+option2+' '+option3+' '+option4+' '+option5+' '+option8
	if dual == False    :
		if check == False:
			os.system(clearcom)
		else:
			os.system(sglcom)
	else                : os.system(dualcom)
		
	secat   = ascii.read(cat)
	return secat, cat
#-------------------------------------------------------------------------#
def psfex(inim, pixscale):
	"""
	PSfextractor
	INPUT   :   (image).fits
	OUTPUT  :   FWHM    [pixel]
				FWHM    [arcsec]
	"""
	import os
   
	#   FILE CHECK
	#	CONFIG FILES (USER BASE PATH)
	psfexconf_prese_conf    = '/home/sonic/Research/yourpy/config/prepsfex.sex'
	psfexconf_prese_param   = '/home/sonic/Research/yourpy/config/prepsfex.param'
	psfexconf_psfex_conf    = '/home/sonic/Research/yourpy/config/default.psfex'
	psfexconf_psfex_conv    = '/home/sonic/Research/yourpy/config/default.conv'
	try:
		comment = '\nPSFex START\n' \
				+ 'IMAGE\t\t: '+inim+'\n' \
				+ 'PRE_CONFIG\t: '+psfexconf_prese_conf+'\n' \
				+ 'PRE_PARAM\t: '+psfexconf_prese_param+'\n' \
				+ 'CONFIG\t\t: '+psfexconf_psfex_conf+'\n' \
				+ 'CONV\t\t: '+psfexconf_psfex_conv
		print(comment)
	except:
		comment = 'CHECK psfexconf_prese/psfexconf_prese_param/psfexconf_psfex_conf/psfexconf_psfex_conv OR OTHERS.'
		print(comment)

	#   FILE NAME
	cat     = inim[:-5]+'.cat'
	xml     = inim[:-5]+'.xml'
	snap    = 'snap_'+inim+'[100:125,100:125]'
	psf     = 'psf-'+inim
	#   OPTION
	presecom1   = psfexconf_prese_conf+" "+inim
	presecom2   = " -CATALOG_NAME "+cat
	presecom3   = " -FILTER_NAME " + psfexconf_psfex_conv + " -PARAMETERS_NAME " + psfexconf_prese_param
	#   COMMAND
	presecom    = "sex -c "+presecom1+presecom2+presecom3
	psfexcom    = "psfex -c "+psfexconf_psfex_conf+" "+cat
	os.system(presecom)
	os.system(psfexcom) 
	os.system('cp psfex.xml '+xml)
	#   SNAP IMAGE
	imcopycom   = 'imcopy '+snap+' '+psf
	print(imcopycom);   os.system(imcopycom)
	#   FWHM [pixel], FWHM [arcsec]
	fwhm_pix    = psfexxml(xml)
	fwhm_arcsec = round(fwhm_pix*pixscale, 3)
	comment     = '\n' \
				+ 'FILE NAME'+'\t'+': '+inim+'\n' \
				+ 'FWHM value'+'\t'+': '+str(fwhm_pix)+'\t'+'[pixel]'+'\n' \
				+ '\t'+'\t'+': '+str(fwhm_arcsec)+'\t'+'[arcsec]'+'\n'
	print(comment)
	return fwhm_pix, fwhm_arcsec
#------------------------------------------------------------
def sexcom(inim, param_insex, dualmode=False):
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

	
	sexcom_normal = 'sex -c {} {} '.format(param_sex['CONF_NAME'], inim)
	sexcom_dual = 'sex -c {} {} '.format(param_sex['CONF_NAME'], inim)
	for key in param_sex.keys():
		if key != 'CONF_NAME':
			sexcom_normal += '-{} {} '.format(key, param_sex[key])

	# print(sexcom_normal)
	# os.system(sexcom_normal)
	return sexcom_normal

def psfexxml(xmlfile):
	"""
	INPUT   :   .xml
	OUTPUT  :   FWHM    [pixel]
	"""
	from astropy.io.votable import parse
	from astropy.table import Table, Column, MaskedColumn
	votable     = parse(xmlfile)
	table       = votable.get_first_table()
	data        = table.array
	#   EXTRACT FWHM [pixel]
	fwhm        = data['FWHM_Mean'][0]
	fwhm        = round(fwhm, 3)
	return fwhm
#-------------------------------------------------------------------------#
#def matching(incat, refcat, sep=2.0):
def matching(intbl, reftbl, inra, indec, refra, refdec, sep=2.0):
	"""
	MATCHING TWO CATALOG WITH RA, Dec COORD. WITH python
	INPUT   :   SE catalog, SDSS catalog file name, sepertation [arcsec]
	OUTPUT  :   MATCED CATALOG FILE & TABLE
	"""
	import numpy as np
	import astropy.units as u
	# from astropy.table import Table, Column
	from astropy.coordinates import SkyCoord
	from astropy.io import ascii

	incoord		= SkyCoord(inra, indec, unit=(u.deg, u.deg))
	refcoord	= SkyCoord(refra, refdec, unit=(u.deg, u.deg))

	#   INDEX FOR REF.TABLE
	indx, d2d, d3d  = incoord.match_to_catalog_sky(refcoord)
	mreftbl			= reftbl[indx]
	mreftbl['sep']	= d2d
	mergetbl		= intbl
	for col in mreftbl.colnames:
		mergetbl[col]	= mreftbl[col]
	indx_sep		= np.where(mergetbl['sep']*3600.<sep)
	mtbl			= mergetbl[indx_sep]
	#mtbl.write(mergename, format='ascii', overwrite=True)
	return mtbl
#-------------------------------------------------------------------------#
def star4zp(intbl, inmagerkey, refmagkey, refmagerkey, 
	refmaglower=14., refmagupper=17., refmagerupper=0.05,
	inmagerupper=0.1, flagcut=0, verbose=False, plot=False, plotout='zphist.png'):
	"""
	SELECT STARS FOR USING ZEROPOINT CALCULATION
	INPUT   :   TABLE, IMAGE MAG.ERR KEYWORD, REF.MAG. KEYWORD, REF.MAG.ERR KEYWORD
	OUTPUT  :   NEW TABLE
	"""
	# import numpy as np

	indx_flag = np.where( intbl['FLAGS'] <= flagcut )
	indx_refmag = np.where(	(intbl[refmagkey] <= refmagupper) &
							(intbl[refmagkey] >= refmaglower))
	indx_refmager = np.where( intbl[refmagerkey] <= refmagerupper )
	indx_mager = np.where( intbl[inmagerkey] <= inmagerupper )

	indx_all = np.where((intbl['FLAGS'] <= flagcut) & 
						(intbl[refmagkey] < refmagupper) & 
						(intbl[refmagkey] > refmaglower) & 
						(intbl[refmagerkey] < refmagerupper) &
						(intbl[inmagerkey] < inmagerupper) )

	newtbl  = intbl[indx_all]
	comment = '-'*50+'\n' \
			+ 'ALL\t\t\t\t\t: {}\n'.format(len(intbl)) \
			+ '-'*50+'\n' \
			+ 'FLAG(<={})\t\t\t\t: {}\n'.format(flagcut, len(indx_flag[0])) \
			+ '{} REF. MAGCUT ({}-{})\t\t: {}\n'.format(refmagkey, refmaglower, refmagupper, len(indx_refmag[0])) \
			+ '{} REF. MAGERR CUT < {}\t\t: {}\n'.format(refmagerkey, refmagerupper, len(indx_refmager[0])) \
			+ '{} CUT < {}\t\t: {}\n'.format(inmagerkey, inmagerupper, len(indx_mager[0])) \
			+ '='*50+'\n' \
			+ 'TOTAL #\t\t\t\t\t: {}\n'.format(len(newtbl)) \
			+ '='*50
	if verbose != False:
		print(comment)
	'''
	if plot!=False:
		# import matplotlib.pyplot as plt
		plt.close('all')
		plt.rc('font', family='serif')
		fig = plt.figure()
		x = 1080 / fig.dpi
		y = 1080 / fig.dpi
		fig.set_figwidth(x)
		fig.set_figheight(y)
		#-------------------------------------------------------------
		plt.subplot(221)
		plt.title('FLAGS')

		bins = np.arange(np.min(intbl['FLAGS']), np.max(intbl['FLAGS'])+0.2, 0.2)
		plt.hist(intbl['FLAGS'], bins=bins, color='tomato', alpha=0.5,)
				# label='{}/{}'.format(len(intbl[indx_flag]), len(intbl)))
		plt.hist(intbl['FLAGS'][indx_flag], bins=bins, color='dodgerblue', alpha=0.75,
				label='{}/{}'.format(len(intbl[indx_flag]), len(intbl)))
		plt.minorticks_on()
		plt.grid(color='grey', alpha=0.5, linestyle='--')
		plt.legend()
		#-------------------------------------------------------------
		plt.subplot(222)
		plt.title('{} MAG'.format(refmagkey))
		bins = np.arange(np.min(intbl[refmagkey]), np.max(intbl[refmagkey])+0.1, 0.1)
		plt.hist(intbl[refmagkey], bins=bins, color='tomato', alpha=0.5)

		plt.hist(intbl[refmagkey][indx_refmag], bins=bins, color='dodgerblue', alpha=0.75,
				label='{}/{}'.format(len(intbl[indx_refmag]), len(intbl)))
		plt.axvline(x=refmaglower, linestyle='--', color='k',)
		plt.axvline(x=refmagupper, linestyle='--', color='k',)

		plt.minorticks_on()
		plt.grid(color='grey', alpha=0.5, linestyle='--')
		plt.legend()
		#-------------------------------------------------------------
		plt.subplot(223)
		plt.title('{}'.format(refmagerkey))
		bins = np.arange(np.min(intbl[refmagerkey]), np.max(intbl[refmagerkey])+0.001, 0.001)
		plt.hist(intbl[refmagerkey], bins=bins, color='tomato', alpha=0.5)

		plt.hist(intbl[refmagerkey][indx_refmager], bins=bins, color='dodgerblue', alpha=0.75,
				label='{}/{}'.format(len(intbl[indx_refmager]), len(intbl)))
		plt.axvline(x=refmagerupper, linestyle='--', color='k',)
		plt.minorticks_on()
		plt.grid(color='grey', alpha=0.5, linestyle='--')
		plt.legend()
		#-------------------------------------------------------------
		plt.subplot(224)
		plt.title('{}'.format(inmagerkey))
		bins = np.arange(np.min(intbl[inmagerkey]), np.max(intbl[inmagerkey])+0.01, 0.01)
		plt.hist(intbl[inmagerkey], bins=bins, color='tomato', alpha=0.5)
		plt.hist(intbl[inmagerkey][indx_mager], bins=bins, color='dodgerblue', alpha=0.75,
				label='{}/{}'.format(len(intbl[indx_mager]), len(intbl)))
		plt.axvline(x=inmagerupper, linestyle='--', color='k',)
		# plt.xlim([0.0, 1.0])
		plt.minorticks_on()
		plt.grid(color='grey', alpha=0.5, linestyle='--')
		plt.legend()
		#-------------------------------------------------------------
		plt.tight_layout()
		# plt.savefig(plotout, dpi=300)
		plt.savefig(plotout)
	'''
	return newtbl
#-------------------------------------------------------------------------#
# @jit
def zpcal(intbl, inmagkey, inmagerkey, refmagkey, refmagerkey, sigma=2.0, method='default'):
	"""
	ZERO POINT CALCULATION
	3 SIGMA CLIPPING (MEDIAN)

	import matplotlib.pyplot as plt
	from numpy import median
	import numpy as np
	"""
	from astropy.stats import sigma_clip
	#	REMOVE BLANK ROW (=99)	
	indx_avail      = np.where( (intbl[inmagkey] != 99) & (intbl[refmagkey] != 99) )
	intbl           = intbl[indx_avail]
	zplist          = np.copy(intbl[refmagkey] - intbl[inmagkey])
	intbl['zp']		= zplist
	#	SIGMA CLIPPING
	zplist_clip     = sigma_clip(zplist, sigma=sigma, maxiters=None, cenfunc=np.median, copy=False)
	# zplist_clip     = sigma_clip(zplist, sigma=sigma, maxiters=None, cenfunc=bn_median, copy=False)
	indx_alive      = np.where( zplist_clip.mask == False )
	indx_exile      = np.where( zplist_clip.mask == True )
	#	RE-DEF. ZP LIST AND INDEXING CLIPPED & NON-CLIPPED
	intbl_alive     = intbl[indx_alive]
	intbl_exile     = intbl[indx_exile]
	#	ZP & ZP ERR. CALC.
	if method == 'default':
		zp              = np.median(np.copy(intbl_alive['zp']))
		zper			= np.std(np.copy(intbl_alive['zp']))
	elif method == 'weightedmean':
		print(method)
		mager = sqsum(intbl_alive[inmagerkey], intbl_alive[refmagerkey])
		w0 = 1/mager
		w = w0/np.sum(w0)
		zp = np.sum(w*intbl_alive['zp'])/np.sum(w)
		zper = 1/np.sqrt(np.sum(w))
	return zp, zper, intbl_alive, intbl_exile
#-------------------------------------------------------------------------#
# def zpplot(outname, alltbl, otbl, xtbl, inmagkey, inmagerkey, refmagkey, refmagerkey, zp, zper, refmaglower, refmagupper):
def zpplot(outname, otbl, xtbl, inmagkey, inmagerkey, refmagkey, refmagerkey, zp, zper, refmaglower, refmagupper):
	# import numpy as np
	# import matplotlib.pyplot as plt
	#   FILE NAME
	plt.close('all')
	plt.rc('font', family='serif')
	fig = plt.figure()
	x = 1920 / 2 / fig.dpi
	y = 1080 / 2 / fig.dpi
	# x = 1920 / 4 / fig.dpi
	# y = 1080 / 4 / fig.dpi
	fig.set_figwidth(x)
	fig.set_figheight(y)

	alltbl = vstack([otbl, xtbl])
	# plt.rcParams.update({'font.size': 16})
	plt.axhline(
				zp,
				linewidth=2, linestyle='-', color='gray',
				label= r'ZP={}$\pm${}'.format(round(zp, 3), round(zper, 3))
				)
	plt.fill_between(
						# [np.min(otbl[refmagkey])-0.05, np.max(otbl[refmagkey])+0.05],
						[np.min(alltbl[refmagkey]), np.max(alltbl[refmagkey])],
						zp-zper, zp+zper,
						color='silver', alpha=0.3
					)
	# plt.errorbar(	
	# 				alltbl[refmagkey], alltbl['zp'],
	# 				yerr=tool.sqsum(alltbl[inmagerkey], alltbl[refmagerkey]),
	# 				c='silver', ms=6, marker='.', ls='',
	# 				capsize=2.5, capthick=1,
	# 				label='All stars ({})'.format(len(alltbl)), alpha=0.25,
	# 			)
	plt.errorbar(	
					otbl[refmagkey], otbl['zp'],
					yerr=tool.sqsum(otbl[inmagerkey], otbl[refmagerkey]),
					c='dodgerblue', ms=6, marker='o', ls='',
					capsize=5, capthick=1,
					label='Used stars ({})'.format(len(otbl)), alpha=0.5,
				)
	# plt.scatter( xtbl[refmagkey], xtbl['zp'], color='tomato', s=50, marker='x', linewidth=1, alpha=1.0, label='CLIPPED ({})'.format(len(xtbl)) )
	plt.errorbar(	
					xtbl[refmagkey], xtbl['zp'],
					yerr=tool.sqsum(xtbl[inmagerkey], xtbl[refmagerkey]),
					c='tomato', ms=6, marker='x', ls='',
					capsize=5, capthick=1,
					label='Clipped stars ({})'.format(len(xtbl)), alpha=0.5,
				)
	plt.axhline(
				zp+zper,
				linewidth=2, linestyle='-.', color='gray', alpha=0.5,
				label='Clip upper ({})'.format(len(xtbl[xtbl['zp']>zp+zper])),
				)
	plt.axhline(
				zp-zper,
				linewidth=2, linestyle='--', color='gray', alpha=0.5,
				label='Clip lower ({})'.format(len(xtbl[xtbl['zp']<zp-zper])),
				)
	plt.axvline(x=refmaglower, linestyle='--', color='k', alpha=0.5)
	plt.axvline(x=refmagupper, linestyle='--', color='k', alpha=0.5)
	# plt.xlim(np.min(otbl[refmagkey])-0.05, np.max(otbl[refmagkey])+0.05)
	# plt.ylim(zp-0.5, zp+0.5)
	#	SETTING
	plt.title(outname, {'fontsize': 12})
	plt.gca().invert_yaxis()
	plt.xlabel('REF.MAG.', {'color': 'black', 'fontsize': 20})
	plt.ylabel('ZP [AB]', {'color': 'black', 'fontsize': 20})
	plt.legend(loc='best', prop={'size': 14}, edgecolor=None)
	plt.tight_layout()
	plt.minorticks_on()
	plt.savefig(outname)
	#	PRINT
	print('MAG TYP     : '+inmagkey)
	print('ZP          = '+str(round(zp, 3)))
	print('ZP ERR      = '+str(round(zper, 3)))
	print('STD.NUMB    = '+str(int(len(otbl))))
	print('REJ.NUMB    = '+str(int(len(xtbl))))
	print('CLIP UPPER  = {}'.format(len(xtbl[xtbl['zp']>zp+zper])))
	print('CTIP LOWER  = {}'.format(len(xtbl[xtbl['zp']<zp-zper])))
#-------------------------------------------------------------------------#
def bkgest_mask(inim):
	'''
	'''
	import numpy as np
	from astropy.io import fits
	from photutils import make_source_mask
	from numpy import mean,median
	from astropy.stats import sigma_clipped_stats
	
	data    =fits.getdata(inim)
	# mask    = make_source_mask(data, snr=3, npixels=5, dilate_size=11)
	mask    = make_source_mask(data, nsigma=3, npixels=5, dilate_size=11)
	mean, median, std = sigma_clipped_stats(data, sigma=3.0, mask=mask)
	return mean, median, std
#-------------------------------------------------------------------------#
def limitmag(N, zp, aper, skysigma):			# 3? 5?, zp, diameter [pixel], skysigma
	import numpy as np
	R           = float(aper)/2.				# to radius
	braket      = N*skysigma*np.sqrt(np.pi*(R**2))
	upperlimit  = float(zp)-2.5*np.log10(braket)
	return round(upperlimit, 3)
#-------------------------------------------------------------------------#
'''
def targetfind(ra1, de1, ra2, de2, sep):

	import numpy as np
	dist	= np.sqrt( (ra1-ra2)**2. + (de1-de2)**2. )
	indx	= np.where( (dist == np.min(dist)) &
						(dist < sep/3600.) )
	return indx
'''
#-------------------------------------------------------------------------#
def plotshow(inim, numb_list, xim_list, yim_list, outname='default', add=None, numb_addlist=None, xim_addlist=None, yim_addlist=None, invert=False):
	'''
	PLOT IMAGE AND SHOW DESINATED OBJECTS
	'''
	import numpy as np
	import matplotlib.pyplot as plt
	from astropy.io import fits
	from matplotlib.colors import LogNorm
	from matplotlib.patches import Circle
	from astropy.visualization import (MinMaxInterval, SqrtStretch, ImageNormalize)
	from astropy.visualization import ZScaleInterval, LinearStretch
	from astropy.wcs import WCS
	if outname == 'default':
		outname		= inim[:-5]+'.png'
	else:
		pass
	data, hdr	= fits.getdata(inim, header=True)
	if invert == True:
		data = -1*data
	#fig, ax		= plt.subplots(1)

	plt.close('all')
	plt.rc('font', family='serif')
	fig = plt.figure()
	x = 1080 / 2 / fig.dpi
	y = 1080 / 2 / fig.dpi
	fig.set_figwidth(x)
	fig.set_figheight(y)

	wcs			= WCS(hdr)
	norm_zscale	= ImageNormalize(data, interval=ZScaleInterval(), stretch=LinearStretch())
	# fig			= plt.figure() 
	ax			= plt.subplot(projection=wcs)
	im			= ax.imshow(data, cmap='gray', origin='lower', norm=norm_zscale) 

	for xx, yy in zip(xim_list, yim_list):
		circ = Circle((xx, yy), 15, color='gold', fill=None, linewidth='0.3')
		ax.add_patch(circ)
	for i, txt in enumerate(numb_list):
		xim		= xim_list[i]
		yim		= yim_list[i]
		ax.text(xim+7.5, yim+7.5, str(txt), color='gold', fontsize=5)
	if add != None:
		for xx, yy in zip(xim_addlist, yim_addlist):
			'''
			circ = Circle((xx, yy), 15, color='tomato', fill=None, linewidth='0.5')
			ax.add_patch(circ)
			'''
			ax.scatter(xx, yy, color='tomato', marker='o', alpha=0.1, s=10)
		for i, txt in enumerate(numb_addlist):
			xim		= xim_addlist[i]
			yim		= yim_addlist[i]
			ax.text(xim+7.5, yim+7.5, str(txt), color='tomato', fontsize=5)
	else:
		pass
	ax.grid('both', linestyle='--', color='silver', alpha=0.5)
	ax.set_xlabel('R.A.', fontsize=20)
	ax.set_ylabel('Dec.', fontsize=20)
	# plt.tight_layout()
	# plt.minorticks_on()
	fig.savefig(outname, facecolor='w', edgecolor='w',
		orientation='portrait', papertype=None, format=None,
		transparent=False, bbox_inches=None, pad_inches=0.1,
		metadata=None)
#-------------------------------------------------------------------------#
def sedualcom(inim, gain, pixscale, seeing, det_sigma=1.5, backsize=str(64), backfiltersize=str(3), detect='detection.fits'):
	"""
	SourceEXtractor
	APERTURE    3", 5", 7",
				1.0seeing, 1.2seeing ,1.5seeing ,1.7seeing ,2.0seeeing
	INPUT   :   (image).fits
				aperture    []
				seeing_fwhm [pixel]
	OUTPUT  :   no return
				.cat
	"""
	import numpy as np
	import os
	from astropy.io import ascii
	#   FILE CHECK
	#	CONFIG FILES (USER BASE PATH)
	sharepath       = '/home/sonic/Research/yourpy/config'
	configfile      = sharepath+'/targetphot.sex'
	paramfile       = sharepath+'/targetphot.param'
	nnwfile		    = sharepath+'/targetphot.nnw'
	convfile	    = sharepath+'/targetphot.conv'
	try:
		comment = '\nSourceEXtractor (DUAL MODE) START\n' \
				+ 'IMAGE\t\t: '+inim+'\n' \
				+ 'GAIN\t\t: '+str(gain)+'\n' \
				+ 'PIXSCALE\t: '+str(pixscale)+'\n' \
				+ 'DETECTION SIGMA\t: '+str(det_sigma)+'\n' \
				+ 'PARAM\t\t: '+paramfile+'\n' \
				+ 'BACKSIZE\t: '+backsize+'\n' \
				+ 'BACKFILTER\t: '+backfiltersize+'\n' \
				+ 'CONFIG\t\t: '+configfile+'\n' \
				+ 'NNW\t\t: '+nnwfile+'\n' \
				+ 'CONVOLVE\t: '+convfile
		print(comment)
	except:
		comment = 'CHECK configfile/paramfile/nnewfile/convfile or others.'
		print(comment)
	#   FILE NAME
	oriim   = inim[2:]
	cat     = inim[:-5]+'.dual.cat'
	seg     = inim[:-5]+'.seg.fits'
	bkg     = inim[:-5]+'.bkg.fits'
	sub     = inim[:-5]+'.sub.fits'
	# psf     = inim[2:-5]+'.psf'
	psf		= inim[:-5]+'.psf'
	aper    = inim[:-5]+'.aper.fits'

	#   BASIC INFO.
	det_area        = 5
	#det_thresh      = det_sigma/np.sqrt(det_area)
	det_thresh      = det_sigma
	detecminarea    = str(det_area)
	detectthresh    = str(det_thresh)
	# seeing, fwhm_arcsec = psfex(oriim, pixscale)
	# seeing, fwhm_arcsec = psfex(inim, pixscale)
	#pixscale        = pixscalecalc(inim)
	#   OPTION
	aperture        = '%.2f'%(3./pixscale)+','+'%.2f'%(5./pixscale)+','+'%.2f'%(7./pixscale)+','+'%.2f'%(seeing)+','+'%.2f'%(1.2*seeing)+','+'%.2f'%(1.5*seeing)+','+'%.2f'%(1.7*seeing)+','+'%.2f'%(2.0*seeing)

	option1 = ' -CATALOG_NAME '+cat+' -PARAMETERS_NAME '+paramfile
	option2 = ' -DETECT_MINAREA '+detecminarea+' -DETECT_THRESH '+detectthresh \
			+' -FILTER Y '+'-FILTER_NAME '+convfile
	option3 = ' -PHOT_APERTURES '+aperture \
			+' -GAIN '+'%.1f'%(gain)+' -PIXEL_SCALE '+'%.1f'%(pixscale)
	option4 = ' -SEEING_FWHM '+'%.3f'%(seeing)+' -STARNNW_NAME '+nnwfile
	option5 = ' -BACK_SIZE '+ backsize \
			+ ' -BACK_FILTERSIZE '+ backfiltersize+' -BACKPHOTO_TYPE LOCAL'
	option6 = ' -CHECKIMAGE_TYPE SEGMENTATION,APERTURES,BACKGROUND,-BACKGROUND'
	option7 = ' -CHECKIMAGE_NAME '+seg+','+aper+','+','+bkg+','+sub
	option8 = ' -PSF_NAME '+psf
	#   COMMAND
	#   detect = detection.fits is fine image show target significantly and have good seeing and many stars i
	dualcom ='sex -c '+configfile+' '+detect+' , '+inim+' -CATALOG_NAME '+cat+' -PARAMETERS_NAME '+paramfile+ ' '+option2+' '+option3+' '+option4+' '+option5+' '+option6+' '+option7+' '+option8

	os.system(dualcom)
	secat   = ascii.read(cat)
	# return secat, cat, seeing, fwhm_arcsec
	return secat, cat
#-------------------------------------------------------------------------#
def targetfind(tra, tdec, refra, refdec, sep):
	import astropy.units as u
	from astropy.coordinates import SkyCoord
	targ_coord	= SkyCoord(tra, tdec, unit=(u.deg, u.deg))
	phot_coord	= SkyCoord(refra, refdec, unit=(u.deg, u.deg))
	indx, d2d, d3d	= targ_coord.match_to_catalog_sky(phot_coord)
	return indx, d2d, d3d
#-------------------------------------------------------------------------#
def apass2med(incat, outcat, sedcat='/home/paek/table/stellar_sed_template_phot4med.dat'):
	from astropy.modeling import models, fitting
	from astropy.io import ascii
	#============================================================
	#	Function
	#------------------------------------------------------------
	def chi_sq(obs, exp):
		'''
		obs : observation value
		exp : expectation value
		'''
		return np.sum( ((obs-exp)**2)/exp )
	#------------------------------------------------------------
	#	Path
	#------------------------------------------------------------
	# path_plot = f'{path_base}/plot'
	#------------------------------------------------------------
	#	Table
	#------------------------------------------------------------
	reftbl = ascii.read(incat)
	sedtbl = ascii.read(sedcat)
	#------------------------------------------------------------
	#	Advance Preparation
	#------------------------------------------------------------
	#	med
	# filterlist_med = ['m575', 'm625', 'm675', 'm725', 'm775']
	filterlist_med = ['m425', 'm475', 'm525', 'm575', 'm625', 'm675', 'm725', 'm775', 'm825', 'm875', 'm925', 'm975', 'm1025', 'n6780', 'n6830',]
	#	APASS
	filterlist_apass = ['B', 'V', 'g', 'r', 'i']
	#	Vega --> AB
	reftbl['Bmag'] = reftbl['Bmag']+(-0.09)
	reftbl['Vmag'] = reftbl['Vmag']+(-0.02)
	#	Generate table space on the reference catalog
	for filte in filterlist_med:
		reftbl[f'{filte}mag'] = 0.0
		reftbl[f'e_{filte}mag'] = 0.0
	reftbl['n_sed'] = 0
	reftbl['mktype'] = ' '*10
	#	Dummy column
	sedtbl['chisq'] = 0.0
	sedtbl['const'] = 0.0
	#------------------------------------------------------------
	#	Main body
	#------------------------------------------------------------
	# i = 2
	print('-'*60)
	print(f'#\tConvert APASS catalog (stars:{len(reftbl)})')
	print('-'*60)
	for i in range(len(reftbl)):
		n_apass = reftbl['NUMBER'][i]
		ms_apass = np.array([reftbl[f'{filte}mag'][i] for filte in filterlist_apass])
		e_ms_apass = np.array([reftbl[f'e_{filte}mag'][i] for filte in filterlist_apass])
		#	Variable for fitting
		y = ms_apass
		yerr = e_ms_apass
		# j = 60
		for j in range(len(sedtbl)):
			ms_sed = np.array([sedtbl[f'{filte}'][j].item() for filte in filterlist_apass])
			#	Variable for fitting
			x = ms_sed
			#	Initial parameter
			m_dif = np.median(y-x)
			#------------------------------------------------------------
			#	Fitting
			#------------------------------------------------------------
			# initialize a linear fitter
			fit = fitting.LinearLSQFitter()
			# initialize a linear model
			line_init = models.Linear1D()
			line_init.slope.fixed = True
			# fit the data with the fitter
			# fitted_line = fit(line_init, x, y, weights=e_ms_apass)
			fitted_line = fit(line_init, x, y, )
			#------------------------------------------------------------
			#	Fitting results
			#------------------------------------------------------------
			chisq = chi_sq(ms_apass, fitted_line(x))
			sedtbl['chisq'][j] = chisq
			sedtbl['const'][j] = fitted_line.intercept.value
		#------------------------------------------------------------
		#	Pick optimized index
		#------------------------------------------------------------
		indx_opt = np.where(sedtbl['chisq']==np.min(sedtbl['chisq']))
		#------------------------------------------------------------
		#	Optimized values
		#------------------------------------------------------------
		chisq_opt = sedtbl['chisq'][indx_opt].item()
		const_opt = sedtbl['const'][indx_opt].item()
		ms_sed = np.array([sedtbl[f'{filte}'][indx_opt].item() for filte in filterlist_apass])
		x = ms_sed
		'''
		#------------------------------------------------------------
		# 	Plot
		#------------------------------------------------------------
		plt.close()
		fig, (ax1, ax2) = plt.subplots(2, 1)
		#	axis 1
		ax1.plot(x, y, 'ko', label='APASS')
		ax1.plot(x, x+const_opt, 'r-', label=f'Fitted Model, chisq={round(chisq_opt, 3)}')
		ax1.legend()
		ax1.set_title(f'APASS : {i}')

		ax1.set_ylabel(r'$\rm m_{APASS}$ [mag]')
		#	axis 2
		ax2.plot(x, y - (x+const_opt), 'ko')
		ax2.axhline(y=0, ls='--', color='grey')
		ax2.set_xlabel(r'$\rm m_{spec}$ [mag]')
		ax2.set_ylabel('Residual')
		ax2.set_ylim([-0.5, +0.5])

		plt.savefig(f'{path_plot}/plot_{i}.png')#, overwrite=True)

		#	Put med-band photometries on the APASS catalog
		for filte in filterlist_med:
			reftbl[f'{filte}mag'][i] = sedtbl[f'{filte}'][indx_opt].item()+const_opt
		reftbl['n_sed'][i] = sedtbl['number'][indx_opt].item()
		reftbl['mktype'][i] = sedtbl['mktype'][indx_opt].item()
		'''
		#	Put med-band photometries on the APASS catalog
		for filte in filterlist_med:
			reftbl[f'{filte}mag'][i] = sedtbl[f'{filte}'][indx_opt].item()+const_opt
		reftbl['n_sed'][i] = sedtbl['number'][indx_opt].item()
		reftbl['mktype'][i] = sedtbl['mktype'][indx_opt].item()
	
	reftbl.write(f'{outcat}', format='ascii.tab', overwrite=True)
	return reftbl

