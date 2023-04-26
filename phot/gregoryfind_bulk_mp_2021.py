#!/home/paek/anaconda3/bin/python3.7
#============================================================
#   Transient Finder
#	
#	21.05.27	Created by Gregory S.H. Paek
#============================================================
import os, glob, sys
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
# from datetime import date
from astropy.visualization.wcsaxes import SphericalCircle
import warnings
from astroquery.imcce import Skybot
warnings.filterwarnings("ignore")
# from itertools import product
from itertools import repeat
import multiprocessing
from matplotlib.patches import Circle, PathPatch
from astropy.visualization import SqrtStretch, LinearStretch, LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import ZScaleInterval, MinMaxInterval
#------------------------------------------------------------
#	Function
#------------------------------------------------------------
def inverse(inim, outim):
	data, hdr = fits.getdata(inim, header=True)
	# w = WCS(inim)
	invdata = data*(-1)
	fits.writeto(outim, invdata, header=hdr, overwrite=True)
#------------------------------------------------------------
def routine_se(inim, outcat, aperpix, seeing, conf_sex, conf_param, conf_conv, conf_nnw):
	param_insex = dict(
						#------------------------------
						#	CATALOG
						#------------------------------
						CATALOG_NAME = outcat,
						#------------------------------
						#	CONFIG FILES
						#------------------------------
						CONF_NAME = conf_sex,
						PARAMETERS_NAME = conf_param,
						FILTER_NAME = conf_conv,    
						STARNNW_NAME = conf_nnw,
						#------------------------------
						#	PHOTOMETRY
						#------------------------------
						#	DIAMETER
						#	OPT.APER, (SEEING x2), x3, x4, x5
						#	MAG_APER	OPT.APER
						#	MAG_APER_1	OPT.GAUSSIAN.APER
						#	MAG_APER_2	SEEINGx2
						#	...
						PHOT_APERTURES = str(aperpix),
						SATUR_LEVEL  = '65000.0',
						# GAIN = str(gain.value),
						# PIXEL_SCALE = str(pixscale.value),
						#------------------------------
						#	STAR/GALAXY SEPARATION
						#------------------------------
						SEEING_FWHM = str(seeing),
						)
	com = gpphot.sexcom(inim, param_insex)
	print(com)
	os.system(com)
#------------------------------------------------------------
def get_mad(data):
	return np.median(np.absolute(data - np.median(data, axis=0)), axis=0)
#------------------------------------------------------------
def plot_snapshot(data, wcs, peeing, outpng, save=True):
	plt.close('all')
	plt.rc('font', family='serif')
	fig = plt.figure(figsize=(1, 1))
	fig.set_size_inches(1. * data.shape[0] / data.shape[1], 1, forward = False)
	x = 720 / fig.dpi
	y = 720 / fig.dpi
	fig.set_figwidth(x)
	fig.set_figheight(y)
	#	No axes
	# ax = plt.subplot(projection=wcs)
	ax = plt.Axes(fig, [0., 0., 1., 1.])
	ax.set_axis_off()
	fig.add_axes(ax)

	from astropy.visualization.stretch import LinearStretch
	#	Sci
	data[np.isnan(data)] = 0.0
	transform = LinearStretch()+ZScaleInterval()
	bdata = transform(data)
	# pylab.subplot(131)
	ax.imshow(bdata, cmap="gray", origin="lower")

	#	Circle
	circle = Circle(
		(data.shape[0]/2., data.shape[1]/2.),
		2*peeing,
		edgecolor='yellow',
		lw=3,
		facecolor=None,
		fill=False
	)

	ax.add_patch(circle)

	#	RA, Dec direction
	ra0, dec0 = wcs.all_pix2world(0, 0, 1)
	ra1, dec1 = wcs.all_pix2world(data.shape[0], data.shape[1], 1)
	if ra0>ra1:
		pass
	elif ra0<ra1:
		ax.invert_xaxis()
	if dec0>dec1:
		ax.invert_yaxis()
	elif dec0<dec1:
		pass
	if save:
		plt.savefig(outpng, dpi=100)#, overwrite=True)
	else:
		pass
#------------------------------------------------------------
from astropy.nddata import Cutout2D
def snapshot(tctbl, i, cutsize=2.0):
	#	Images
	n = tctbl['NUMBER'][i].item()
	inim, hcim, hdim = tctbl['inim'][i], tctbl['hcim'][i], tctbl['hdim'][i]
	#	Poistion of transient candidate
	tra = tctbl['ALPHA_J2000'][i].item()
	tdec = tctbl['DELTA_J2000'][i].item()
	position = SkyCoord(tra, tdec, frame='icrs', unit='deg')
	#	Seeing
	seeing = tctbl['seeing'][i].item()
	# peeing = tctbl['peeing'][i].item()
	size = u.Quantity((cutsize, cutsize), u.arcmin)

	for image, kind in zip([inim, hcim, hdim], ['new', 'ref', 'sub']):
		hdu = fits.open(image)[0]
		wcs = WCS(hdu.header)
		peeing = hdu.header['PEEING']
		# Make the cutout, including the WCS
		cutout = Cutout2D(
			data=hdu.data,
			position=position,
			size=size,
			wcs=wcs,
			#	pad setting
			mode='partial',
			fill_value=float(hdu.header['SKYVAL']),
			)
		data = 	cutout.data
		# Put the cutout image in the FITS HDU
		hdu.data = cutout.data
		# Update the FITS header with the cutout WCS
		hdu.header.update(cutout.wcs.to_header())
		# Write the cutout to a new FITS file
		outim = f'{os.path.splitext(hdim)[0]}.{n}.{kind}{os.path.splitext(hdim)[1]}'
		outpng = f'{os.path.splitext(hdim)[0]}.{n}.{kind}.png'
		#	Save postage stamp *.png & *.fits
		hdu.writeto(outim, overwrite=True)
		plot_snapshot(data, wcs, peeing, outpng, save=True)
#------------------------------------------------------------
# from astropy.nddata import Cutout2D
# def snapshot(tctbl, i, cutsize=2.0):
# 	#	Images
# 	n = tctbl['NUMBER'][i].item()
# 	inim, hcim, hdim = tctbl['inim'][i], tctbl['hcim'][i], tctbl['hdim'][i]
# 	#	Poistion of transient candidate
# 	tra = tctbl['ALPHA_J2000'][i].item()
# 	tdec = tctbl['DELTA_J2000'][i].item()
# 	position = SkyCoord(tra, tdec, frame='icrs', unit='deg')
# 	#	Seeing
# 	seeing = tctbl['seeing'][i].item()
# 	# peeing = tctbl['peeing'][i].item()
# 	size = u.Quantity((cutsize, cutsize), u.arcmin)
# 	size_small = u.Quantity((seeing, seeing), u.arcsec)

# 	for image in [inim, hcim, hdim]:
# 		hdu = fits.open(image)[0]
# 		wcs = WCS(hdu.header)
# 		peeing = hdu.header['PEEING']
# 		# Make the cutout, including the WCS
# 		cutout = Cutout2D(hdu.data, position=position, size=size, wcs=wcs)
# 		scutout = Cutout2D(hdu.data, position=position, size=size_small, wcs=wcs)
# 		data = 	cutout.data
# 		sdata = scutout.data
# 		# Put the cutout image in the FITS HDU
# 		hdu.data = cutout.data
# 		# Update the FITS header with the cutout WCS
# 		hdu.header.update(cutout.wcs.to_header())
# 		#	Ref, Subt, Sci images
# 		if 'hcCalib' in image:
# 			kind = 'ref'
# 			'''
# 			norm = ImageNormalize(	
# 				vmin=np.median(data),
# 				vmax=np.median(data)*1e2,
# 				stretch=LogStretch(),
# 				interval=ZScaleInterval(),
# 				)
# 			'''
# 			#	Use scale & limit of sci image 
# 			norm = norm_save
			
# 		elif 'hdCalib' in image:
# 			kind = 'sub'
# 			norm = ImageNormalize(	
# 				vmin=np.median(data),
# 				vmax=np.max(sdata),
# 				stretch=LogStretch(),
# 				interval=ZScaleInterval(),
# 				)	
# 		else:
# 			kind = 'new'
# 			norm = ImageNormalize(	
# 				vmin=np.median(data),
# 				vmax=np.max(sdata),
# 				stretch=LogStretch(),
# 				interval=ZScaleInterval(),
# 				)		
# 			norm_save = norm	
		
# 		# Write the cutout to a new FITS file
# 		# import os 
# 		outim = f'{os.path.splitext(hdim)[0]}.{n}.{kind}{os.path.splitext(hdim)[1]}'
# 		outpng = f'{os.path.splitext(hdim)[0]}.{n}.{kind}.png'
# 		hdu.writeto(outim, overwrite=True)
# 		#------------------------------------------------------------
# 		#	Plot
# 		#------------------------------------------------------------
# 		plt.close('all')
# 		plt.rc('font', family='serif')
# 		fig = plt.figure(figsize=(1, 1))
# 		fig.set_size_inches(1. * data.shape[0] / data.shape[1], 1, forward = False)
# 		x = 720 / fig.dpi
# 		y = 720 / fig.dpi
# 		fig.set_figwidth(x)
# 		fig.set_figheight(y)
# 		#	No axes
# 		ax = plt.Axes(fig, [0., 0., 1., 1.])
# 		ax.set_axis_off()
# 		fig.add_axes(ax)

# 		ax.imshow(
# 			data, 
# 			origin='lower',
# 			norm=norm,
# 			cmap='gist_gray'
# 			)
# 		#	Circle
# 		# from matplotlib.patches import Circle, PathPatch
# 		circle = Circle(
# 			(data.shape[0]/2., data.shape[1]/2.),
# 			2*peeing,
# 			edgecolor='yellow',
# 			lw=3,
# 			facecolor=None,
# 			fill=False
# 		)
# 		ax.add_patch(circle)
# 		#	x, y range
# 		xl, xr = ax.set_xlim()
# 		yl, yu = ax.set_ylim()
# 		ax.set_xlim([xr, xl])
# 		ax.set_ylim([yu, yl])
# 		plt.savefig(outpng, dpi=200, overwrite=True)
#------------------------------------------------------------
def findloc(inim):
	path_table = "/home/paek/qsopy/tables"
	loctbl = ascii.read(f'{path_table}/obs.location.smnet.dat')
	#	code == 500 : default
	code = 500
	for i, obs in enumerate(loctbl['obs']):
		if obs in inim:
			code = loctbl['Code'][i].item()
			break
	return code
#------------------------------------------------------------
def routine(tstbl, i, cutsize=0.5):
	hdim = tstbl['hdim'][i]
	hdcat = tstbl['hdcat'][i]
	hcim = tstbl['hcim'][i]
	refim = tstbl['refim'][i]
	scicat = tstbl['scicat'][i]
	fovval = float(tstbl['fovval'][i])

	print(f"{'-'*60}\n#\tInput\n{'-'*60}")
	print(f'hdim : {os.path.basename(hdim)}')	
	print(f'hdcat : {os.path.basename(hdcat)}')	
	print(f'hcim : {os.path.basename(hcim)}')	
	print(f'scicat : {os.path.basename(scicat)}')	
	print(f'fovval : {fovval}')	
	print('-'*60)

	#	Image information
	inim = hdim.replace('hd', '')
	part = hdim.split('-')
	#	Header
	data, hdr = fits.getdata(hdim, header=True)
	w = WCS(hdim)
	filte = hdr['FILTER']
	#	Positional information
	xcent, ycent = data.shape
	c_cent = w.pixel_to_world(xcent, ycent)
	seeing = hdr['SEEING']
	aperpix = hdr['APERPIX']
	epoch = Time(hdr['DATE-OBS'], format='isot')
	#	Name
	invhdim = hdim.replace('hd', 'invhd')
	invhcim = hcim.replace('hc', 'invhc')
	invhdcat = invhdim.replace('fits', 'cat')
	invhccat = invhcim.replace('fits', 'cat')
	#	Inverse images
	inverse(hdim, invhdim)
	inverse(hcim, invhcim)
	#	SEtractor for inverse images
	routine_se(invhdim, invhdcat, aperpix, seeing, conf_sex, conf_param, conf_conv, conf_nnw)
	routine_se(invhcim, invhccat, aperpix, seeing, conf_sex, conf_param, conf_conv, conf_nnw)
	#	Catalog --> Table
	hdtbl = ascii.read(hdcat)
	hdtbl['inim'] = inim
	hdtbl['hcim'] = hcim
	hdtbl['hdim'] = hdim
	hdtbl['snr'] = 1/hdtbl['magerr_aper']
	hdtbl['seeing'] = seeing
	hdtbl['ratio_seeing'] = hdtbl['FWHM_WORLD']/hdtbl['seeing']
	hdtbl['ratio_ellipticity'] = hdtbl['ELLIPTICITY']/hdr['ELLIP']
	hdtbl['ratio_elongation'] = hdtbl['ELONGATION']/hdr['ELONG']
	#	Coordinate
	c_hd = SkyCoord(hdtbl['ALPHA_J2000'], hdtbl['DELTA_J2000'], unit=u.deg)
	#------------------------------------------------------------
	#	Flagging
	#------------------------------------------------------------
	#	Generate blank flag
	numbers = np.arange(0, 12+1, 1)	#	flag 0-11
	for num in numbers: hdtbl[f'flag_{num}'] = False
	hdtbl['flag'] = False
	#------------------------------------------------------------
	#	flag 0
	#------------------------------------------------------------
	#	Skybot query
	try:
		code = findloc(inim)
		
		sbtbl = Skybot.cone_search(c_cent, fovval*u.arcmin, epoch, location=code)
		c_sb = SkyCoord(sbtbl['RA'], sbtbl['DEC'])
		sbtbl['sep'] = c_cent.separation(c_sb).to(u.arcmin)
		#	Skybot matching
		indx_sb, sep_sb, _ = c_hd.match_to_catalog_sky(c_sb)
		hdtbl['flag_0'][
			# (sep_sb.arcsec<seeing*10)
			(sep_sb.arcsec<5)
			] = True

		#	Table combine
		"""
		from astropy.table import vstack, hstack
		vstack([hstack(hdtbl[(sep_sb.arcsec<5)], sbtbl[indx_sb][(sep_sb.arcsec<5)]), hdtbl[(sep_sb.arcsec>=5)]])
		"""

	except:
		print(f'No solar system object was found in the requested FOV ({fovval} arcmin)')
		pass
	#------------------------------------------------------------
	#	flag 1
	#------------------------------------------------------------
	invhdtbl = ascii.read(invhdcat)
	if len(invhdtbl)>0:
		#	Coordinate
		c_invhd = SkyCoord(invhdtbl['ALPHA_J2000'], invhdtbl['DELTA_J2000'], unit=u.deg)
		#	Matching with inverted images
		indx_invhd, sep_invhd, _ = c_hd.match_to_catalog_sky(c_invhd)
		hdtbl['flag_1'][
			(sep_invhd.arcsec<seeing)
			] = True
	else:
		print('Inverted subtraction image has no source. ==> pass flag1')
		pass
	#------------------------------------------------------------
	#	flag 2
	#------------------------------------------------------------
	invhctbl = ascii.read(invhccat)
	if len(invhctbl)>0:
		#	Coordinate
		c_invhc = SkyCoord(invhctbl['ALPHA_J2000'], invhctbl['DELTA_J2000'], unit=u.deg)
		#	Matching with inverted images
		indx_invhc, sep_invhc, _ = c_hd.match_to_catalog_sky(c_invhc)
		hdtbl['flag_2'][
			(sep_invhc.arcsec<seeing)
			] = True
	else:
		print('Inverted reference image has no source. ==> pass flag2')
		pass
	#------------------------------------------------------------
	#	SEtractor criterion
	#------------------------------------------------------------
	#	flag 3
	#------------------------------------------------------------
	#	Sources @edge
	frac = 0.9
	hdtbl['flag_3'][
		((hdtbl['X_IMAGE']<xcent-xcent*frac) |
		(hdtbl['X_IMAGE']>xcent+xcent*frac) |
		(hdtbl['Y_IMAGE']<ycent-ycent*frac) |
		(hdtbl['Y_IMAGE']>ycent+ycent*frac))
		] = True
	#------------------------------------------------------------
	#	flag 4
	#------------------------------------------------------------
	#	More than 5 sigma signal
	hdtbl['flag_4'][
		(hdtbl['mag_aper']>hdr['ul5_1'])
		# ((hdtbl['magerr_aper_1']>0.2) |
		# (hdtbl['mag_aper_1']>hdr['ul5_1']) |
		# (hdtbl['mag_aper_2']>hdr['ul5_1']) |
		# (hdtbl['mag_aper_3']>hdr['ul5_1']) |
		# (hdtbl['mag_aper_4']>hdr['ul5_1']) |
		# (hdtbl['mag_aper_5']>hdr['ul5_1']))
		] = True
	#	Empirical criterion
	#------------------------------------------------------------
	#	flag 5
	#------------------------------------------------------------
	hdtbl['flag_5'][
		(hdtbl['ratio_ellipticity'] > 5)
		] = True
	#------------------------------------------------------------
	#	flag 6
	#------------------------------------------------------------
	hdtbl['flag_6'][
		(hdtbl['FLAGS'] > 4.0)
		] = True
	#------------------------------------------------------------
	#	flag 7
	#------------------------------------------------------------
	hdtbl['flag_7'][
		(hdtbl['FWHM_WORLD']>seeing*3.0) |
		(hdtbl['FWHM_WORLD']<seeing*0.5)
		] = True
	#------------------------------------------------------------
	#	flag 8
	#------------------------------------------------------------
	hdtbl['flag_8'][
		(hdtbl['BACKGROUND']<-50) |
		(hdtbl['BACKGROUND']>+50)
		] = True
	#------------------------------------------------------------
	#	flag 9
	#------------------------------------------------------------
	scitbl = ascii.read(scicat)
	scitbl = scitbl[
		(scitbl['FLAGS']==0) &
		(scitbl['CLASS_STAR']>0.5)
	]

	aperdict = {
		'mag_aper':'SNR_curve',
		'mag_aper_1':'Best_Aperture',
		'mag_aper_2':'2seeing',
		'mag_aper_3':'3seeing',
		'mag_aper_4':'3arcsec',
		'mag_aper_5':'5arcsec',	
	}
		
	key0 = 'mag_aper_1'
	key1 = 'mag_aper_3'
	#	Sci. sources magnitude diff.
	indelm = scitbl[key0] - scitbl[key1]
	#	Subt. sources magnitude diff.
	hddelm = hdtbl[key0] - hdtbl[key1]
	hdtbl['del_mag'] = hddelm
	#	MED & MAD
	indelm_med = np.median(indelm)
	indelm_mad = get_mad(indelm)
	hdtbl['del_mag_med'] = indelm_med
	hdtbl['del_mag_mad'] = indelm_mad
	hdtbl['N_del_mag_mad'] = np.abs((hdtbl['del_mag']-hdtbl['del_mag_med'])/hdtbl['del_mag_mad'])
	#	out
	n = 10
	indx_out = np.where(
		(hddelm<indelm_med-indelm_mad*n) |
		(hddelm>indelm_med+indelm_mad*n)
		)
	hdtbl['flag_9'][indx_out] = True
	#------------------------------------------------------------
	#	flag 10+11
	#------------------------------------------------------------
	peeing = hdr['PEEING']
	skysig = hdr['SKYSIG']

	nbadlist = []
	ratiobadlist = []
	nnulllist = []
	if 'KCT' in inim:
		f = 0.05	# Set tighter criterion 
	else:
		f = 0.3
	for i in range(len(hdtbl)):
		tx, ty = hdtbl['X_IMAGE'][i], hdtbl['Y_IMAGE'][i]
		bkg = hdtbl['BACKGROUND'][i]
		#	Snapshot
		tsize = peeing
		y0, y1 = int(ty-tsize), int(ty+tsize)
		x0, x1 = int(tx-tsize), int(tx+tsize)
		cdata = data[y0:y1, x0:x1]
		# plt.close()
		# plt.imshow(cdata)
		crt = bkg - skysig
		cutline = cdata.size*f
		nbad = len(cdata[cdata<crt])
		try:
			ratiobad = nbad/cdata.size
		except:
			ratiobad = -99.0
		nnull = len(np.where(cdata == 1e-30)[0])
		#	Dipole
		if nbad > cutline:
			hdtbl['flag_10'][i] = True
		#	HOTPANTS Null value
		if nnull != 0:
			hdtbl['flag_11'][i] = True
		nbadlist.append(nbad)
		ratiobadlist.append(ratiobad)
		nnulllist.append(nnull)

	hdtbl['n_bad'] = nbadlist
	hdtbl['ratio_bad'] = ratiobadlist
	hdtbl['n_null'] = nnulllist
	#------------------------------------------------------------
	#	flag 12
	#------------------------------------------------------------
	# x, y = fits.getdata(hcim).shape
	# w_ref = WCS(hcim)
	# xim, yim = w_ref.world_to_pixel(c_hd)
	# indx_nosci = np.where(
	# 	(xim < 0) | (xim > x) | (yim < 0) | (yim > y)
	# )
	# hdtbl['flag_12'][indx_nosci] = True
	x, y = fits.getdata(refim).shape
	w_ref = WCS(refim)
	xim, yim = w_ref.world_to_pixel(c_hd)
	hdtbl['x_refim'] = xim
	hdtbl['y_refim'] = yim
	indx_nosci = np.where(
		(xim < 0) | (xim > x) | (yim < 0) | (yim > y)
	)
	hdtbl['flag_12'][indx_nosci] = True
	#------------------------------------------------------------
	#	Final flag
	#------------------------------------------------------------
	flag = hdtbl['flag']
	n_all = len(hdtbl)
	for n in numbers:
		tmptbl = hdtbl[hdtbl[f'flag_{n}']==True] 
		print(f'flag=={n} : {len(tmptbl)} {int(100*len(tmptbl)/n_all)}%')
		flag = flag + hdtbl[f'flag_{n}']
	hdtbl['flag'] = flag
	indx_sb = np.where(hdtbl['flag_0']==True)
	hdtbl['flag'][indx_sb] = False

	outcat = hdcat.replace('phot_sub.cat', 'transients.cat')
	hdtbl.write(outcat, format='ascii.tab', overwrite=True)
	print('-'*60)
	bgstbl = hdtbl[hdtbl[f'flag']==True]
	tctbl = hdtbl[hdtbl[f'flag']==False]
	print(f'Filtered sources\t: {len(bgstbl)} {int(100*len(bgstbl)/n_all)}%')
	print(f'Transient Candidates\t: {len(tctbl)} {int(100*len(tctbl)/n_all)}%')
	#------------------------------------------------------------
	#	Snapshot maker
	#------------------------------------------------------------
	print(f"#\tSnapshot maker ({len(tctbl)})")
	# cutsize=0.5
	if len(tctbl) > 0:
		for i in range(len(tctbl)):
			snapshot(tctbl, i, cutsize)
	else:
		print('No transient candidates.')
#------------------------------------------------------------
#	Path & Configuration
path_config = '/home/paek/config'
conv, nnw, param, sex = 'transient.conv', 'transient.nnw', 'transient.param', 'transient.sex'
conf_sex = '{}/{}'.format(path_config, sex)
conf_param = '{}/{}'.format(path_config, param)
conf_nnw = '{}/{}'.format(path_config, nnw)
conf_conv = '{}/{}'.format(path_config, conv)
#	Input
out_tstbl = sys.argv[1]
ncores = int(sys.argv[2])

tstbl = ascii.read(out_tstbl)
# cutsize=0.5
# cutsize = 1.0
cutsize = 2.0

print(f"{'='*60}\n#\tTransient Search Process\n{'-'*60}")
print(f'{len(tstbl)} images with {ncores} cores')
print('='*60)

#	Multi-processing
with multiprocessing.Pool(processes=ncores) as pool:
	results = pool.starmap(routine, zip(repeat(tstbl), range(len(tstbl)), repeat(cutsize)))
