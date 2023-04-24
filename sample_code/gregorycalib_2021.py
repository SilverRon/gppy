#============================================================
#   CALIBTRATION ROUTINE FOR IMSNG TELESCOPES
#	LOAO, DOAO, SOAO, CBNUO
#	20.10.11	Created by Gregory S.H. Paek
#============================================================
from __future__ import print_function, division, absolute_import
# from timeit import default_timer as timer
# from numba import jit
# from pyraf import iraf
import os, sys, glob, subprocess
import numpy as np
import astropy.io.ascii as ascii
import matplotlib.pyplot as plt
plt.ioff()
from astropy.nddata import CCDData
from imsng import calib
from imsng import tool_tbd
from astropy.io import fits
from astropy.table import Table, vstack
from astropy import units as u
from ccdproc import ImageFileCollection
import warnings
warnings.filterwarnings(action='ignore')
# from itertools import product
from itertools import repeat
import multiprocessing
import time
start_localtime = time.strftime('%Y-%m-%d %H:%M:%S (%Z)', time.localtime())
#============================================================
#	USER SETTING
#============================================================
try:
	obs = (sys.argv[1]).upper()
except:
	# obs = (input('OBSERVATORY? (LOAO/DOAO/SOAO/CBNUO)\t: ')).upper()
	obs = input('''Observatory(_ccd) to run
--------------------
LOAO
DOAO
SOAO
CBNUO
KCT_ASI1600MM
KCT_STX16803
RASA36
KHAO
MDFTS
SAO_C361K
---------------------
:''').upper()
print('# Observatory : {}'.format(obs.upper()))
try:
	ncores = int(sys.argv[2])
except:
	ncores = 8
# ncores = 9
# ncores = 8
# ncores = 4
# ncores = 2
#------------------------------------------------------------
#	PATH
#------------------------------------------------------------
path_base = '/data3/paek/factory'
path_gal = '/data6/IMSNG/IMSNGgalaxies'
path_refcat = '/data3/paek/factory/ref_frames/LOAO'
#------------------------------------------------------------
path_config = '/home/paek/config'
path_default_gphot = '/home/paek/config/gphot.{}.config'.format(obs.lower())
path_mframe = '{}/master_frames'.format(path_base)
path_calib = '{}/calib'.format(path_base)
#------------------------------------------------------------
#	Codes
path_phot_sg = '/home/paek/qsopy/phot/gregoryphot_2021.py'
path_phot_mp = '/home/paek/qsopy/phot/gregoryphot_mp_2021.py'
path_phot_sub = '/home/paek/qsopy/phot/gregoryphot_sub_2021.py'
# path_find = '/home/paek/qsopy/phot/gregoryfind_2021.py'
# path_find = '/home/paek/qsopy/phot/gregoryfind_bulk_2021.py'
path_find = '/home/paek/qsopy/phot/gregoryfind_bulk_mp_2021.py'
#------------------------------------------------------------
# if obs == 'DOAO':
# 	path_raw = '{}/{}'.format(path_base, obs.upper())
# elif 'KCT' in obs:
if 'KCT_ASI1600MM' in obs:
	path_raw = '/data3/IMSNG/KCT/obsdata'
else:
	path_raw = '/data6/obsdata/{}'.format(obs.upper())
rawlist = sorted(glob.glob(path_raw+'/2*'))
#------------------------------------------------------------
path_ref = '{}/ref_frames/{}'.format(path_base, obs.upper())
path_factory = '{}/{}'.format(path_base, obs.lower())
# path_save = '/data3/IMSNG/{}'.format(obs.upper())
path_save = f'/data6/bkgdata/{obs.upper()}'
path_log = '/home/paek/log/{}.log'.format(obs.lower())
#------------------------------------------------------------
path_obs = '/home/paek/table/obs.dat'
path_changehdr = '/home/paek/table/changehdr.dat'
path_alltarget = '/home/paek/table/alltarget.dat'
ccdinfo = tool_tbd.getccdinfo(obs, path_obs)
#------------------------------------------------------------
logtbl = ascii.read(path_log)
datalist = np.copy(logtbl['date'])
obstbl = ascii.read(path_obs)
hdrtbl = ascii.read(path_changehdr)
alltbl = ascii.read(path_alltarget)
#============================================================
#	MAIN BODY
#============================================================
newlist = [i for i in rawlist if (i not in datalist) & (i+'/' not in datalist)]
if len(newlist) == 0:
	print('No new data')
	sys.exit()
else:
	print(newlist)
# newlist = ['/data3/paek/factory/loao/2020_1215']
# path = newlist[0]
# path = newlist[-1]
# path = '/data6/obsdata/RASA36/20210509'
# for path in newlist:
path = newlist[0]
tdict = dict()
starttime = time.time()
path_data = '{}/{}'.format(path_factory, os.path.basename(path))
#	Remove old folder and re-copy folder
rmcom = 'rm -rf {}'.format(path_data)
print(rmcom)
os.system(rmcom)
cpcom = 'cp -r {} {}'.format(path, path_data)
print(cpcom)
os.system(cpcom)
#------------------------------------------------------------
#	Process summary status
#------------------------------------------------------------
protbl = Table()
protbl['process'] = ['master_frame', 'pre_process', 'astrometry', 'cr_removal', 'defringe', 'photometry', 'image_stack', 'photometry_com', 'subtraction', 'photometry_sub', 'transient_search', 'total']
protbl['status'] = False
protbl['time'] = 0.0 * u.second
#------------------------------------------------------------
#	CCD TYPE
#------------------------------------------------------------
ic0 = ImageFileCollection(path_data, keywords='*')
ic0.summary.write('{}/hdr.raw.dat'.format(path_data), format='ascii.tab', overwrite=True) 
#	Exceptions
#	DOAO CCD
if obs == 'DOAO':
	instrume = ic0.summary['instrume'][0]
	if instrume == 'Apogee USB/Net':obs = 'DOAO_APOGEE'
	elif instrume == '':			obs = 'DOAO_FLI'
	elif instrume == 'FLI':			obs = 'DOAO_FLI'
	elif instrume == 'Moravian C4-16000':		obs = 'DOAO_C416K'
else:
	pass
obsinfo = calib.getobsinfo(obs, obstbl)
if obs == 'SAO_C361K':
	xbinning = ic0.summary['xbinning'][0]
	if xbinning > 1:
		print(f'{obs} : BINNINGx{xbinning}')
		obsinfo['pixscale'] = obsinfo['pixscale']*xbinning
	if ic0.summary['instrume'][0] == 'ASCOM Camera Driver for FLI Kepler':
		obsinfo['pixscale'] = 0.311
		obsinfo['fov'] = 0.25*60.
		obsinfo['gain'] = float(ic0.summary['egain'][0])
#	KHAO Binning
if obs == 'KHAO':
	xbinning = ic0.summary['xbinning'][0]
	if xbinning > 1:
		obsinfo['pixscale'] = obsinfo['pixscale']*2
#	LOAO
if obs == 'LOAO':
	instrume = ic0.summary['instrume'][0]
	if 'Finger Lakes' in instrume:
		obs = 'LOAO_FLI'
		obsinfo = calib.getobsinfo(obs, obstbl)
#	RASA Mode
if obs == 'RASA36':
	if 'hdr' in path:
		mode = 'hdr'
		badmode = True
	elif 'high' in path:
		mode = 'high'
		badmode = True
	else:
		mode = 'high'
		badmode = True
		pass
		# biasim = ic0.summary['file'][ic0.summary['imagetyp']=='Bias Frame'][0]
		# biaslevel = np.median(fits.getdata(f'{path_data}/{biasim}').flatten())
		# if biaslevel > 87:
		# 	mode = 'hdr'
		# else:
		# 	mode = 'high'
		# badmode = False
	# print(f'RASA36 [{mode} mode (Bad mode : {badmode})]')# : bias level = {biaslevel}')
	print(f'Master Frame Mode:{mode} [Bad Mode:{badmode}]')


calib.correcthdr_routine(path_data, hdrtbl, obs)
# objfilterlist, objexptimelist, flatfilterlist, darkexptimelist, obstime = calib.correcthdr_routine(path_data, hdrtbl)
ic1 = ImageFileCollection(path_data, keywords='*')
ic1.summary.write('{}/hdr.cor.dat'.format(path_data), format='ascii.tab', overwrite=True)
try:
	nobj = len(ic1.filter(imagetyp='OBJECT').summary)
except:
	nobj = len(ic1.filter(imagetyp='object').summary)	
#------------------------------------------------------------
#	Slack messages
#------------------------------------------------------------	
path_keys = '/home/paek/table'
keytbl = ascii.read(f'{path_keys}/keys.dat')
OAuth_Token = keytbl['key'][keytbl['name']=='slack'].item()

channel = '#pipeline'
text = f'[Pipeline/{obs}] Start Processing {os.path.basename(path)} Data ({nobj} objects) with {ncores} cores'

param_slack = dict(
	token = OAuth_Token,
	channel = channel,
	text = text,
)

tool_tbd.slack_bot(**param_slack)
#============================================================
#	Making master frames (BIAS, DARK, FLAT)
#============================================================
st = time.time()
#------------------------------------------------------------
#	BIAS
#------------------------------------------------------------
try:
	biasnumb = len(ic1.filter(imagetyp='Bias').summary)
except:
	biasnumb = 0
# if len(ic1.filter(imagetyp='Bias').summary) != 0:
if biasnumb != 0:
	mzero = calib.master_zero(ic1, fig=False)
	# print(zeroim)
	date = fits.getheader(f'{path_data}/zero.fits')['date-obs'][:10].replace('-', '')
	if obs == 'RASA36':
		zeroim = f'{path_mframe}/{obs}/zero/{date}-zero_{mode}.fits'
	else:
		zeroim = f'{path_mframe}/{obs}/zero/{date}-zero.fits'
	cpcom = f'cp {path_data}/zero.fits {zeroim}'
	print(cpcom)
	os.system(cpcom)
	plt.close('all')
else:
	#	IF THERE IS NO FLAT FRAMES, BORROW FROM CLOSEST OTHER DATE
	print('\nNO BIAS FRAMES\n')
	if obs == 'RASA36':
		pastzero = np.array(glob.glob(f'{path_mframe}/{obs}/zero/*zero_{mode}.fits'))
	else:
		pastzero = np.array(glob.glob(f'{path_mframe}/{obs}/zero/*zero.fits'))
	#	CALCULATE CLOSEST ONE FROM TIME DIFFERENCE
	deltime = []
	for date in pastzero:
		zeromjd = calib.isot_to_mjd((os.path.basename(date)).split('-')[0])
		deltime.append(np.abs(ic1.summary['mjd'][0]-zeromjd))
	indx_closet = np.where(deltime == np.min(deltime))
	tmpzero = path_data+'/'+os.path.basename(np.asscalar(pastzero[indx_closet]))
	cpcom = 'cp {} {}'.format(np.asscalar(pastzero[indx_closet]), tmpzero)
	print(cpcom)
	os.system(cpcom)
	# if obs != 'KCT':
	if 'KCT_ASI1600MM' in obs:
		#KCT Exception
		mzero = CCDData.read(tmpzero, hdu=0, unit='adu')
	elif obs == 'RASA36':
		if (mode == 'high') & (badmode == True):
			mzero = CCDData.read(tmpzero, hdu=0).multiply(20)
			print('[Bad mode] Multiply 20 on the high mode bias.')
		else:
			mzero = CCDData.read(tmpzero, hdu=0)
	elif obs == 'LSGT':
		mzero = CCDData.read(tmpzero, hdu=0, unit='adu')
	else:
		mzero = CCDData.read(tmpzero, hdu=0)#, unit='adu')
	mzero.meta['FILENAME'] = os.path.basename(tmpzero)
'''
#	RASA36 discreminator for high/hdr modes with bias level
if obs == 'RASA36':
	biaslevel = np.median(mzero)
	if biaslevel > 87:
		mode = 'hdr'
	else:
		mode = 'high'
	print(f'RASA36 [{mode} mode] : bias level = {biaslevel}')
'''
#------------------------------------------------------------
#	DARK (ITERATION FOR EACH EXPOSURE TIMES)
#------------------------------------------------------------
try:
	darkexptimelist = sorted(list(set(ic1.filter(imagetyp='dark').summary['exptime'])))
	darknumb = len(darkexptimelist)
except:
	darknumb = 0
darkdict = dict()
# if len(darkexptimelist) != 0:
if darknumb != 0:
	dark_process = True
	for i, exptime in enumerate(darkexptimelist):
		print('PRE PROCESS FOR DARK ({} sec)\t[{}/{}]'.format(exptime, i+1, len(darkexptimelist)))
		mdark = calib.master_dark(ic1, mzero=mzero, exptime=exptime, fig=False)
		darkdict['{}'.format(int(exptime))] = mdark

		date = fits.getheader(f'{path_data}/dark-{int(exptime)}.fits')['date-obs'][:10].replace('-', '')
		if obs == 'RASA36':
			darkim = f'{path_mframe}/{obs}/dark/{int(exptime)}-{date}-dark_{mode}.fits'
		else:
			darkim = f'{path_mframe}/{obs}/dark/{int(exptime)}-{date}-dark.fits'
		# print(zeroim)
		cpcom = 'cp {}/dark-{}.fits {}'.format(path_data, int(exptime), darkim)
		print(cpcom)
		os.system(cpcom)
		plt.close('all')
else:
	#	Borrow
	print('\nNO DARK FRAMES\n')
	objexptimelist = sorted(list(set(ic1.filter(imagetyp='object').summary['exptime'])))
	exptime = objexptimelist[-1]
	# pastdark = np.array(glob.glob('{}/{}/dark/{}*dark*.fits'.format(path_mframe, obs, int(exptime))))
	if obs == 'RASA36':
		pastdark = np.array(glob.glob(f'{path_mframe}/{obs}/dark/{int(exptime)}*dark_{mode}.fits'))
	else:
		pastdark = np.array(glob.glob(f'{path_mframe}/{obs}/dark/{int(exptime)}*dark.fits'))

	if len(pastdark) == 0:
		pastdark = np.array(glob.glob('{}/{}/dark/*dark*.fits'.format(path_mframe, obs)))
	else:
		pass
	#	CALCULATE CLOSEST ONE FROM TIME DIFFERENCE
	deltime = []
	delexptime = []
	darkexptimes = []
	for date in pastdark:
		# darkmjd = calib.isot_to_mjd((os.path.basename(date)).split('-')[0])
		darkmjd = calib.isot_to_mjd((os.path.basename(date)).split('-')[1])
		darkexptime = int( os.path.basename(date).split('-')[0] )
		# darkexptime = delexptime.append(int( os.path.basename(date).split('-')[1] ))
		darkexptimes.append(darkexptime)
		deltime.append(np.abs(ic1.summary['mjd'][0]-darkmjd))
	if 'KCT' in obs:
		indx_closet = np.where(
			(np.abs(np.array(darkexptimes)-exptime) == np.min(np.abs(np.array(darkexptimes)-exptime)))
		)
	else:
		indx_closet = np.where(
			(deltime == np.min(deltime)) &
			(darkexptimes == np.max(darkexptimes))
		)
	if len(indx_closet[0]) == 0:
		indx_closet = np.where(
			(deltime == np.min(deltime))
		)
	else:
		pass
	# tmpdark = path_data+'/'+os.path.basename(pastdark[indx_closet].item())
	# tmpdark = '{}/{}'.format(path_data, os.path.basename(pastdark[indx_closet[0]].item()))
	# tmpdark = pastdark[indx_closet[0]].item()
	tmpdark = pastdark[indx_closet][-1]
	exptime = int(fits.getheader(tmpdark)['exptime'])

	# cpcom = 'cp {} {}/dark-{}.fits'.format(tmpdark, path_data, int(exptime))
	cpcom = 'cp {} {}'.format(tmpdark, path_data, int(exptime))
	print(cpcom)
	os.system(cpcom)
	if 'KCT' in obs:
		#KCT Exception
		mdark = CCDData.read(tmpdark, hdu=0, unit='adu')
	elif obs == 'RASA36':
		if (mode == 'high') & (badmode == True):
			mdark = CCDData.read(tmpdark, hdu=0).multiply(20)
			print('[Bad mode] Multiply 20 on the high mode dark.')
		else:
			mdark = CCDData.read(tmpdark, hdu=0)
	elif obs == 'LSGT':
		mdark = CCDData.read(tmpdark, hdu=0, unit='adu')
	else:
		mdark = CCDData.read(tmpdark, hdu=0)#, unit='adu')
	mdark.meta['FILENAME'] = os.path.basename(tmpdark)
	mdark.meta['EXPTIME'] = exptime
	darkdict['{}'.format(int(exptime))] = mdark
#------------------------------------------------------------
#	FLAT (ITERATION FOR EACH FILTERS)
#------------------------------------------------------------
flatdict = dict()
try:
	flatfilterlist = list(set(ic1.filter(imagetyp='flat').summary['filter']))
	for i, filte in enumerate(flatfilterlist):
		# print(i, filte)
		print('MAKING MASTER FLAT IN {}-BAND'.format(filte))
		mflat = calib.master_flat(ic1, mzero, filte, mdark=mdark, fig=True)
		flatdict[filte] = mflat

		date = fits.getheader(f'{path_data}/dark-{int(exptime)}.fits')['date-obs'][:10].replace('-', '')
		if obs == 'RASA36':
			flatim = f'{path_mframe}/{obs}/flat/{date}-n{filte}_{mode}.fits'
		else:
			flatim = f'{path_mframe}/{obs}/flat/{date}-n{filte}.fits'
		cpcom = f'cp {path_data}/n{filte}.fits {flatim}'
		print(cpcom)
		os.system(cpcom)
		plt.close('all')
except:
	print('No flat calibration image.')
	# flatdict['None'] = None
	pass
# tdict['masterframe'] = time.time() - st
protbl['status'][protbl['process']=='master_frame'] = True
protbl['time'][protbl['process']=='master_frame'] = int(time.time() - st)
#------------------------------------------------------------
#	OBJECT CALIBARTION (ZERO, DARK, FLAT)
#------------------------------------------------------------
st_ = time.time()
comment     = '='*60+'\n' \
			+ 'OBJECT CALIBRATION\n' \
			+ '='*60+'\n'
print(comment)
objfilterlist = sorted(list(set(ic1.filter(imagetyp='object').summary['filter'])))
objexptimelist = sorted(list(set(ic1.filter(imagetyp='object').summary['exptime'])))
for i, filte in enumerate(objfilterlist):
	print('PRE PROCESS FOR {} FILTER OBJECT\t[{}/{}]'.format(filte, i+1, len(objfilterlist)))
	if filte in flatdict.keys():
		mflat = flatdict[filte]
	else:
		print('\nNO {} FLAT FRAMES\n'.format(filte))
		#	CALCULATE CLOSEST ONE FROM TIME DIFFERENCE
		deltime = []
		if obs != 'RASA36':
			pastflat = np.array(glob.glob('{}/{}/flat/*n{}*.fits'.format(path_mframe, obs, filte)))
			for date in pastflat:
				flatmjd = calib.isot_to_mjd((os.path.basename(date)).split('-')[0])
				deltime.append(np.abs(ic1.summary['mjd'][0]-flatmjd))
		elif obs == 'RASA36':
			pastflat = np.array(glob.glob(f'{path_mframe}/{obs}/flat/*_{mode}-n{filte}*.fits'))
			for date in pastflat:
				flatmjd = fits.getheader(date)['MJD']
				deltime.append(np.abs(ic1.summary['mjd'][0]-flatmjd))

		indx_closet = np.where(deltime == np.min(deltime))
		tmpflat = '{}/{}'.format(path_data, os.path.basename(pastflat[indx_closet][0].item()))
		# tmpflat = pastflat[indx_closet][0].item()
		cpcom = 'cp {} {}'.format(pastflat[indx_closet][0].item(), tmpflat)
		print(cpcom)
		os.system(cpcom)
		if ('KCT' not in obs) & (obs != 'LSGT'):
			mflat = CCDData.read(tmpflat, hdu=0)#, unit='adu')
		elif obs == 'LSGT':
			mflat = CCDData.read(tmpflat, hdu=0, unit='adu')
		else:
			#KCT Exception
			mflat = CCDData.read(tmpflat, hdu=0, unit='adu')
			mflat.meta['FILENAME'] = os.path.basename(tmpflat)
		mflat.meta['FILENAME'] = os.path.basename(tmpflat)
		flatdict[filte] = mflat
	for expt in objexptimelist:
		if str(int(expt)) in darkdict.keys():
			mdark = darkdict[str(int(expt))]
		else:
			mdark = darkdict[list(darkdict.keys())[-1]]
		calib.calibration(ic1, mzero, mflat, filte, mdark=mdark)
# tdict['objectcorrection'] = time.time() - st - tdict[list(tdict.keys())[-1]]
protbl['status'][protbl['process']=='pre_process'] = True
protbl['time'][protbl['process']=='pre_process'] = int(time.time() - st_)
#	Corrected image list
fzimlist = []
for ims in ('{}/fz*.fits'.format(path_data), '{}/fz*.fit'.format(path_data), '{}/fz*.fts'.format(path_data)):
	fzimlist.extend(sorted(glob.glob(ims)))
# print("Skip master frame & object correction")
# #	Corrected image list
# fzimlist = []
# for ims in ('{}/*.fits'.format(path_data), '{}/*.fit'.format(path_data), '{}/*.fts'.format(path_data)):
# 	fzimlist.extend(sorted(glob.glob(ims)))
# pass
#------------------------------------------------------------
#	ASTROMETRY
#------------------------------------------------------------
st_ = time.time()
print('ASTROMETRY START')
print('='*60)
astrometryfailist = []
# fzimlist = sorted(glob.glob(path_data+'/fz*.fits'))
astimlist = []
astotlist = []
astralist = []
astdelist = []

for inim in fzimlist:
	obj = (fits.getheader(inim)['object']).upper()
	if (obj in alltbl['obj']):
		indx_target = np.where(obj == alltbl['obj'])[0][0]
		ra, dec = alltbl['ra'][indx_target].item(), alltbl['dec'][indx_target].item()
		astimlist.append(inim)
		astralist.append(ra)
		astdelist.append(dec)
	elif (obs == 'RASA36') & (obj == 'GRB221009A'):
		"""
		#	Generate default.sex
		os.system(f'sex -d > {path_data}/default.sex')
		
		#	Modify DETECT_THRESH
		f = open(f'{path_data}/default.sex', 'r')
		lines = f.readlines()
		f.close()

		for ll, l in enumerate(lines):
			if 'DETECT_THRESH' in l:
				lines[ll] = 'DETECT_THRESH 10 # EDIT\n'
		nlines = ''.join(lines)
		# print(nlines)

		#	Re-write default.sex 
		with open(f'{path_data}/default.sex', 'w') as f:
			f.write(nlines)
		"""		
		# os.system(f"cp {path_config}/default.* {path_data}")
		os.system(f"cp {path_config}/default.* .")
		
		ra, dec = 288.245, 19.739
		radius = 1.898
		loscl = 2.36*0.95
		upscl = 2.36*1.05
		outname = f'{os.path.dirname(inim)}/a{os.path.basename(inim)}'
		cpulimit=60
		com = f'solve-field {inim} --scale-unit arcsecperpix --scale-low {loscl} --scale-high {upscl} --no-plots --new-fits {outname} --overwrite --use-sextractor --cpulimit {cpulimit} --ra {ra} --dec {dec} --radius {radius}'
		print(com)
		os.system(com)
	else:
		astotlist.append(inim)
#	Astrometry (IMSNG field)
if __name__ == '__main__':
	with multiprocessing.Pool(processes=ncores) as pool:
		results = pool.starmap(calib.astrometry, zip(astimlist, repeat(obsinfo['pixscale']), astralist, astdelist, repeat(obsinfo['fov']/60.), repeat(15)))
#	Astrometry (non IMSNG field)
if __name__ == '__main__':
	with multiprocessing.Pool(processes=ncores) as pool:
		results = pool.starmap(calib.astrometry, zip(astotlist, repeat(obsinfo['pixscale']), repeat(None), repeat(None), repeat(None), repeat(60)))
#	Astrometry (failed IMSNG field)
astfailist = []
for inim in astimlist:
	if (os.path.exists('{}/a{}'.format(path_data, os.path.basename(inim))) == False):
		astfailist.append(inim)
if __name__ == '__main__':
	with multiprocessing.Pool(processes=ncores) as pool:
		results = pool.starmap(calib.astrometry, zip(astfailist, repeat(obsinfo['pixscale']), repeat(None), repeat(None), repeat(None), repeat(60)))
for inim in astfailist:
	if (os.path.exists('{}/a{}'.format(path_data, os.path.basename(inim))) == False):
		astrometryfailist.append('{}/a{}'.format(path_data, os.path.basename(inim)))

os.system('rm '+path_data+'/*.axy '+path_data+'/*.corr '+path_data+'/*.xyls '+path_data+'/*.match '+path_data+'/*.rdls '+path_data+'/*.solved '+path_data+'/*.wcs ')
print('ASTROMETRY COMPLETE\n'+'='*60)
# tdict['astronomy'] = time.time() - st - tdict[list(tdict.keys())[-1]]
protbl['status'][protbl['process']=='astrometry'] = True
protbl['time'][protbl['process']=='astrometry'] = int(time.time() - st_)
#------------------------------------------------------------
#	Quick seeing measurement with SE & Cosmic ray removal
#------------------------------------------------------------
st_ = time.time()
print('Quick seeing measurement with SE & Cosmic ray removal')
print('='*60)
gain = ccdinfo['gain'].value
rdnoise = ccdinfo['rdnoise']

# afzimlist = sorted(glob.glob(path_data+'/afz*.fits'))
afzimlist = []
for ims in ('{}/a*.fits'.format(path_data), '{}/a*.fit'.format(path_data), '{}/a*.fts'.format(path_data)):
	afzimlist.extend(sorted(glob.glob(ims)))
outimlist = []
for i, inim in enumerate(afzimlist):
	outim = '{}/cr{}'.format(os.path.dirname(inim), os.path.basename(inim))
	outimlist.append(outim)
if ('KCT' not in obs) & ('RASA36' not in obs) & ('LOAO_FLI' not in obs):
	#	Seeing measurement
	if __name__ == '__main__':
		with multiprocessing.Pool(processes=ncores) as pool:
			results = pool.starmap(tool_tbd.SE_seeing, zip(afzimlist, repeat(obs), repeat(path_obs), repeat(path_config), repeat(3*u.arcsecond), repeat(0.95), repeat(True)))
	#	Remove cosmic-ray
	if __name__ == '__main__':
		with multiprocessing.Pool(processes=ncores) as pool:
			results = pool.starmap(calib.cr_removal, zip(afzimlist, outimlist, repeat(gain), repeat(rdnoise)))
else:
	print('Skip Seeing measurement & CR remove processes for {}'.format(obs))
	for inim, outim in zip(afzimlist, outimlist):
		cpcom = 'cp {} {}'.format(inim, outim)
		print(cpcom)
		os.system(cpcom)
protbl['status'][protbl['process']=='cr_removal'] = True
protbl['time'][protbl['process']=='cr_removal'] = int(time.time() - st_)
#------------------------------------------------------------
#	FILE NAME CHANGE
#------------------------------------------------------------
# for inim in sorted(glob.glob(path_data+'/afz*.fits')):
# for inim in sorted(glob.glob(path_data+'/crafz*.fits')):
# 	obj = (fits.getheader(inim)['object']).upper()
# 	#	consider IMSNG galaxy
# 	if len(astrometryfailist) != 0:
# 		print('Astrometry fail list:', astrometryfailist)
# 		if (obj in alltbl['obj']) & (inim in astrometryfailist):
# 			robj, sep = tool_tbd.imsng_name_correction(inim, alltbl, radius=obsinfo['fov']*u.arcmin)
# 		else:
# 			pass
# 	else:
# 		pass
# 	calib.fnamechange(inim, obs)

fov = obsinfo['fov']*u.arcmin 
crafzimlist = []
for ims in ('{}/cra*.fits'.format(path_data), '{}/cra*.fit'.format(path_data), '{}/cra*.fts'.format(path_data)):
	crafzimlist.extend(sorted(glob.glob(ims)))
# for inim in sorted(glob.glob('{}/crafz*.fits'.format(path_data))): 
for inim in crafzimlist:
	obj = fits.getheader(inim)['object'] 
	#	Modify incorrect object header
	if (inim.replace('crafz', 'afz') in astrometryfailist) & (obj in alltbl['obj']): 
		robj, sep = tool_tbd.imsng_name_correction(inim, alltbl, radius=fov) 
	else:
		pass
	calib.fnamechange(inim, obs)

caliblist = sorted(glob.glob(path_data+'/Calib*.fits'))
ic_cal = ImageFileCollection(path_data, glob_include='Calib*0.fits', keywords='*')
os.system('chmod 777 {}'.format(path_data))
os.system('chmod 777 {}/*'.format(path_data))
#	Calib-*.fits TO SAVE PATH
f = open(path_data+'/object.txt', 'a')
f.write('obs obj dateobs filter exptime\n')
for inim in caliblist:
	img = os.path.basename(inim)
	part = img.split('-')
	line = '{} {} {} {} {}\n'.format(part[1], part[2], part[3]+'T'+part[4], part[5], part[6])
	print(line)
	f.write(line)
f.close()

#	DATA FOLDER TO SAVE PATH
# os.system('rm {}/afz*.fits {}/fz*.fits'.format(path_data, path_data))
os.system(f'rm {path_data}/*fz*.f*')
os.system(f'rm -rf {path_save}/{os.path.basename(path_data)}')
plt.close('all')
#------------------------------------------------------------
#	Defringe for LOAO I-band images
#------------------------------------------------------------
st_ = time.time()
if (obs == 'LOAO') & ('I' in ic_cal.filter(imagetyp='object').summary['filter']):
	dfim = '/home/paek/qsopy/fringe/LOAO/fringe_i_ori.fits'
	dfdat = '/home/paek/qsopy/fringe/LOAO/fringe_i.dat'
	dfimlist = []
	for inim in ic_cal.filter(imagetyp='object', filter='I').summary['file']:
		# dfimlist.append(calib.defringe(str(inim), dfim, dfdat))
		dfedim = calib.defringe(str(inim), dfim, dfdat)
		mvcom = 'mv {} {}'.format(dfedim, inim)
		print(mvcom)
		os.system(mvcom)
	# tdict['defringe'] = time.time() - st - tdict[list(tdict.keys())[-1]]
else:
	print('No images to defringe')
	pass
protbl['status'][protbl['process']=='defringe'] = True
protbl['time'][protbl['process']=='defringe'] = int(time.time() - st_)
#------------------------------------------------------------
print('='*60)
print('Calibration IS DONE.\t('+str(int(time.time() - starttime))+' sec)')
print('='*60)
#------------------------------------------------------------
#	Photometry for single images
#------------------------------------------------------------
st_ = time.time()
print('#\tPhotometry')
# path_data = '{}/{}'.format(path_factory, os.path.basename(path))
path_infile = '{}/{}'.format(path_data, os.path.basename(path_default_gphot))
path_new_gphot = '{}/gphot.config'.format(os.path.dirname(path_infile))
#	Read default photometry configuration
os.system('cp {} {}'.format(path_default_gphot, path_new_gphot))
f = open(path_default_gphot, 'r')
lines = f.read().splitlines()
f.close()
#	Write photometry configuration
g = open(path_new_gphot, 'w')
for line in lines:
	if 'imkey' in line:
		line = '{}\t{}/C*0.fits'.format('imkey', path_data)
	else:
		pass
	g.write(line+'\n')
g.close()
#	Execute
if obs == 'DOAO':
	path_phot = path_phot_sg
else:
	path_phot = path_phot_mp
com = 'python {} {} {}'.format(path_phot, path_data, ncores)
print(com)
os.system(com)
# tdict['photometry'] = time.time() - st - tdict[list(tdict.keys())[-1]]
protbl['status'][protbl['process']=='photometry'] = True
protbl['time'][protbl['process']=='photometry'] = int(time.time() - st_)

#------------------------------------------------------------
#	Image registering & combine
#------------------------------------------------------------
st_ = time.time()
print('IMAGE REGISTERING & COMBINE')
combined_images = []
step = (1/24/60)*60	# 1 hour
ic_cal_phot = ImageFileCollection(path_data, glob_include='Calib*0.fits', keywords='*')
calist = sorted(glob.glob('{}/Calib*.fits'.format(path_data)))
objlist = sorted(list(set(ic_cal_phot.summary['object'])))
filterlist = sorted(list(set(ic_cal_phot.summary['filter'])))
# obj = 'NGC3147'
# filte = 'R'
for obj in objlist:
	for filte in filterlist:
		imlist_tmp = sorted(glob.glob('{}/Calib*-{}-*-{}-*.fits'.format(path_data, obj, filte)))
		if len(imlist_tmp) == 0:
			pass
		elif len(imlist_tmp) == 1:
			inim = imlist_tmp[0]
			comim = inim.replace('.fits', '.com.fits')
			cpcom = f'cp {inim} {comim}'
			print(cpcom)
			os.system(cpcom)
		else:
			print(obj, filte)
			# ic_part = sorted(glob.glob('{}/Calib*{}*{}*.fits'.format(path_data, obj, filte)))

			jds = np.array([fits.getheader(inim)['jd'] for inim in imlist_tmp])
			delts = jds - np.min(jds)

			grouplist = []
			grouplists = []

			i = 0
			for i in range(len(delts)):
				#	Initial setting
				if i == 0:
					t0 = delts[i]
				#	Add last group to grouplists
				elif i == len(delts)-1:
					grouplists.append(grouplist)
				t1 = delts[i]
				# print(t0, t1)
				dif = np.abs(t0-t1)

				if dif < step:
					grouplist.append(imlist_tmp[i])
				#	Generate new group
				else:
					grouplists.append(grouplist)
					grouplist = [imlist_tmp[i]]	
				t0 = t1
				
			for group in grouplists:
				print('-'*60)
				if len(group) > 1:
					images_to_align = group
					ref_image = images_to_align[0]
					# print(images_to_align, ref_image)
					for inim in images_to_align: print(inim)
					try:
						outim = tool_tbd.imcombine_routine(images_to_align, ref_image)
						combined_images.append(outim)
					except:
						print('Fail to image align & combine routine.')
						print(images_to_align)
						pass
				else:
					print('There is only one image.')
					combined_images.append(group[0])


rmcom = 'rm {}/*Calib*gregister.fits'.format(path_data)
print(rmcom)
os.system(rmcom)
# tdict['imagecombine'] = time.time() - st - tdict[list(tdict.keys())[-1]]
protbl['status'][protbl['process']=='image_stack'] = True
protbl['time'][protbl['process']=='image_stack'] = int(time.time() - st_)
#------------------------------------------------------------
#	Photometry for combined images
#------------------------------------------------------------
st_ = time.time()
#	Write photometry configuration
h = open(path_new_gphot, 'w')
for line in lines:
	if 'imkey' in line:
		line = '{}\t{}/C*com.fits'.format('imkey', path_data)
	else:
		pass
	h.write(line+'\n')
h.close()
#	Execute
path_phot = path_phot_mp
com = 'python {} {} {}'.format(path_phot, path_data, ncores)
print(com)
os.system(com)
# tdict['photometry_com'] = time.time() - st - tdict[list(tdict.keys())[-1]]
protbl['status'][protbl['process']=='photometry_com'] = True
protbl['time'][protbl['process']=='photometry_com'] = int(time.time() - st_)

ic_com_phot = ImageFileCollection(path_data, glob_include='Calib*com.fits', keywords='*')	
#	Summary
print('Draw observation summary plots')
# for filte in list(set(ic_cal_phot.summary['filter'])):
for filte in filterlist:
	try:
		tool_tbd.obs_summary(filte, ic_cal_phot, ic_com_phot, path_save=path_data)
	except:
		print('Fail to make summary plots.')
		pass
	plt.close('all')
#------------------------------------------------------------
#	Image subtraction
#------------------------------------------------------------
print('IMAGE SUBTRACTION')
subtracted_images = []
ds9comlist = []
for inim in combined_images:
	hdr = fits.getheader(inim)
	# obs = os.path.basename(inim).split('-')[1]
	# obs = 'LOAO'
	obj = hdr['object']
	filte = hdr['filter']
	path_refim = '/data3/paek/factory/ref_frames/{}'.format(obs)
	refimlist = glob.glob('{}/Ref*{}*{}*.fits'.format(path_refim, obj, filte))
	if len(refimlist) > 0:
		refim = refimlist[0]

		if obs not in ['LSGT', 'DOAO', 'RASA36', 'SAO_C361K',]:
			subim, ds9com = tool_tbd.subtraction_routine(inim, refim)
		else:
			subim, ds9com = tool_tbd.subtraction_routine2(inim, refim)
			if os.path.getsize(subim) != 0:
				rmcom = f"rm {subim}"
				print(rmcom)
				os.system(rmcom)
				subim, ds9com = tool_tbd.subtraction_routine(inim, refim)
			else:
				pass

		subtracted_images.append(subim)
		ds9comlist.append(ds9com)
	else:
		print('There is no reference image for {}'.format(os.path.basename(inim)))
		pass
rmcom = 'rm {}/*Ref*gregister.fits'.format(path_data)
print(rmcom)
os.system(rmcom)
# tdict['subtraction'] = time.time() - st - tdict[list(tdict.keys())[-1]]
protbl['status'][protbl['process']=='subtraction'] = True
protbl['time'][protbl['process']=='subtraction'] = int(time.time() - st_)
#------------------------------------------------------------
#	Remove bad subtraction image
#------------------------------------------------------------
import glob
imlist = sorted(glob.glob(f'{os.path.dirname(inim)}/hd*m.fits')) 
if len(imlist) > 0:
	from astropy.table import Table
	hdimtbl = Table()
	hdimtbl['image'] = imlist
	hdimtbl['n'] = 0
	hdimtbl['zero'] = 0
	hdimtbl['1e-30'] = 0
	hdimtbl['nan'] = 0
	hdimtbl['flag'] = 0

	import numpy as np 
	for ii, inim in enumerate(imlist): 
		data = fits.getdata(inim).flatten() 
		indx_zero = np.where(data == 0) 
		indx_nan = np.isnan(data)
		indx_fail = np.where(data == 1e-30) 
		
		n = len(data)
		nzero = len(data[indx_zero])
		nfail = len(data[indx_fail])
		nnan = len(data[indx_nan])

		#	
		hdimtbl['n'][ii] = n
		hdimtbl['zero'][ii] = nzero
		hdimtbl['1e-30'][ii] = nfail
		hdimtbl['nan'][ii] = nnan
		
		if (nzero/n > 0.5) | (nnan/n > 0.5) | (nnan/n > 0.3):
			flag = 0
		else:
			flag = 1

		hdimtbl['flag'][ii] = flag

		print(ii, inim)
		print(f'zero {1e2*len(data[indx_zero])/n:1.3f}%') 
		print(f'nan {1e2*len(data[indx_nan])/n:1.3f}%') 
		print(f'1e-30 {1e2*len(data[indx_fail])/n:1.3f}%') 
		print()
	#	Make a bad subt. image list
	import os
	hdimtbl.write(f'{os.path.dirname(inim)}/bad.subtraction.txt', format='ascii.tab', overwrite=True)

	for inim, flag in zip(hdimtbl['image'], hdimtbl['flag']):
		if flag == 0:
			#	Bad subt. image
			os.system(f'mv {inim} {inim}.bad')
else:
	print('There is no subtracted images.')

#------------------------------------------------------------
#	Photometry for subtracted images
#------------------------------------------------------------
st_ = time.time()
#	Write photometry configuration
s = open(path_new_gphot, 'w')
for line in lines:
	if 'imkey' in line:
		# line = '{}\t{}/hd*com.fits'.format('imkey', path_data)
		line = '{}\t{}/hd*.fits'.format('imkey', path_data)
	else:
		pass
	if 'photfraction' in line:
		line = '{}\t{}'.format('photfraction', 1.0)
	else:
		pass
	if 'DETECT_MINAREA' in line:
		line = '{}\t{}'.format('DETECT_MINAREA', 10)
	else:
		pass
	if 'DETECT_THRESH' in line:
		line = '{}\t{}'.format('DETECT_THRESH', 1.25)
	else:
		pass
	s.write(line+'\n')
s.close()
#	Execute
hdimlist = sorted(glob.glob('{}/hd*.fits'.format(path_data)))
if len(hdimlist) > 0:
	com = 'python {} {}'.format(path_phot_sub, path_data)
	print(com)
	os.system(com)
	# tdict['photometry_sub'] = time.time() - st - tdict[list(tdict.keys())[-1]]
else:
	print('No subtracted image.')
	pass
protbl['status'][protbl['process']=='photometry_sub'] = True
protbl['time'][protbl['process']=='photometry_sub'] = int(time.time() - st_)
#------------------------------------------------------------
#	Transient Search
#------------------------------------------------------------
st_ = time.time()
fovval = fov.value
#	Input table for transient search
tstbl = Table()
# hdimlist = sorted(glob.glob(f'{path_data}/hd*com.fits'))
hdimlist = sorted(glob.glob(f'{path_data}/hd*.fits'))
if len(hdimlist) != 0:
	tstbl['hdim'] = hdimlist
	tskeys = ['hdcat', 'hcim', 'inim', 'scicat', 'refim']
	for key in tskeys:
		tstbl[key] = ' '*300
	tstbl['fovval'] = fovval

	for i, hdim in enumerate(hdimlist):
		hdcat = hdim.replace('.fits','.phot_sub.cat')
		hcim = hdim.replace('hdCalib', 'hcCalib')
		inim = hdim.replace('hdCalib', 'Calib')
		scicat = inim.replace('.fits', '.phot.cat')

		hdr = fits.getheader(hdim)
		obj = hdr['object']
		filte = hdr['filter']
		path_refim = f'/data3/paek/factory/ref_frames/{obs}'
		refimlist = glob.glob(f'{path_refim}/Ref*{obj}*{filte}*.fits')
		refim = refimlist[0]


		for key, im in zip(tskeys, [hdcat, hcim, inim, scicat, refim]):
			tstbl[key][i] = im

	out_tstbl = f'{path_data}/transient_search.txt'
	tstbl.write(out_tstbl, format='ascii.tab', overwrite=True)

	com = f'python {path_find} {out_tstbl} {ncores}'
	print(com)
	subprocess.call(com, shell=True)		


protbl['status'][protbl['process']=='transient_search'] = True
protbl['time'][protbl['process']=='transient_search'] = int(time.time() - st_)
#------------------------------------------------------------
#	Summary file
#------------------------------------------------------------
protbl['status'][protbl['process']=='total'] = True
protbl['time'][protbl['process']=='total'] = int(time.time() - st)	
protbl.write('{}/obs.summary.log'.format(path_data), format='ascii.tab', overwrite=True)
print(protbl)
#	Write data summary
f = open(path_data+'/obs.summary.log', 'a')
end_localtime = time.strftime('%Y-%m-%d %H:%M:%S (%Z)', time.localtime())
f.write('Pipelne start\t: {}\n'.format(start_localtime))
f.write('Pipelne end\t: {}\n'.format(end_localtime))
try:
	f.write('='*60+'\n')
	f.write('PATH :{}\n'.format(path))
	f.write('OBJECT NUMBER # :{}\n'.format(len(ic_cal.summary)))
	objkind = sorted(set(ic_cal.summary['object']))
	f.write('OBJECTS # : {}\n'.format(objkind))
	for obj in objkind:
		f.write('-'*60+'\n')
		for filte in list(set(ic_cal.summary['filter'])):
			indx_tmp = ic_cal.files_filtered(filter=filte, object=obj)
			if len(indx_tmp) > 0:
				f.write('{}\t{}\n'.format(obj, filte))
except:
	pass
f.close()

#------------------------------------------------------------
#	Slack message
#------------------------------------------------------------	
total_time = round(protbl['time'][protbl['process']=='total'].item()/60., 1)

channel = '#pipeline'
text = f'[Pipeline/{obs}] End Processing {os.path.basename(path)} Data ({nobj} objects) with {ncores} cores taking {total_time} mins'

param_slack = dict(
	token = OAuth_Token,
	channel = channel,
	text = text,
)

tool_tbd.slack_bot(**param_slack)
#------------------------------------------------------------
# tdict['total'] = time.time() - st
# tool_tbd.dict2table(tdict, '{}/{}'.format(path_data, 'times.log'))
#------------------------------------------------------------
#	Finish the process
#------------------------------------------------------------
# rmcom = 'rm {}/hc*gregister.fits {}/inv*.fits'.format(path_data, path_data)
rmcom = 'rm {}/inv*.*'.format(path_data, path_data)
print(rmcom)
os.system(rmcom)
tails = ['.transients.', '.new.', '.ref.', '.sub.', '']
for obj in objlist:
	for filte in filterlist:
		for tail in tails:
			# tail = 'transients'
			# obj = 'NGC3147'
			# filte = 'B'

			pathto = f'{path_gal}/{obj}/{obs}/{filte}'
			files = f'{path_data}/*Calib*-{obj}-*-{filte}-*{tail}*'
			nfiles = len(glob.glob(files))
			# print(files, nfiles)
			# if nfiles >0:
			# 	print(obj, filte, pathto, files, glob.glob(files)[-1])
			if nfiles !=0:
				#	Save path
				if tail == '':
					pathto = f'{path_gal}/{obj}/{obs}/{filte}'
				else:
					pathto = f'{path_gal}/{obj}/{obs}/{filte}/transients'
				#	Make path
				if (not os.path.exists(pathto)):
					os.makedirs(pathto)
				mvcom = f'mv {files} {pathto}'
				print(mvcom)
				os.system(mvcom)
#	Image transfer
mvcom = f'mv {path_data} {path_save}'
os.system(mvcom)
#	WRITE LOG
f = open(path_log, 'a')
# f.write(path_raw+'/'+os.path.basename(path_data)+'\n')
# f.write('{}/{}\n'.format(path_raw, os.path.basename(path_data)))
f.write(f'{path_raw}/{os.path.basename(path_data)}\n')
f.close()
#============================================================
#	Time
#============================================================
delt = time.time() - starttime
dimen = 'seconds'
if delt > 60.:
	delt = delt/60.
	dimen = 'mins'
if delt > 60.:
	delt = delt/60.
	dimen = 'hours'
print('ALL PROCESS IS DONE.\t({} {})'.format(round(delt, 3), dimen))
print(obs, newlist)
