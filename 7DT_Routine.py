# %% [markdown]
# # IMSNG_Rougine.ipynb
# - Automatically Process the image data from GECKO facilities and search for transients in the subtracted images with `gpPy`
# - Author: Gregory S.H. Paek (23.04.24)

# %% [markdown]
# ## Library

# %%
from __future__ import print_function, division, absolute_import
import os, sys, glob, subprocess
import numpy as np
import astropy.io.ascii as ascii
import matplotlib.pyplot as plt
plt.ioff()
from astropy.nddata import CCDData
from preprocess import calib
from util import tool
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

# %%
#	plot setting
import matplotlib.pyplot as plt
import matplotlib as mpl
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "last_expr"

mpl.rcParams["axes.titlesize"] = 14
mpl.rcParams["axes.labelsize"] = 20
plt.rcParams['savefig.dpi'] = 500
plt.rc('font', family='serif')

# %% [markdown]
# ## Ready

# %%
try:
	obs = (sys.argv[1]).upper()
except:
	obs = input(f"7DT## (e.g. 7DT01):").upper()

# obs = '7DT01'
print(f'# Observatory : {obs.upper()}')

try:
	ncores = int(sys.argv[2])
except:
	ncores = 2
print(f"- Number of Cores: {ncores}")

# %% [markdown]
# ### Path

# %%
path_base = '/data4/gecko/factory'
path_ref = f'{path_base}/ref_frames/{obs.upper()}'
path_factory = f'{path_base}/{obs.lower()}'
path_save = f'/data6/bkgdata/{obs.upper()}'
path_log = f'/home/paek/log/{obs.lower()}.log'
path_keys = '/home/paek/table'
#------------------------------------------------------------
path_gal = '/data6/IMSNG/IMSNGgalaxies'
path_refcat = '/data4/gecko/factory/ref_frames/LOAO'
#------------------------------------------------------------
# path_config = '/home/paek/config'
path_config = './config'
path_default_gphot = f'{path_config}/gphot.{obs.lower()}.config'
path_mframe = f'/data3/paek/factory/master_frames'
path_calib = f'{path_base}/calib'
#------------------------------------------------------------
#	Codes
#------------------------------------------------------------
path_phot_sg = './phot/gregoryphot_2021.py'
path_phot_mp = './phot/gregoryphot_mp_2021.py'
path_phot_sub = './phot/gregoryphot_sub_2021.py'
path_find = './phot/gregoryfind_bulk_mp_2021.py'
#------------------------------------------------------------
path_obsdata = '/data6/obsdata'
path_raw = f'{path_obsdata}/{obs.upper()}'
rawlist = sorted(glob.glob(path_raw+'/2*'))
#------------------------------------------------------------
path_obs = '/home/paek/table/obs.dat'
path_changehdr = '/home/paek/table/changehdr.dat'
path_alltarget = '/home/paek/table/alltarget.dat'
ccdinfo = tool.getccdinfo(obs, path_obs)

# %% [markdown]
# - Folder setting

# %%
if not os.path.exists(path_base):
	os.makedirs(path_base)
if not os.path.exists(path_ref):
	os.makedirs(path_ref)
if not os.path.exists(path_factory):
	os.makedirs(path_factory)
if not os.path.exists(path_save):
	os.makedirs(path_save)

# %% [markdown]
# ### Table

# %%
logtbl = ascii.read(path_log)
datalist = np.copy(logtbl['date'])
obstbl = ascii.read(path_obs)
hdrtbl = ascii.read(path_changehdr)
alltbl = ascii.read(path_alltarget)
keytbl = ascii.read(f'{path_keys}/keys.dat')

# %% [markdown]
# ## Process Summary Status

# %%
protbl = Table()
protbl['process'] = ['master_frame', 'pre_process', 'astrometry', 'cr_removal', 'defringe', 'photometry', 'image_stack', 'photometry_com', 'subtraction', 'photometry_sub', 'transient_search', 'total']
protbl['status'] = False
protbl['time'] = 0.0 * u.second

# %% [markdown]
# ## Main Body

# %%
newlist = [i for i in rawlist if (i not in datalist) & (i+'/' not in datalist)]
if len(newlist) == 0:
	print('No new data')
	sys.exit()
else:
	print(newlist)

# path = newlist[0]
path = newlist[-1]
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

# %% [markdown]
# ### Header Correction

# %%
ic0 = ImageFileCollection(path_data, keywords='*')
ic0.summary.write(f'{path_data}/hdr.raw.dat', format='ascii.tab', overwrite=True) 
obsinfo = calib.getobsinfo(obs, obstbl)
calib.correcthdr_routine(path_data, hdrtbl, obs)
print("Correction Done")
# objfilterlist, objexptimelist, flatfilterlist, darkexptimelist, obstime = calib.correcthdr_routine(path_data, hdrtbl)
ic1 = ImageFileCollection(path_data, keywords='*')
ic1.summary.write('{}/hdr.cor.dat'.format(path_data), format='ascii.tab', overwrite=True)

try:
	nobj = len(ic1.filter(imagetyp='OBJECT').summary)
except:
	try:
		nobj = len(ic1.filter(imagetyp='object').summary)	
	except:
		nobj = 0

# %% [markdown]
# ### Marking the `GECKO` data

# %%
# testobj = 'S190425z'
project = "7DT"
obsmode = "MONITORING" # Default
if 'OBJECT' in ic1.summary.keys():
	for obj in ic1.filter(imagetyp='object').summary['object']:
		if 'MS' in obj[:2]: # MS230425 (test event)
			print(obj)
			project = "GECKO"
			obsmode = "TEST"
		elif 'S2' in obj[:2]: # S230425 (super event)
			print(obj)
			project = "GECKO"
			obsmode = "FOLLOWUP" # Follow-up
		else:
			pass
else:
	pass
print(f"[{project}] {obsmode}")

# %% [markdown]
# - Slack notification

# %%
OAuth_Token = keytbl['key'][keytbl['name']=='slack'].item()

channel = '#pipeline'
text = f'[`gpPy`/{project}-{obsmode}] Start Processing {obs} {os.path.basename(path)} Data ({nobj} objects) with {ncores} cores'

param_slack = dict(
	token = OAuth_Token,
	channel = channel,
	text = text,
)

tool.slack_bot(**param_slack)

# %% [markdown]
# ### Master Frame

# %% [markdown]
# - Bias

# %%
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
	zeroim = f'{path_mframe}/{obs}/zero/{date}-zero.fits'
	if not os.path.exists(os.path.dirname(zeroim)):
		os.makedirs(os.path.dirname(zeroim))
	cpcom = f'cp {path_data}/zero.fits {zeroim}'
	print(cpcom)
	os.system(cpcom)
	plt.close('all')
else:
	#	IF THERE IS NO FLAT FRAMES, BORROW FROM CLOSEST OTHER DATE
	print('\nNO BIAS FRAMES\n')
	pastzero = np.array(glob.glob(f'{path_mframe}/{obs}/zero/*zero.fits'))
	#	CALCULATE CLOSEST ONE FROM TIME DIFFERENCE
	deltime = []
	for date in pastzero:
		zeromjd = calib.isot_to_mjd((os.path.basename(date)).split('-')[0])
		deltime.append(np.abs(ic1.summary['mjd'][0]-zeromjd))
	indx_closet = np.where(deltime == np.min(deltime))
	# tmpzero = path_data+'/'+os.path.basename(np.asscalar(pastzero[indx_closet]))
	tmpzero = f"{path_data}/{os.path.basename(pastzero[indx_closet][0])}"
	# cpcom = 'cp {} {}'.format(np.asscalar(pastzero[indx_closet]), tmpzero)
	cpcom = f'cp {pastzero[indx_closet][0]} {tmpzero}'
	print(cpcom)
	os.system(cpcom)
	mzero = CCDData.read(tmpzero, hdu=0)#, unit='adu')
	mzero.meta['FILENAME'] = os.path.basename(tmpzero)
delt_bias = time.time() - st
print(f"{delt_bias:.3f} sec")

# %% [markdown]
# - Dark

# %%
#------------------------------------------------------------
#	DARK (ITERATION FOR EACH EXPOSURE TIMES)
#------------------------------------------------------------
st = time.time()

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

	cpcom = 'cp {} {}/dark-{}.fits'.format(tmpdark, path_data, int(exptime))
	# cpcom = 'cp {} {}'.format(tmpdark, path_data, int(exptime))
	print(cpcom)
	os.system(cpcom)

	mdark = CCDData.read(tmpdark, hdu=0)#, unit='adu')
	mdark.meta['FILENAME'] = os.path.basename(tmpdark)
	mdark.meta['EXPTIME'] = exptime
	darkdict['{}'.format(int(exptime))] = mdark
delt_dark = time.time() - st
print(f"{delt_dark:.3f} sec")
# %% [markdown]
# - Flat

# %%
flatdict = dict()
try:
	flatfilterlist = list(set(ic1.filter(imagetyp='flat').summary['filter']))
	for i, filte in enumerate(flatfilterlist):
		# print(i, filte)
		print('MAKING MASTER FLAT IN {}-BAND'.format(filte))
		mflat = calib.master_flat(ic1, mzero, filte, mdark=mdark, fig=True)
		flatdict[filte] = mflat

		date = fits.getheader(f'{path_data}/dark-{int(exptime)}.fits')['date-obs'][:10].replace('-', '')
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
		pastflat = np.array(glob.glob('{}/{}/flat/*n{}*.fits'.format(path_mframe, obs, filte)))
		for date in pastflat:
			flatmjd = calib.isot_to_mjd((os.path.basename(date)).split('-')[0])
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

# %% [markdown]
# ### WCS Calculation

# %%
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

# %% [markdown]
# ### Cosmic-ray Removal

# %%
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
if ('KCT' not in obs) & ('RASA36' not in obs) & ('LOAO_FLI' not in obs) & ('LSGT_ASI1600MM' != obs) & ('DNSM' != obs) & ('7DT01' != obs):
	#	Seeing measurement
	if __name__ == '__main__':
		with multiprocessing.Pool(processes=ncores) as pool:
			results = pool.starmap(tool.SE_seeing, zip(afzimlist, repeat(obs), repeat(path_obs), repeat(path_config), repeat(3*u.arcsecond), repeat(0.95), repeat(True)))
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

# %% [markdown]
# ### Rename to IMSNG/GECKO Convention

# %%
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

# %%
for calibim in caliblist:

	center, vertices = tool.get_wcs_coordinates(calibim)

	fits.setval(calibim, "RACENT", value=round(center[0].item(), 3), comment="RA CENTER [deg]")	
	fits.setval(calibim, "DECCENT", value=round(center[1].item(), 3), comment="DEC CENTER [deg]")	
	
	for ii, (_ra, _dec) in enumerate(vertices):
		# print(ii, _ra, _dec)
		fits.setval(calibim, f"RAPOLY{ii}", value=round(_ra, 3), comment=f"RA POLYGON {ii} [deg]")	
		fits.setval(calibim, f"DEPOLY{ii}", value=round(_dec, 3), comment=f"DEC POLYGON {ii} [deg]")
		

# %% [markdown]
# ### Defringe
# - Only for LOAO z, I, Y-bands

# %%
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

# %% [markdown]
# ## Photometry

# %%
st_ = time.time()
print('#\tPhotometry')
path_infile = f'{path_data}/{os.path.basename(path_default_gphot)}'
path_new_gphot = f'{os.path.dirname(path_infile)}/gphot.config'

#	Copy default photometry configuration
cpcom = f'cp {path_default_gphot} {path_new_gphot}'
print(cpcom)
os.system(cpcom)

#	Read default photometry configuration
f = open(path_default_gphot, 'r')
lines = f.read().splitlines()
f.close()

#	Write photometry configuration
g = open(path_new_gphot, 'w')
for line in lines:
	if 'imkey' in line:
		line = f'imkey\t{path_data}/C*0.fits'
	else:
		pass
	g.write(line+'\n')
g.close()

if obs == 'DOAO':
	path_phot = path_phot_sg
else:
	path_phot = path_phot_mp

#	Execute
com = f'python {path_phot} {path_data} {ncores}'
# com = f'python {path_phot} {path_data} 1'
print(com)
os.system(com)

protbl['status'][protbl['process']=='photometry'] = True
protbl['time'][protbl['process']=='photometry'] = int(time.time() - st_)

# %% [markdown]
# ## Image registering & combine

# %%
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
					ref_image = group[0]
					images_to_align = group[1:]
					for inim in group:
						print(inim)
					# _data, _hdr = fits.getdata(ref_image, header=True)
					# ximage, yimage = _data.shape
					# racent, decent = _hdr['RACENT'], _hdr['DECCENT']
					try:
						# outim = tool.swarpcomb(
						# 	images_to_align,
						# 	gain=obsinfo['gain'].value,
						# 	pixscale=obsinfo['pixscale'],
						# 	racent=racent,
						# 	decent=decent,
						# 	ximage=ximage,
						# 	yimage=yimage,
						# 	listname='obj.list',
						# 	path_save='.',
						# 	keys_to_get=[
						# 		'OBJECT',
						# 		'FILTER',
						# 		'RACENT',
						# 		'DECCENT',
						# 		'RAPOLY0',
						# 		'DEPOLY0',
						# 		'RAPOLY1',
						# 		'DEPOLY1',
						# 		'RAPOLY2',
						# 		'DEPOLY2',
						# 		'RAPOLY3',
						# 		'DEPOLY3',
						# 		]
						# 	)
						outim = tool.imcombine_routine(images_to_align, ref_image)
						combined_images.append(outim)
					except:
						print('Fail to image align & combine routine.')
						print(images_to_align)
						pass
				else:
					print('There is only one image.')
					combined_images.append(group[0])

protbl['status'][protbl['process']=='image_stack'] = True
protbl['time'][protbl['process']=='image_stack'] = int(time.time() - st_)

# %%
# images_to_align = group
# ref_image = images_to_align[0]

# outim = tool.imcombine_routine(images_to_align, ref_image)

# %% [markdown]
# ## Photometry for combined images

# %%
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
		tool.obs_summary(filte, ic_cal_phot, ic_com_phot, path_save=path_data)
	except:
		print('Fail to make summary plots.')
		pass
	plt.close('all')

# %% [markdown]
# #	Image subtraction
# 

# %%
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

		# subim, ds9com = tool.subtraction_routine3(inim, refim)

		# if False:
		if obs not in ['LSGT', 'DOAO', 'RASA36', 'SAO_C361K',]:
			subim, ds9com = tool.subtraction_routine(inim, refim)
		else:
			subim, ds9com = tool.subtraction_routine2(inim, refim)
			if os.path.getsize(subim) != 0:
				rmcom = f"rm {subim}"
				print(rmcom)
				os.system(rmcom)
				subim, ds9com = tool.subtraction_routine(inim, refim)
			else:
				pass
		if subim != None:
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

# %% [markdown]
# ##	Photometry for subtracted images
# 

# %%
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

# %% [markdown]
# ##	Transient Search
# 

# %%
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

# %% [markdown]
# #	Summary file

# %%
#------------------------------------------------------------
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

# %% [markdown]
# ## File Transfer

# %%
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

# %% [markdown]
# ##	Slack message

# %%
total_time = round(protbl['time'][protbl['process']=='total'].item()/60., 1)

channel = '#pipeline'
text = f'[`gpPy`/{project}-{obsmode}] Processing Complete {obs} {os.path.basename(path)} Data ({nobj} objects) with {ncores} cores taking {total_time} mins'

param_slack = dict(
	token = OAuth_Token,
	channel = channel,
	text = text,
)

tool.slack_bot(**param_slack)


