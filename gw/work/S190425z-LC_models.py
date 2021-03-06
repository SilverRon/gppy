#	S190425z MODEL LC
#	Kasen+17
#	Tanaka+??
#	GW170817 LC
#------------------------------------------------------------
#	2019.08.21	CREATED BY Gregory S.H. Paek
#============================================================
from matplotlib import pyplot as plt
import ligo.skymap.plot
from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
import os, sys, glob
import numpy as np
import matplotlib.pyplot as plt 
from astropy.table import Table, Column, MaskedColumn, vstack
from astropy.io import ascii
from astropy.time import Time
#============================================================
#	FUNCTION
#============================================================
def asciiread(path_ascii):
	intbl = ascii.read(path_ascii)
	intbl.meta['path'] = path_ascii
	intbl.meta['name'] = os.path.basename(path_ascii)
	return intbl
#------------------------------------------------------------
def splinetable(gwtbl, filterlist):
	tblist = []
	for filte in filterlist:
		tmptbl = gwtbl[(gwtbl['Band']==filte) & (gwtbl['l_Mag']!='>')]
		tblist.append(tmptbl)
	return tblist
#------------------------------------------------------------
def kasen2table(kasentbl, filte, dist0, dist0er):
	kasentbl = kasentbl[kasentbl['day']>=0.0]
	appmag, appmagerr = abs2app(kasentbl[filte], np.zeros(len(kasentbl)), dist0, dist0er)
	filtes = np.array([filte for i in range(len(kasentbl))])
	newtbl = Table([np.copy(kasentbl['day']), filtes, appmag, appmagerr],
					names=('Phase', 'Band','Mag0', 'e_Mag0'))
	return newtbl
#------------------------------------------------------------
def tanaka2table(tanakatbl, filte, dist0, dist0er, z):
	tanakatbl = tanakatbl[tanakatbl['redshift']==z]
	appmag, appmagerr = abs2app(tanakatbl[filte], np.zeros(len(tanakatbl)), dist0, dist0er)
	filtes = np.array([filte for i in range(len(tanakatbl))])
	newtbl = Table([np.copy(tanakatbl['day']), filtes, appmag, appmagerr],
					names=('Phase', 'Band','Mag0', 'e_Mag0'))
	return newtbl
#------------------------------------------------------------
def abs2app(mag, magerr, dist, disterr):
	appmag = mag+5*np.log10(dist)-5
	appmagerr = np.sqrt(magerr**2 + ((5*disterr)/(dist*np.log(10)))**2)
	return appmag, appmagerr
#------------------------------------------------------------
def app2distscaling(mag, magerr, dist, disterr, dist0, dist0err):
	'''
	mag		: APPARENT MAGNITUE
	magerr	: ERROR
	dist	: DISTANCE
	disterr	: DISTANCE ERROR
	dist0	: DISTANCE TO SCALE
	dist0err: DISTANCE ERROR TO SCALE
	'''
	mag0 = mag - 5*np.log10(dist/dist0)
	mag0err = np.sqrt((magerr)**2 + ((5*disterr)/(dist*np.log(10)))
	                  ** 2 + ((5*dist0err)/(dist0*np.log(10)))**2)
	return mag0, mag0err
#------------------------------------------------------------
def drawLC(obstbl, obs, param_obs):
	plotbl = obstbl[obstbl['obs']==obs.upper()]
	param_plot = dict(	capsize=10, capthick=1,
						marker='v', ms=20, mec='k', ecolor='k',
						linestyle='', alpha=1.0)
	if 'kmtnet-' in obs:
		observat = obs.split('-')[1]
		plt.errorbar(plotbl['Phase'], plotbl['ul'], yerr=plotbl['ulstd'], **param_obs[observat.lower()], **param_plot)
	else:
		plt.errorbar(plotbl['Phase'], plotbl['ul'], yerr=plotbl['ulstd'], **param_obs[obs.lower()], **param_plot)
#------------------------------------------------------------
def drawLC2(obstbl, obs, param_obs):
	plotbl = obstbl[obstbl['obs']==obs.upper()]
	'''
	param_plot = dict(	capsize=10, capthick=1,
						marker='v', ms=20, 
						mec='k', ecolor='k',
						# mec=param_obs[obs]['mfc'], ecolor=param_obs[obs]['mfc'],
						linestyle='', alpha=0.75)
	'''
	if 'kmtnet-' in obs:
		observat = obs.split('-')[1]
		param_plot = dict(	capsize=10, capthick=1,
							marker='v', ms=20, 
							# mec='k', ecolor='k',
							mec=param_obs[observat]['mfc'], ecolor=param_obs[observat]['mfc'],
							linestyle='', alpha=0.75)
		'''
		plt.errorbar(plotbl['Phase'], plotbl['ul'],
					yerr=[plotbl['yer1'], plotbl['yer2']], xerr=[plotbl['xer1'], plotbl['xer2']],
					**param_obs[observat.lower()], **param_plot)
		'''
		plt.errorbar(plotbl['Phase'], plotbl['ul'],
					xerr=[plotbl['xer1'], plotbl['xer2']],
					**param_obs[observat.lower()], **param_plot)
	else:
		param_plot = dict(	capsize=10, capthick=1,
							marker='v', ms=20, 
							# mec='k', ecolor='k',
							mec=param_obs[obs]['mfc'], ecolor=param_obs[obs]['mfc'],
							linestyle='', alpha=0.75)
		'''
		plt.errorbar(plotbl['Phase'], plotbl['ul'],
					yerr=[plotbl['yer1'], plotbl['yer2']], xerr=[plotbl['xer1'], plotbl['xer2']],
					**param_obs[obs.lower()], **param_plot)
		'''
		plt.errorbar(plotbl['Phase'], plotbl['ul'],
					xerr=[plotbl['xer1'], plotbl['xer2']],
					**param_obs[obs.lower()], **param_plot)



#------------------------------------------------------------
#	PATH
#------------------------------------------------------------
path_save = '/data1/S190425z/1.result/figure2020'
path_base = '/data1/S190425z/1.result/table'
path_gw170817 = path_base+'/GW170817_Villar+17_LC.dat'
# path_obs = '/data1/S190425z/1.result/Update/phot_upd_sources_GECKO.dat'
# path_obs = '/data1/S190425z/1.result/table/obs_all.dat'
# path_obs = '/data1/S190425z/1.result/table/obs_all_kmtnetsplit.dat'
path_obs = '/data1/S190425z/1.result/table/obs_complete_20191009.dat'
#------------------------------------------------------------
path_kasen_blue = path_base+'/knova_d1_n10_m0.025_vk0.30_Xlan1e-4.0.h5_z0.009787_mag.dat'
path_kasen_red = path_base+'/knova_d1_n10_m0.040_vk0.15_Xlan1e-1.5.h5_z0.009787_mag.dat'
#------------------------------------------------------------
path_tanaka_apr4_1215 = path_base+'/knova_APR4-1215_Mej0.9e-2_LCabs.dat'	#	Table1
path_tanaka_apr4_1314 = path_base+'/knova_APR4-1314_Mej0.8e-2_LCabs.dat'	#	Table2
path_tanaka_h4_1215 = path_base+'/knova_H4-1215_Mej0.4e-2_LCabs.dat'		#	Table3
path_tanaka_h4_1314 = path_base+'/knova_H4-1314_Mej0.07e-2_LCabs.dat'		#	Table4
path_tanaka_sly_135 = path_base+'/knova_Sly-135_Mej2.0e-2_LCabs.dat'		#	Table5
#------------------------------------------------------------
path_piro = path_base+'/cocoon_Piro+18.dat'
#------------------------------------------------------------
#	SETTING
#------------------------------------------------------------
dist = 37.7*1e6								#	[Mpc]
dister = 8.7*1e6							#	[Mpc]
dist0 = 156*1e6								#	[Mpc]
dist0er = 41*1e6							#	[Mpc]
t0 = Time('2019-04-25T08:18:05.017', format='isot', scale='utc')
#------------------------------------------------------------
#	PARAM FOR PLOT
param_obs = dict(	loao =	dict(mfc='gold',
								label='LOAO'),
					sao =	dict(mfc='hotpink',
								label='SAO'),
					lsgt =	dict(mfc='yellowgreen',
								label='LSGT'),
					sso = dict(mfc='teal',
							label='KMTNet-SSO'),
					saao = dict(mfc='cyan',
							label='KMTNet-SAAO'),
					ctio = dict(mfc='deepskyblue',
							label='KMTNet-CTIO'),
					kmtnet= dict(mfc='k',
								label='KMTNet'),
					squean= dict(mfc='purple',
								label='SQUEAN'),
					ukirt= dict(mfc='slategrey',
								label='UKIRT')
				)
#------------------------------------------------------------
obstbl = asciiread(path_obs)
# obstbl['Phase'] = obstbl['jd']-t0.jd
# gwtbl = ascii.read(path_gw170817)			#	GW170817 LC
gwtbl = asciiread(path_gw170817)
#------------------------------------------------------------
# kbltbl = ascii.read(path_kasen_blue)		#	Kasen+17 BLUE KN
# krdtbl = ascii.read(path_kasen_red)		#	Kasen+17 RED KN
kbltbl = asciiread(path_kasen_blue)
krdtbl = asciiread(path_kasen_red); krdtbl['day'] = krdtbl['t']
#------------------------------------------------------------
#	DIFFERENT EOS (TANAKA)
#------------------------------------------------------------
# t1tbl = ascii.read(path_tanaka_apr4_1215)
# t2tbl = ascii.read(path_tanaka_apr4_1314)
t1tbl = asciiread(path_tanaka_apr4_1215)
t2tbl = asciiread(path_tanaka_apr4_1314)
t3tbl = asciiread(path_tanaka_h4_1215)
t4tbl = asciiread(path_tanaka_h4_1314)
t5tbl = asciiread(path_tanaka_sly_135)
#------------------------------------------------------------
#	COCOON MODEL (Piro) - DIFFERENT FORMAT WITH OTHERS (VERTICAL)
#------------------------------------------------------------
ctbl = asciiread(path_piro)
#------------------------------------------------------------
#	GW170817 LC TABLE -> FOR EACH FILTERS(BANDS)
tblist = splinetable(gwtbl, ['r', 'R', 'i', 'J', 'Ks'])
comtbl = vstack(tblist)
comtbl.write(path_save+'/gw170817_LC.dat', format='ascii', overwrite=True)
comtbl = ascii.read(path_save+'/gw170817_LC.dat')
comtbl['Mag0'], comtbl['e_Mag0'] = app2distscaling(comtbl['Mag'], comtbl['e_Mag'], dist, dister, dist0, dist0er) 
comtbl.meta['path'] = path_gw170817
comtbl.meta['name'] = 'GW170817'
#------------------------------------------------------------
#	Kasen MODEL (BLUE & RED KN)
bltblist = []
rdtblist = []
for filte in ['r', 'i', 'J', 'Ks']:
	bltblist.append(kasen2table(kbltbl, filte, dist0, dist0er))
	rdtblist.append(kasen2table(krdtbl, filte, dist0, dist0er))
bluetbl = vstack(bltblist)
redtbl = vstack(rdtblist)

# bluetbl.meta = kbltbl.meta
# redtbl.meta = krdtbl.meta
bluetbl.meta['name'] = 'Blue KN'
redtbl.meta['name'] = 'Red KN'
#------------------------------------------------------------
#	Tanaka MODEL 
t1tblist = []
t2tblist = []
t3tblist = []
t4tblist = []
t5tblist = []
for filte in ['r', 'i', 'J', 'K']:
	t1tblist.append(tanaka2table(t1tbl, filte, dist0, dist0er, z=0.01))
	t2tblist.append(tanaka2table(t2tbl, filte, dist0, dist0er, z=0.01))
	t3tblist.append(tanaka2table(t3tbl, filte, dist0, dist0er, z=0.01))
	t4tblist.append(tanaka2table(t4tbl, filte, dist0, dist0er, z=0.01))
	t5tblist.append(tanaka2table(t5tbl, filte, dist0, dist0er, z=0.01))
tt1tbl = vstack(t1tblist)
tt2tbl = vstack(t2tblist)
tt3tbl = vstack(t3tblist)
tt4tbl = vstack(t4tblist)
tt5tbl = vstack(t5tblist)

tt1tbl.meta = t1tbl.meta
tt2tbl.meta = t2tbl.meta
tt3tbl.meta = t3tbl.meta
tt4tbl.meta = t4tbl.meta
tt5tbl.meta = t5tbl.meta

#------------------------------------------------------------
#	PLOT
#------------------------------------------------------------
# band = 'i'
band = 'r'
plt.close('all')
fig = plt.figure(figsize=(10, 8))
ax = plt.subplot(111)
# for intbl in [bluetbl, redtbl, tt1tbl, tt2tbl, tt3tbl, tt4tbl, tt5tbl, comtbl]:
# for intbl in [bluetbl, redtbl, comtbl, tt1tbl, tt2tbl, tt3tbl, tt4tbl, tt5tbl]:
#------------------------------------------------------------
#	KN MODEL PLOT
#------------------------------------------------------------
# for intbl in [bluetbl, redtbl]:
for intbl in [bluetbl]:
	indx = np.where((intbl['Band']==band)&(intbl['Phase']<10.0))
	if len(indx[0]) == 0:
		indx = np.where((intbl['Band']==band+'s')&(intbl['Phase']<10.0))
	pltbl = intbl[indx]
	if intbl.meta['name'] == 'GW170817':
		ax.errorbar(intbl[indx]['Phase'], intbl[indx]['Mag0'], yerr=intbl[indx]['e_Mag0'],
					marker='o', color='grey', alpha=0.5,
					capsize=5, capthick=1,
					linestyle='None', label=intbl.meta['name'])
	if 'KN' in intbl.meta['name']:
		if 'Blue' in intbl.meta['name']:
			c, fc = 'blue', 'dodgerblue'
		elif 'Red' in intbl.meta['name']:
			c, fc = 'red', 'tomato'

		
		ax.plot(pltbl['Phase'], pltbl['Mag0']+pltbl['e_Mag0'],
				# y2=pltbl['Mag0']-pltbl['e_Mag0'],
				color=fc,
				alpha=0.5,
				)
		ax.plot(pltbl['Phase'], pltbl['Mag0']-pltbl['e_Mag0'],
				color=fc,
				alpha=0.5,
				)
		ax.fill_between(pltbl['Phase'],y1=pltbl['Mag0']+pltbl['e_Mag0'],
										y2=pltbl['Mag0']-pltbl['e_Mag0'],
										facecolor=fc,
										alpha=0.3,
										interpolate=True,
										label=intbl.meta['name'])
		
		# plt.plot(intbl[indx]['Phase'], intbl[indx]['Mag0']+intbl[indx]['e_Mag0'], color=fc, alpha=0.5)
		# plt.plot(intbl[indx]['Phase'], intbl[indx]['Mag0']-intbl[indx]['e_Mag0'], color=fc, alpha=0.5)
		# ax.plot(intbl[indx]['Phase'], intbl[indx]['Mag0'], color=c, alpha=0.5)
	elif 'knova' in intbl.meta['name']:
		ax.fill_between(pltbl['Phase'],y1=pltbl['Mag0']+pltbl['e_Mag0'],
										y2=pltbl['Mag0']-pltbl['e_Mag0'],
										# facecolor=fc,
										alpha=0.3,
										interpolate=True,
										label=(intbl.meta['name']).split('_')[1])
		# ax.plot(intbl[indx]['Phase'], intbl[indx]['Mag0'], alpha=0.5)
		# plt.plot(intbl[indx]['Phase'], intbl[indx]['Mag0'], color=c, alpha=0.5)
#------------------------------------------------------------
#	COCOON MODEL PLOT
#------------------------------------------------------------
ctbl0 = ctbl[(ctbl['filter']==band)&(ctbl['delmjd']<0.5)]
if band == 'K':
	ctbl0 = ctbl[(ctbl['filter']=='J')&(ctbl['delmjd']<3.0)]
cphase = ctbl0['delmjd']
cmag, cmager = abs2app(ctbl0['absmag'], 0, dist0, dist0er)
ax.plot(cphase, cmag+cmager,
		color='crimson',
		alpha=0.5,
		)
ax.plot(cphase, cmag-cmager,
		color='crimson',
		alpha=0.5,
		)
ax.fill_between(cphase, y1=cmag+cmager,
						y2=cmag-cmager,
						facecolor='crimson',
						alpha=0.3,
						interpolate=True,
						label='Cocoon')
# ax.plot(cphase, cmag, color='blue', alpha=0.5)

#------------------------------------------------------------
#	EOS MODEL PLOT
#------------------------------------------------------------
"""
for intbl in [tt1tbl, tt2tbl, tt3tbl, tt4tbl, tt5tbl]:
	indx = np.where((intbl['Band']==band)&(intbl['Phase']<10.0))
	if len(indx[0]) == 0:
		indx = np.where((intbl['Band']==band+'s')&(intbl['Phase']<10.0))
	pltbl = intbl[indx]
	if intbl.meta['name'] == 'GW170817':
		ax.errorbar(intbl[indx]['Phase'], intbl[indx]['Mag0'], yerr=intbl[indx]['e_Mag0'],
					marker='o', color='grey', alpha=0.5,
					capsize=5, capthick=1,
					linestyle='None', label=intbl.meta['name'])
	if 'KN' in intbl.meta['name']:
		if 'Blue' in intbl.meta['name']:
			c, fc = 'blue', 'dodgerblue'
		elif 'Red' in intbl.meta['name']:
			c, fc = 'red', 'tomato'		
		ax.fill_between(pltbl['Phase'],y1=pltbl['Mag0']+pltbl['e_Mag0'],
										y2=pltbl['Mag0']-pltbl['e_Mag0'],
										facecolor=fc,
										alpha=0.3,
										interpolate=True,
										label=intbl.meta['name'])
	elif 'knova' in intbl.meta['name']:
		'''
		ax.fill_between(pltbl['Phase'],y1=pltbl['Mag0']+pltbl['e_Mag0'],
										y2=pltbl['Mag0']-pltbl['e_Mag0'],
										# facecolor=fc,
										alpha=0.5,
										interpolate=True,
										label=(intbl.meta['name']).split('_')[1])
		'''
		# ax.plot(pltbl['Phase'], pltbl['Mag0'], alpha=0.5, label='_nolabel_')
		ax.plot(pltbl['Phase'], pltbl['Mag0'], alpha=0.75, label=(intbl.meta['name']).split('_')[1])
		# ax.errorbar(pltbl['Phase']+10, pltbl['Mag0'], pltbl['e_Mag0'], alpha=0.5, label=(intbl.meta['name']).split('_')[1])
"""
#------------------------------------------------------------
#	OBSERVATION
#------------------------------------------------------------
if band == 'r':
	# obslist = ['kmtnet', 'loao', 'lsgt', 'sao', 'kmtnet-sso',  'kmtnet-saao',  'kmtnet-ctio']
	obslist = ['loao', 'lsgt', 'sao', 'kmtnet-sso',  'kmtnet-saao',  'kmtnet-ctio']

elif band == 'i':
	obslist = ['squean']
elif band == 'K':
	obslist = ['ukirt']
for obs in obslist:
	drawLC2(obstbl, obs, param_obs)

#	GW170817 DISCOVERY TIME
ax.axvline(x=0.48, color='grey', alpha=0.5, linestyle='--')
#------------------------------------------------------------
'''
chartBox = ax.get_position()
ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.8, chartBox.height])
ax.legend(loc='upper center', bbox_to_anchor=(1.15, 1.00), shadow=False, ncol=1, fontsize=15, frameon=False)
'''
#------------------------------------------------------------
plt.gca().invert_yaxis()
# plt.title('{}-band'.format(band), fontsize=20)
# plt.xlim([0, 7])
if band == 'r':
	plt.xlim([-0.1, 4])
	plt.ylim([24.0, 19])
else:
	plt.xlim([-0.5, 4])
	plt.ylim([26.5, 18.5])
plt.xlabel('Phase [days]', fontsize=20)
plt.ylabel('Apparent magnitude', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(np.arange(18, 30, 1), fontsize=20)
# plt.legend(fontsize=15, loc='upper right')
plt.legend(fontsize=20, loc='lower right', framealpha=1.0)
plt.minorticks_on()
#------------------------------------------------------------
ax0 = ax.twinx()
ax.tick_params(which='both', direction='in', labelsize=15)
ax.minorticks_on()
ax0.tick_params(which='both', direction='in', labelsize=15)
ax0.set_ylim(29-5*np.log10(dist0)+5, 18-5*np.log10(dist0)+5)
ax0.minorticks_on()
ax0.set_ylabel('Absolute magnitude', fontsize=20, rotation=270, labelpad=20)
#------------------------------------------------------------
plt.tight_layout()
plt.savefig('{}/S190425z_obs+model_in{}band.png'.format(path_save, band), dpi=500, overwrite=True)