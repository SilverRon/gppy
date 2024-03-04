#   SELECT GW HOST GALAXY CANDIDATES
#	2019.02.10	MADE BY Gregory S.H. Paek
#	2019.08.29	UPDATED BY Gregory S.H. Paek
#============================================================#
import os, glob, sys
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
from astropy.table import Table, vstack, hstack, Column
from astropy.io import ascii
from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy import special
import time
from astropy.table import Table, Column, MaskedColumn, vstack
import astropy.units as u
from astropy.coordinates import SkyCoord
from imsng import tool
from ligo.skymap.postprocess import find_greedy_credible_levels
import astropy.utils.data
import ligo.skymap.plot
from scipy.stats import norm
import scipy.stats
#============================================================#
#    FUNCTION
#------------------------------------------------------------#
def sqsum(a, b):
	'''
	SQUARE SUM
	USEFUL TO CALC. ERROR
	'''
	import numpy as np
	return np.sqrt(a**2.+b**2.)
#------------------------------------------------------------#
def calc_MLR(s, colour, z):
	import numpy as np
	MLR        = 10.**(s*colour+z)
	return MLR
#------------------------------------------------------------#
def calc_MLR_error(s, colour, z, band0_err, band1_err, MLR_err):
	import numpy as np
	term0        = np.log(10)*s*(10**(s*colour+z))*band0_err
	term1        = np.log(10)*(-s)*(10**(s*colour+z))*band1_err
	MLR_err0    = sqsum(term0, term1)
	total_err    = sqsum(MLR_err0, MLR_err)
	total_err.name    = 'MLR_err'
	return total_err
#------------------------------------------------------------#
def calc_mass_bol(M_abs_gal, MLR):
	import numpy as np
	L_sun_bol        = 3.828*10.**26                    [W]
	M_abs_sun_bol    = 4.74                            [mag]
	mass_gal    = MLR*(10.**(-0.4*(M_abs_gal-M_abs_sun_bol)))
	#mass_gal.name    = 'stellar mass'
	print('RETURN STELLAR MASS [SOLAR MASS]\n'+'1.\t[kg]\n'+'2.\t[solar mass]\n'+'3.\t[log10([solar mass])]')
	return mass_gal, mass_gal/mass_sun, np.log10(mass_gal/mass_sun)
#------------------------------------------------------------#
def calc_mass(M_abs_gal, MLR, L_sun, M_abs_sun):
	import numpy as np
	L_sun_bol        = 3.828*10.**26                    #    [W]
	MLR_sun            = L_sun/L_sun_bol
	mass_gal        = MLR*MLR_sun*(10.**(-0.4*(M_abs_gal-M_abs_sun)))
	print('RETURN STELLAR MASS [SOLAR MASS]\n'+'1.\t[solar mass]\n'+'2.\t[log10([solar mass])]')
	return mass_gal, np.log10(mass_gal)
#------------------------------------------------------------#
'''
def calc_mass_error(M_abs_gal, MLR, L_sun, M_abs_sun, M_abs_gal_err, MLR_err):
	import numpy as np
	L_sun_bol        = 3.828*10.**26                    [W]
	MLR_sun            = L_sun/L_sun_bol
	term0            = MLR_sun*(10.**(-0.4*(M_abs_gal-M_abs_sun)))
	term1            = sqsum(MLR_err, MLR*(-0.4)*np.log(10)*M_abs_gal_err)
	total_err        = term0*term1
	total_err.name    = 'mass_err'
	return total_err, np.log10(total_err)
'''
#------------------------------------------------------------#
def calc_mass_error(M_abs_gal, MLR, L_sun, M_abs_sun, M_abs_gal_err, MLR_err):
	import numpy as np
	L_sun_bol        = 3.828*10.**26                    #    [W]
	MLR_sun            = L_sun/L_sun_bol
	total_err        = MLR_err*MLR_sun*(10.**(-0.4*(M_abs_gal-M_abs_sun)))
	total_err.name    = 'mass_err'
	return total_err, np.log10(total_err)
#------------------------------------------------------------#
def calc_lumi(M_abs):
	'''
	CALCULATE LUMINOSITY FOR DESINATED BAND WITH B BAND.
	'''
	import numpy as np
	L_sun_B        = 1.9*10.**26                    #    [W]
	M_abs_B        = 5.497                            #    [AB]
	L            = L_sun_B*(10.**((M_abs_B-M_abs)/2.5))
	return L
#------------------------------------------------------------#
def calc_mass_K(M_abs_gal, M_abs_sun):
	import numpy as np
	mass_gal    = 10.**((M_abs_sun-M_abs_gal)/2.5)
	return mass_gal
#------------------------------------------------------------#
def prob_fin(map_prob, P_2D, D_gal, D_gal_sig, D_gw, D_gw_sig):
	'''
	P_gal    =    (normalization factor)    *    (Mappeli+18 prob) *
				(P_2D from bayester)    *    (1-error fuction for distance)
	'''
	import numpy as np
	from scipy.special import erf
	term    = (1.-erf((np.abs(D_gal-D_gw))/(sqsum(D_gal_sig, D_gw_sig)**2)))
	P_gal    = map_prob * P_2D * term
	nP_gal    = (P_gal)/np.sum(P_gal)
	return nP_gal
#------------------------------------------------------------#
#   FUNCTION
def replaceNULL(x):
	if x == 'null':
		return -99
	return x
#--------------------------------------------------------------------------------#
def matching(inra, inde, refra, refde):
	"""
	MATCHING TWO CATALOG WITH RA, Dec COORD. WITH python
	INPUT   :   SE catalog, SDSS catalog file name, sepertation [arcsec]
	OUTPUT  :   MATCED CATALOG FILE & TABLE
	"""
	#	MODULE
	import numpy as np
	import astropy.units as u
	from astropy.table import Table, Column
	from astropy.coordinates import SkyCoord
	from astropy.io import ascii
	#comment		= '\nMATCHING START\t'+str(len(inra))+', '+str(len(refra));print(comment)
	#	SKY COORD.
	coo_intbl   = SkyCoord(inra, inde, unit=(u.deg, u.deg))
	coo_reftbl  = SkyCoord(refra, refde, unit=(u.deg, u.deg))
	#   INDEX FOR REF.TABLE
	indx, d2d, d3d  = coo_intbl.match_to_catalog_sky(coo_reftbl, nthneighbor=1)
	comment		= 'MATHCING END';print(comment)
	return indx, d2d
#--------------------------------------------------------------------------------#
def probcut(healpixfits, prob_cut=0.9):
	import healpy as hp
	import numpy as np
	from astropy.table import Table, Column, MaskedColumn
	comment     = '='*80+'\n' \
				+ 'Reading healpix file ...'
	print(comment)
	hpx, hdr    = hp.read_map(healpixfits, h=True, verbose=False)
	hdr         = dict(hdr)
	gwdist      = hdr['DISTMEAN']
	gwdiststd   = hdr['DISTSTD']
	#   h       : GIVE ME HEADER INFO
	#   verbose : DON'T PRINT INFO MESSAGE
	#   AFTER TYPING '%matplotlib' ON THE PYTHON PROMPT.
	#   SHOW MAP
	#   hp.mollview(hpx)
	npix        = hdr['NAXIS2']         # i.e. 50331648
	nside       = hp.npix2nside(npix)   # i.e. 2048
	ordering    = hdr['ORDERING']
	#   SORTING
	#   i.e.
	#   array([  1.72849272e-121,   1.72849272e-121,   1.72849272e-121, ...,
	#           6.41252562e-005,   6.41299517e-005,   6.41368262e-005])
	hpx_sort    = (-1)*np.sort((-1)*hpx)
	indx_sort   = np.argsort((-1)*hpx)
	stack       = np.cumsum(hpx_sort)
	#   MAKE PROB. CUT INDEX (SMALLER THAN PROB.) FOR SORTED DATA
	#   prob_cut    = 0.9
	indx_prob   = np.where(stack <= prob_cut)
	ipix_nparr  = np.linspace(0, npix, npix, dtype=int)
	#   NEED TO CONVERT NUMPY ARRAY TO LIST (WHY?)
	ipix        = ipix_nparr.tolist()
	#   PIXEL COORD. TO SPHERERICAL COORD.
	if ordering == 'NESTED':
		theta, phi  = hp.pix2ang(nside, ipix, nest=False)
	else :
		theta, phi  = hp.pix2ang(nside, ipix, nest=True)
	#   RADIAN TO DEGREE(RA & DEC FORM)
	ra, de      = np.rad2deg(phi), np.rad2deg(0.5*np.pi-theta)
	#   SORTING TO DECREASING ORDER
	rasort      = ra[indx_sort]
	desort      = de[indx_sort]
	#   EXPORT WITHIN PROB.
	racut       = rasort[indx_prob]
	decut       = desort[indx_prob]
	#   MAKE ASCII FILE: ra[deg], dec[deg]
	#   PUT # in front of column name because of running stilts
	outdata		= Table([racut, decut, hpx_sort[indx_prob]], \
						names=['ra', 'dec', 'prob'])	
	#   CALCULATE MEAN DISTANCE BETW. PIXELS
	alldeg		= 41253.         						#	[deg^2]
	mpix_dist   = np.sqrt(alldeg / npix)*3600.			#	["/pix]
	
	return outdata, mpix_dist*2., gwdist, gwdiststd
#--------------------------------------------------------------------------------#
def generate_target(result, dist_lower, dist_upper):
	from astropy.table import Table, Column
	from astropy.io import ascii
	indx_dist   = np.where(((dist_lower < result['dist']) & (dist_upper > result['dist'])))
	num_within  = len(indx_dist[0])
	comment     = 'Cadidates galaxies\t: '+str(num_within)
	print(comment)
	result_in   = result[indx_dist]
	#   CONVERT 'masked' to -99
	result_in   = result_in.filled(-99)
	ra_result   = result_in['ra']
	dec_result  = result_in['dec']
	name        = []
	ra          = []
	dec         = []
	prob        = []
	dist        = []
	disterr     = []
	Bmag        = []
	Bmagerr     = []
	Kmag        = []
	Kmagerr     = []
	W1          = []
	W1err       = []
	for i in range(len(result_in['GWGCname'])):
		if result_in['GWGCname'][i] != '-99':
			name.append(result_in['GWGCname'][i])
		else:
			if result_in['PGC'][i]  != '-99':
				name.append(result_in['PGC'][i])
			else:
				if result_in['HyperLEDAname'][i]    != '-99':
					name.append(result_in['HyperLEDAname'][i])
				else:
					if result_in['2MASSname'][i]    != '-99':
						name.append('2MASS+'+result_in['2MASSname'][i])
					else:
						name.append('SDSS+'+result_in['SDSS-DR12name'][i])
		coord_new   = SkyCoord(ra=ra_result[i]*u.degree, \
								dec=dec_result[i]*u.degree, frame='icrs')
		#   RA [deg] to RA [hh:mm:ss.s]
		ra_h        = int(coord_new.ra.hms[0])
		ra_m        = int(coord_new.ra.hms[1])
		ra_s        = int(coord_new.ra.hms[2])
		if ra_h < 10: ra_h_str      = '0'+str(ra_h)
		else: ra_h_str              = str(ra_h)     
		if ra_m < 10: ra_m_str      = '0'+str(ra_m)
		else: ra_m_str              = str(ra_m)     
		if ra_s < 10: ra_s_str      = '0'+str(ra_s)
		else: ra_s_str              = str(ra_s)     
		ra_new      = ra_h_str+':'+ra_m_str+':'+ra_s_str
		ra.append(ra_new)

		de_d        = int(coord_new.dec.dms[0])
		de_m        = np.abs(int(coord_new.dec.dms[1]))
		de_s        = round(np.abs(coord_new.dec.dms[2]), 1)
		if de_d > 0 :
			if de_d < 10:
				de_d_str            = '+0'+str(de_d)
			else:
				de_d_str            = '+'+str(de_d)
		else:
			de_d_str                = str(de_d)
		if de_m < 10: de_m_str      = '0'+str(de_m)
		else: de_m_str              = str(de_m)     
		if de_s < 10: de_s_str      = '0'+str(de_s)
		else: de_s_str              = str(de_s)     
		de_new      = de_d_str+':'+de_m_str+':'+de_s_str
		dec.append(de_new)

	name        = np.array(name)
	ra          = np.array(ra)
	dec         = np.array(dec)
	data        = Table([name, ra, dec], \
							names=['#name', 'ra', 'dec'])
	outdat      = 'gw_candidates_qv.dat'
	ascii.write(data, outdat, format='fixed_width', delimiter=' ')
	return result_in
#------------------------------------------------------------#
def addwise(result_in, raname='ra_1', dename='dec_1'):
	import numpy as np
	from astropy.table import Table, Column
	from astropy.io import ascii
	ra      = result_in[raname]
	dec     = result_in[dename]
	w1      = []
	w1sig   = []

	for i in range(len(ra)):
		table_  = Irsa.query_region(SkyCoord(ra[i], dec[i], \
									unit=(u.deg,u.deg), frame='icrs'), \
									catalog='allwise_p3as_psd', radius='0d1m0s')
		wra,wdec= table_['ra'], table_['dec']
		dist    = np.sqrt( (ra[i] - wra)**2. + (dec[i] -wdec)**2. )
		indx_min= np.where( dist == np.min(dist) )
		table   = table_[indx_min]
		try:
			print('add wise info ...')
			w1.append(table['w1mpro'][0])
			w1sig.append(table['w1sigmpro'][0])
		except:
			print('exception!')
			w1.append(-99)
			w1sig.append(-99)
	
	w1      = np.array(w1)
	w1sig   = np.array(w1sig)
	result_in['w1']     = w1
	result_in['w1sig']  = w1sig
	new_result          = result_in
	return new_result
#------------------------------------------------------------#
def app2abs(mag_app, dist, appkey='K', distkey='dist'):
	#appmag  = intbl[appkey]
	#dist    = intbl[distkey]
	import numpy as np
	abslist = []
	for i in range(len(mag_app)):
		if (mag_app[i] != -99) & (dist[i] != -99):
			mag_abs  = mag_app[i]-5.*np.log10(dist[i]*1e6)+5
		else:
			mag_abs = -99
		abslist.append(round(mag_abs, 3))
	return np.array(abslist)
#------------------------------------------------------------#
def prob(gwdist, gwdister, gldist, gldister, appmag, prob):
	'''
	Swope Supernova Survey 2017a (SSS17a),
	the Optical Counterpart to a Gravitational Wave Sourve
	'''
	import numpy as np
	from scipy import special
	P_gal   = []
	for i in range(len(prob)):
		if gldister[i] == 'null':
			gldister[i] = float(0)
		else:
			gldister[i]	= float(gldister[i])

		L   = ((gldist[i])**2.)*10.**(-0.4*appmag[i])
		braket  = np.abs(gldist[i]-gwdist)/((gwdister**2.)+(gldister[i]**2.))
		P       = L*prob[i]*(1.-special.erf(braket))

		P_gal.append(P)
	P_gal   = np.array(P_gal)
	nP_gal  = P_gal*(1./np.sum(P_gal))
	return nP_gal
#------------------------------------------------------------#
def generate_final_list(incat, outdat):
	import numpy as np
	from astropy.table import Table, Column
	from astropy.io import ascii
	from astropy.coordinates import SkyCoord
	import astropy.units as u
	num         = len(incat)
	try:
		indx_TheMost= np.where(incat['K_prob'] == np.max(incat['K_prob']))
	except:
		indx_TheMost= np.where(incat['B_prob'] == np.max(incat['B_prob']))
	comment     = '='*80+'\n' \
				+ 'Generating '+outdat+'\n' \
				+ 'Cadidates galaxies\t: '  +str(num)+'\n'+'-'*80+'\n' \
				+ 'The Highest Prob.\t: '   +str(incat['GWGCname'][indx_TheMost][0])
	print(comment)

	ra_result   = incat['ra']
	dec_result  = incat['dec']
	name        = []
	ra          = []
	dec         = []
	dist        = []
	Bmag        = []
	Kmag        = []

	for i in range(len(incat['GWGCname'])):
		if incat['GWGCname'][i] != '-99':
			name.append(incat['GWGCname'][i])
		else:
			if incat['PGC'][i]  != -99:
				name.append(incat['PGC'][i])
			else:
				if incat['HyperLEDAname'][i]    != '-99':
					name.append(incat['HyperLEDAname'][i])
				else:
					if incat['2MASSname'][i]    != '-99':
						name.append('2MASS+'+incat['2MASSname'][i])
					else:
						name.append('SDSS+'+incat['SDSS-DR12name'][i])
		coord_new   = SkyCoord(ra=ra_result[i]*u.degree, \
								dec=dec_result[i]*u.degree, frame='icrs')
		#   RA [deg] to RA [hh:mm:ss.s]
		ra_h        = int(coord_new.ra.hms[0])
		ra_m        = int(coord_new.ra.hms[1])
		ra_s        = int(coord_new.ra.hms[2])
		if ra_h < 10: ra_h_str      = '0'+str(ra_h)
		else: ra_h_str              = str(ra_h)     
		if ra_m < 10: ra_m_str      = '0'+str(ra_m)
		else: ra_m_str              = str(ra_m)     
		if ra_s < 10: ra_s_str      = '0'+str(ra_s)
		else: ra_s_str              = str(ra_s)     
		ra_new      = ra_h_str+':'+ra_m_str+':'+ra_s_str
		ra.append(ra_new)

		de_d        = int(coord_new.dec.dms[0])
		de_m        = np.abs(int(coord_new.dec.dms[1]))
		de_s        = round(np.abs(coord_new.dec.dms[2]), 1)
		if de_d > 0 :
			if de_d < 10:
				de_d_str            = '+0'+str(de_d)
			else:
				de_d_str            = '+'+str(de_d)
		else:
			de_d_str                = str(de_d)
		if de_m < 10: de_m_str      = '0'+str(de_m)
		else: de_m_str              = str(de_m)     
		if de_s < 10: de_s_str      = '0'+str(de_s)
		else: de_s_str              = str(de_s)     
		de_new      = de_d_str+':'+de_m_str+':'+de_s_str
		dec.append(de_new)
		dist.append(round(incat['dist'][i], 3))
		Bmag.append(round(incat['B'][i], 3))
		Kmag.append(round(incat['K'][i], 3))
	name        = np.array(name)
	ra          = np.array(ra)
	dec         = np.array(dec)
	dist        = np.array(dist)
	Bmag        = np.array(Bmag)
	Kmag        = np.array(Kmag)
	data        = Table([name, ra, dec, dist, Bmag, Kmag], names=['name', 'ra', 'dec', 'dist', 'Bmag', 'Kmag'])
	# outdat      = 'gw_candidates_sorted.dat'
	ascii.write(data, outdat, format='fixed_width', delimiter=' ')
	outcat      = ascii.read(outdat)
	return outcat
#------------------------------------------------------------#
def heal2table(healpixfits, confidence=0.9, eventname='GW_signal', save_path='./', hdr=True, view=True):
	import numpy as np
	import healpy as hp
	from astropy.table import Table
	hpx, hdr    = hp.read_map(healpixfits, h=True, verbose=False)
	hdr         = dict(hdr)

	gwdist      = hdr['DISTMEAN']
	gwdiststd   = hdr['DISTSTD']
	#   h       : GIVE ME HEADER INFO
	#   verbose : DON'T PRINT INFO MESSAGE
	#   AFTER TYPING '%matplotlib' ON THE PYTHON PROMPT.
	npix        = hdr['NAXIS2']         # i.e. 50331648
	nside       = hp.npix2nside(npix)   # i.e. 2048
	ordering    = hdr['ORDERING']
	#   SORTING
	#   i.e.
	#   array([  1.72849272e-121,   1.72849272e-121,   1.72849272e-121, ...,
	#           6.41252562e-005,   6.41299517e-005,   6.41368262e-005])
	hpx_sort    = (-1)*np.sort((-1)*hpx)
	indx_sort   = np.argsort((-1)*hpx)
	stack       = np.cumsum(hpx_sort)
	#   MAKE PROB. CUT INDEX (SMALLER THAN PROB.) FOR SORTED DATA
	#   confidence    = 0.9
	indx_prob   = np.where(stack <= confidence)
	ipix_nparr  = np.linspace(0, npix, npix, dtype=int)
	#   NEED TO CONVERT NUMPY ARRAY TO LIST (WHY?)
	ipix        = ipix_nparr.tolist()
	#   PIXEL COORD. TO SPHERERICAL COORD.
	if ordering == 'NESTED':
		theta, phi  = hp.pix2ang(nside, ipix, nest=False)
	else :
		theta, phi  = hp.pix2ang(nside, ipix, nest=True)
	#   RADIAN TO DEGREE(RA & DEC FORM)
	ra, de      = np.rad2deg(phi), np.rad2deg(0.5*np.pi-theta)
	#   SORTING TO DECREASING ORDER
	rasort      = ra[indx_sort]
	desort      = de[indx_sort]
	#   EXPORT WITHIN PROB.
	racut       = rasort[indx_prob]
	decut       = desort[indx_prob]
	#   MAKE ASCII FILE: ra[deg], dec[deg]
	#   PUT # in front of column name because of running stilts
	healtbl		= Table([racut, decut, hpx_sort[indx_prob]], names=['ra', 'dec', 'P_2D'])
	#   SHOW MAP
	if view == True:
		hp.mollview(hpx, title=eventname)
		plt.savefig(save_path+'/'+eventname+'-skymap.png')
	return healtbl, hdr
#------------------------------------------------------------
def kilonova_mag(gwdist, gwdiststd):
	import numpy as np
	m0		= 17.476	#	[AB] in i-band (t0+10h)
	m0err	= 0.018	
	dist0	= 38.4		#	[MPC]	Im et al. 2017
	dist0err= 8.9
	
	m		= m0+5.*np.log10(gwdist/dist0)
	merr	= np.sqrt( (m0err)**2 + ((5.*gwdiststd)/(gwdist*np.log(10)))**2 + ((5.*dist0err)/(dist0*np.log(10)))**2 )
	return m, merr
#------------------------------------------------------------
def expectedLC(eventname, hdr, save_path):
	import os, glob
	import numpy as np
	from astropy.io import ascii
	from astropy.table import Table, vstack
	import matplotlib.pyplot as plt
	from astropy.time import Time
	from imsng import tool, gw
	path_base	= '/home/gw/Research'
	t0			= Time(hdr['DATE-OBS'], format='isot', scale='utc')
	jd0_170817	= 2457983.02852
	droutbl		= ascii.read(path_base+'/phot_gw170817_Drout.dat')
	pirotbl		= ascii.read(path_base+'/cocoon_Piro+18.dat')
	#------------------------------------------------------------
	#	COCOON MODEL
	prtbl		= pirotbl[pirotbl['filter']=='r']
	prtbl['mag'], prtbl['magerr']	=	tool.abs2app(prtbl['absmag'], 0, hdr['DISTMEAN']*1e6, hdr['DISTSTD']*1e6)
	prtbl['mag']= prtbl['mag']-0.7
	#	GW170817-like
	rtbl		= droutbl[droutbl['filter']=='r']; rtbl['delmjd'].sort()
	rmag, rmagerr= tool.calc_app(rtbl['mag'], rtbl['magerr'], 38.4, 8.9, hdr['DISTMEAN'], hdr['DISTSTD'])
	#------------------------------------------------------------
	#	PLOT 1	: TIME - MAG.	(r-band)
	#------------------------------------------------------------
	plt.close('all')
	plt.rcParams.update({'font.size': 16})
	fig, ax0	= plt.subplots(nrows=1, ncols=1, sharey=False, figsize=(9, 9))
	#------------------------------------------------------------
	#	GW170817-like
	ax0.plot(rtbl['delmjd'], rmag, color='red', alpha=0.5, label='GW170817-like')
	ax0.fill_between(rtbl['delmjd'], rmag-rmagerr, rmag+rmagerr, color='tomato', alpha=0.15, label='_nolegend_')
	#	COCOON MODEL (Piro+2018)
	ax0.plot(prtbl['delmjd'], prtbl['mag'], color='dodgerblue', alpha=0.5, label='Shock Cooling')
	ax0.fill_between(prtbl['delmjd'], prtbl['mag']-prtbl['magerr'], prtbl['mag']+prtbl['magerr'], color='dodgerblue', alpha=0.3, label='_nolegend_')
	#------------------------------------------------------------
	#	SETTING
	#------------------------------------------------------------
	ax0.set(xlabel='Time (Days from merger)', ylabel=r'Magnitude')
	ax0.set_ylim([int(np.max(rmag))+0.5, int(np.min(prtbl['mag']))-0.5])
	ax0.set_xlim([0,2])
	plt.axvline(x=0.48, color='grey', linewidth=2, linestyle='--', label='GW170817 EM discovery')
	ax0.legend(loc='upper right', prop={'size':20})
	plt.title('{0} r-band'.format(eventname))
	plt.tight_layout()
	plt.minorticks_on()
	plt.savefig('{0}/{1}_LC_rband.png'.format(save_path, eventname), overwrite=True)

#------------------------------------------------------------
def read3Dhealpix2candidates(path_healpix, path_catalog, eventname='GW_signal', conflist=[0.5, 0.9], distsigcut=3, header=True):
	hpx, hdr	= hp.read_map(path_healpix, verbose=True, h=True)
	hdr			= dict(hdr)
	prob, distmu, distsigma, distnorm = hp.read_map(path_healpix,
													field=range(4),
													dtype=('f8', 'f8', 'f8', 'f8'))
	npix		= len(prob)
	nside		= hp.npix2nside(npix)
	pixarea		= hp.nside2pixarea(nside)		
	pixarea_deg2= hp.nside2pixarea(nside, degrees=True)
	#------------------------------------------------------------
	gldtbl0		= Table.read(path_catalog, format='ascii')
	gldtbl		= gldtbl0[	(gldtbl0['dist']<=hdr['DISTMEAN']+distsigcut*hdr['DISTSTD'])&
							(gldtbl0['dist']>=hdr['DISTMEAN']-distsigcut*hdr['DISTSTD'])]
	gldcoord	= SkyCoord(ra=gldtbl['ra']*u.deg, dec=gldtbl['dec']*u.deg)
	#	VOID COLUMN IN CATALOG TABLE
	ngld		= np.size(gldtbl)
	probdencol	= Column(np.zeros(ngld, dtype='f4'), name='dP_dV')
	probcol		= Column(np.zeros(ngld, dtype='f4'), name='P')
	probdenAcol	= Column(np.zeros(ngld, dtype='f4'), name='dP_dA')
	probAcol	= Column(np.zeros(ngld, dtype='f4'), name='P_A')
	gldtbl.add_columns([probdencol, probcol, probdenAcol, probAcol])
	#------------------------------------------------------------
	#	CALC hp INDEX FOR EACH GALAXY
	theta	= 0.5 * np.pi - gldcoord.dec.to('rad').value
	phi		= gldcoord.ra.to('rad').value
	ipix	= hp.ang2pix(nside, theta, phi)
	cumP2D	= np.cumsum(prob[np.argsort(-1*prob)])[ipix]
	#	CALC. PROBABILITY (P_2D)
	dp_dA		= prob[ipix]/pixarea
	dp_dA_deg2	= prob[ipix]/pixarea_deg2
	#	CALC. PROBABILITY DENSITY PER VOLUME (P_3D)
	dp_dV	= prob[ipix] * distnorm[ipix] * norm(distmu[ipix],distsigma[ipix]).pdf(gldtbl['dist'])/pixarea 
	gldtbl['dP_dV']	= dp_dV
	gldtbl['dP_dA']	= dp_dA
	#------------------------------------------------------------
	#	CALC. SCORE BASED ON POSITION AND P_3D
	#------------------------------------------------------------
	gldtbl['Prob']	= gldtbl['dist']**2 * 10**(-0.4*gldtbl['K']) * gldtbl['dP_dV']
	#------------------------------------------------------------
	cantbl			= gldtbl[	(gldtbl['K']!=-99.0)&
								(gldtbl['dist']!=-99.0)&
								(gldtbl['Prob']!=0.0)]
	cantbl['Prob']	= cantbl['Prob']/np.sum(cantbl['Prob'])
	cantbl.sort('Prob')
	cantbl.reverse()
	cantbl.meta['event'] = eventname
	cantbl.meta['distcut'] = distsigcut
	cantbl.meta['path_healpix'] = path_healpix
	cantbl.meta['path_catalog'] = path_catalog
	#------------------------------------------------------------
	confareainfo = dict()
	credible_levels = find_greedy_credible_levels(prob)
	for conf in conflist:
		areainconf = np.sum(credible_levels <= conf) * hp.nside2pixarea(nside, degrees=True)
		confareainfo[str(conf)] = areainconf	
	if header == True:
		return cantbl, prob, confareainfo, hdr
	else:
		return cantbl, prob, confareainfo
#------------------------------------------------------------
def plotcumscore(cantbl, probcutlist=[0.5, 0.90, 0.95, 0.99], eventname='GW_signal', path_save='.'):
	plt.close('all')
	plt.plot(np.arange(len(cantbl)), np.cumsum(cantbl['Prob']), 'dodgerblue', label='All({})'.format(len(cantbl)))
	# plt.scatter(np.arange(len(cnatbl)), np.cumsum(cantbl['Prob']), color='tomato', marker='+')
	# for probcut in [0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.99]:
	for probcut in probcutlist:
		subtbl = cantbl[np.cumsum(cantbl['Prob'])<=probcut]
		print('PROB.CUT\t{} : {}'.format(probcut, len(subtbl)))
		plt.axhline(y=probcut, color='tomato', alpha=0.5, linestyle='--', label='{}({})'.format(probcut, len(subtbl)))
	plt.xlim([-5, 2*len(subtbl)])
	plt.xlabel('Cumulative Score', fontsize=15)
	plt.ylabel('Number of objects', fontsize=15)
	plt.xticks(size=15)
	plt.yticks(np.arange(0, 1.1, 0.1), size=15)
	plt.legend(fontsize=15, loc='lower right')
	plt.minorticks_on()
	plt.title(eventname, size=20)
	plt.tight_layout()
	plt.savefig(path_save+'/'+eventname+'-cumscore.png', overwrite=True)
