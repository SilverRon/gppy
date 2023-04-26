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
# from imsng import phot
import sys
sys.path.append('..')
from phot import gpphot
import time
#============================================================
#	FUNCTION
#============================================================
def extractonerow(alltbl, number, apertures):
	onerow = alltbl[alltbl['NUMBER']==number]
	if len(onerow) > 1:
		onerow = Table(onerow[0])
	# onerow = alltbl[0]
	mag, mager = [], []
	flux, fluxer = [], []
	aper = apertures
	snr = []
	# for i in range(len(aper)):
	for i in range(0, 32):
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
	# x, y = aper[1:].value, snr[1:]
	x, y = aper[1:].value, snr[1:]
	# print(len(x), len(y))
	# print(x, y)
	# plt.scatter(x, y)
	# plt.plot(x, y, marker='o', linestyle='')
	plt.plot(x, y, '-', color='grey', alpha=0.5)
	optaper = x[y == np.max(y)]
	plt.axvline(x=optaper)

	# plt.errorbar(aper[1:], magdif, yerr=magdifer,
				# color='grey', alpha=0.5,
				# marker='o', linestyle='')
	
	return magdif, magdifer, optaper
#------------------------------------------------------------
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
#------------------------------------------------------------
def gcurveplot(inim, alltbl, apertures, pixelscale):
	plt.close('all')
	plt.rc('font', family='serif')
	fig = plt.figure()
	# x = 1920 / 4 / fig.dpi
	# y = 1080 / 4 / fig.dpi
	x = 1920 / 2 / fig.dpi
	y = 1080 / 2 / fig.dpi
	fig.set_figwidth(x)
	fig.set_figheight(y)

	optapers = []
	for number in set(alltbl['NUMBER']):
		magdif, magdifer, optaper = extractonerow(alltbl, number, apertures)
		optapers.append(optaper.item())

	# y, yer = meandelmagdif(alltbl)

	# plt.plot(aper[1:], yer, color='dodgerblue')
	# plt.axvline(aper[1:][yer == np.min(yer)], color='tomato', linestyle='--')

	optaper = round(np.mean(optapers), 2)
	plt.axvline(x=optaper, color='tomato', label="""{} pix ({} ")""".format(optaper, round((optaper/pixelscale).value, 2)), alpha=0.5)
	# print(optaper)

	down, up = plt.ylim()
	left, right = plt.xlim()
	plt.xlim([0, right])
	plt.ylim([0, up])
	plt.legend(fontsize=10, loc='best', framealpha=1, markerscale=2)
	plt.tick_params(which='both', direction='in')
	plt.xticks(fontsize=10)
	plt.yticks(fontsize=10)
	# plt.gca().invert_yaxis()
	plt.title(os.path.basename(inim), fontsize=10)
	plt.xlabel(r'Diameter [pixel]', fontsize=10)
	plt.ylabel(r'SNR', fontsize=10)
	plt.grid(color='dimgrey', linestyle='--', linewidth=1, alpha=0.5)
	plt.tight_layout()
	plt.minorticks_on()

	# plt.savefig('{}.gcurve.png'.format(inim[:-5]), overwrite=True)
	plt.savefig('{}.gcurve.png'.format(inim[:-5]),)

	return optaper, optapers
