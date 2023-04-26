#	USEFUL FUNC.(TOOL) IN IMSNG MODULE
#	CREATED IN 19.03.03 BY GREGORY S.H. PAEK
#	UPDATE : 20.01.03
#============================================================
#	MODULE
#------------------------------------------------------------
import os, sys, glob
import numpy as np
from astropy.io import ascii, fits
import astropy.coordinates as coord
import astropy.units as u
from multiprocessing import Process, Pool
import multiprocessing as mp
import time
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
import matplotlib.pyplot as plt
sys.path.append('..')
from phot import gpphot
# from imsng import phot_tbd
#============================================================
def timename():
	'''
	CONVERT 'TIME' TO YYMMDD, HHMMSS FORM.
	INPUT	:	NONE
	OUTPUT	:	STRIG FORM OF 'YYMMDD', 'HHMMSS'
	'''
	import numpy as np
	import time
	now				= time.gmtime(time.time())
	y, m, d			= now.tm_year, now.tm_mon, now.tm_mday
	ho, mi, se		= now.tm_hour, now.tm_min, now.tm_sec
	yy				= str(y)[2:]
	if len(str(m)) < 2:
		mm			= '0'+str(m)
	else:
		mm			= str(m)
	if len(str(d)) < 2:
		dd			= '0'+str(d)
	else:
		dd			= str(d)
	if len(str(ho)) < 2:
		hour		= '0'+str(ho)
	else:
		hour		= str(ho)
	if len(str(mi)) < 2:
		mini		= '0'+str(mi)
	else:
		mini		= str(mi)
	if len(str(se)) < 2:
		sec			= '0'+str(se)
	else:
		sec			= str(se)
	yymmdd			= yy+mm+dd
	hhmmss			= hour+mini+sec
	return yymmdd, hhmmss
#------------------------------------------------------------
def detection(name, ra, dec, time, location):
	import numpy as np
	import os, glob, sys
	from astropy import units as u
	from astropy.time import Time
	from astropy.coordinates import SkyCoord, EarthLocation, AltAz
	from astropy.coordinates import get_sun, get_moon
	from astropy.io import ascii
	from astropy.table import Table, Column
	target      = SkyCoord(ra, dec, unit='deg') # defaults to ICRS frame
	site        = location
	del_midnight= np.linspace(-12, +12, 720) * u.hour
	time_night  = time+del_midnight
	frame_night = AltAz(obstime=time_night, location=site)

	targetaltaz_night   = target.transform_to(frame_night)
	sunaltaz_night      = get_sun(time_night).transform_to(frame_night)

	# indx_set            = np.where( sunaltaz_night.alt > -18 * u.deg )
	indx_rise           = np.where( sunaltaz_night.alt < -18 * u.deg )
	
	sunset              = del_midnight[np.min(indx_rise)]
	sunrise             = del_midnight[np.max(indx_rise)]
	
	del_midnight= np.linspace(sunset.value, sunrise.value, 100) * u.hour
	time_night  = time+del_midnight
	frame_night = AltAz(obstime=time_night, location=site)
	targetaltaz_night   = target.transform_to(frame_night)

	return targetaltaz_night
#------------------------------------------------------------
def rts_maker(filename, savepath, obs, obstbl, intbl, date, hhmmss):
	from astropy.coordinates import Angle
	import numpy as np
	import os, glob, sys
	from astropy import units as u
	from astropy.time import Time
	from astropy.coordinates import SkyCoord, EarthLocation, AltAz
	from astropy.coordinates import get_sun, get_moon
	from astropy.io import ascii
	from astropy.table import Table, Column

	indx_obs	= np.where(obstbl['obs'] == obs)
	lat, lon, height = obstbl['lat'][indx_obs], obstbl['lon'][indx_obs], obstbl['height'][indx_obs]
	utoff, ul	= obstbl['utoff'][indx_obs], obstbl['ul'][indx_obs]

	lat		    = lat		* u.deg		# North
	lon		    = lon		* u.deg		# East
	height		= height	* u.m
	utoff		= utoff		* u.hour
	ul			= 20					# limiting magnitude

	location    = EarthLocation(lat=lat, lon=lon, height=height)

	time        = Time('20'+date[0:2]+'-'+date[2:4]+'-'+date[4:6]+' 00:00:00')+utoff
	risetime	= []
	maxtime		= []
	settime		= []
	namelist	= []
	ralist		= []
	delist		= []

	for i in range(len(intbl)):
		name, ra, dec		= intbl['name'][i], intbl['ra'][i], intbl['dec'][i]
		targetaltaz_night	= detection(name, ra, dec, time, location)
		try:
			alt_max             = np.max(targetaltaz_night.alt)
			alt_max_time		= targetaltaz_night.obstime[targetaltaz_night.alt == alt_max] + utoff
			alt_rise30_time		= targetaltaz_night.obstime[targetaltaz_night.alt >= 25 * u.deg][0] + utoff
			alt_set30_time		= targetaltaz_night.obstime[targetaltaz_night.alt <= 0 * u.deg][0] + utoff

			if alt_max >= 30.0 *u.deg:

				risetime.append(alt_rise30_time.value[0][11:])
				maxtime.append(alt_max_time.value[0][11:])
				settime.append(alt_set30_time.value[0][11:])
				namelist.append(name)
				ralist.append(Angle(str(ra)+'d').to_string(unit=u.hour, sep=':'))
				delist.append(Angle(str(dec)+'d').to_string(unit=u.degree, sep=':'))
		except:
			pass

	risetime	= np.array(risetime)
	maxtime		= np.array(maxtime)
	settime		= np.array(settime)
	namelist	= np.array(namelist)
	ralist		= np.array(ralist)
	delist		= np.array(delist)

	targettbl	= Table(	[namelist, ralist, delist, risetime, maxtime, settime], names=['name', 'ra', 'dec', 'rise', 'transit', 'set'])

	ascii.write(	targettbl,
					output=savepath+date+'/'+date+'-'+hhmmss+'-targetlist-'+obs+'-'+filename+'.txt',
					format='fixed_width',
					delimiter=None,
					overwrite=True)
	'''
	ascii.write(	targettbl,
					output=savepath+date+'/'+date+'-'+hhmmss+'-targetlist-'+obs+'-'+filename+'.txt',
					delimiter=None,
					overwrite=True)
	'''
#------------------------------------------------------------
def sendmail(filename, subject, sendID, sendPW, reciver):
	'''
	Security reference
	https://cpuu.postype.com/post/23066
	Code reference
	https://kimdoky.github.io/python/2017/07/21/smtplib_email.html
	File attach
	https://brunch.co.kr/@jk-lab/31
	'''
	import smtplib
	from email.mime.text import MIMEText
	import codecs
	email_text = codecs.open(filename, 'rb', 'utf-8')
	msg = MIMEText(email_text.read())
	email_text.close()
	msg['Subject']	= subject
	msg['From']		= sendID
	smtp_gmail		= smtplib.SMTP_SSL('smtp.gmail.com', 465)
	smtp_gmail.login(sendID, sendPW)
	smtp_gmail.sendmail(sendID, reciver, msg.as_string())
	smtp_gmail.quit()
	comment	= 'Send '+filename+'\n'+'From '+sendID+' To '+reciver; print(comment)
#------------------------------------------------------------
def send_gmail(subject, contents, fromID, fromPW, toIDs, ccIDs=None, path=None):
	'''
	SEND GMAIL
	Security reference
	https://cpuu.postype.com/post/23066
	Code reference
	https://kimdoky.github.io/python/2017/07/21/smtplib_email.html
	File attach
	https://brunch.co.kr/@jk-lab/31
	'''
	import os
	import smtplib
	from email.mime.base import MIMEBase
	from email.mime.text import MIMEText
	from email.mime.multipart import MIMEMultipart
	from email.header import Header  
	#msg		= MIMEBase('mixed')
	#msg		= MIMEText(contents, 'plain', 'utf-8')
	msg		= MIMEMultipart()
	msg['Subject']	= Header(s=subject, charset="utf-8")
	msg['From']		= fromID
	msg['To']		= toIDs
	if ccIDs != None:
		msg['Cc']		= ccIDs
	msg.attach(MIMEText(contents, 'plain', 'utf-8'))
	#	ATTACH TEXT FILE ON MAIL
	if path != None:
		if type(path) != list:
			filelist	= []
			filelist.append(path)
		else:
			filelist	= path

		for file in filelist:
			part	= MIMEBase("application", "octet-stream")
			part.set_payload(open(file, 'rb').read())
			part.add_header(	'Content-Disposition',
								'attachment; filename="%s"'% os.path.basename(file))
			msg.attach(part)
	


	#	ACCESS TO GMAIL & SEND MAIL
	smtp_gmail		= smtplib.SMTP_SSL('smtp.gmail.com', 465)
	smtp_gmail.login(fromID, fromPW)
	smtp_gmail.sendmail(msg["From"], msg["To"].split(",") + msg["Cc"].split(","), msg.as_string())
	smtp_gmail.quit()
	comment	= 'Send '+str(path)+'\nFrom\t'+fromID+'\nTo'; print(comment); print(toIDs)
#------------------------------------------------------------
def calc_rts(filename, observatory, obsdate, dts, obsinfofile, catname, altlimit=30., moonseperation=30.):
	#cal_visibility.py - Changsu Choi, 2015/02/01
	#calculate rise transit set time and moon distance for observation of IMSNG
	#to be changed : adding observatory, bar plot, adding object, calculating for all nights of 2015  

	# Usage : python cal_visibility.py obsdate dts
	# python cal_visibility.py 2015/03/21 1
	# python 3 ported, 2019-03-05, changsu choi

	import mskpy.observing as obs
	from astropy.coordinates import Angle
	import astropy.units as u
	from astropy.io import ascii
	from astropy.table import Table, vstack
	import ephem
	import numpy as np
	import string
	import os, sys
	import astropy.coordinates as coord


	#	altitite and moon seperation parameter
	#moon serperation is a little bit close (2~3 deg)
	#altlimit		= 25.
	#moonseperation	= 30.
	#observatory		= 'LOAO'
	#obsdate			= '2019/10/06'
	#dts				= '0'

	observatory		= str(observatory)
	obsdate			= str(obsdate)
	dts				= str(dts)
	moonsepcut		= 360.-moonseperation
	
	#obs_info		= ascii.read("obs_info.dat")
	obs_info		= ascii.read(obsinfofile)

	if type(catname)==str	: tdata	= ascii.read(catname)
	else					: tdata = catname
	#catname			= 'targetlist_test.dat'
	
	#	Obseravatory information
	indx_obs		= np.where(observatory == obs_info['obs'])
	obsname			= obs_info['obs'][indx_obs][0]
	obslat			= obs_info['lat'][indx_obs][0]
	obslon			= obs_info['lon'][indx_obs][0]
	obstz			= obs_info['utoff'][indx_obs][0]

	observ			= ephem.Observer()
	observ.date		= obsdate+' 01:00:00'
	observ.lon		= str(obslon)
	observ.lat		= str(obslat)
	observ.elevation= obs_info['height'][indx_obs][0]

	#	Day Time Saving
	if int(dts)	==0:
		#print ('No day Time saving')
		dts	= float(dts)
	else:
		#print ('Ok then I will plus 1 hr to local time.')
		dts	= float(dts)

	#	objects from catalog file
	objname			= tdata['name']
	ra				= tdata['ra']
	dec				= tdata['dec']
	prior			= tdata['sort']

	radd			= ra
	decdd			= dec

	#	Moon distance and information
	mcoord			= ephem.Moon()
	mcoord.compute(obsdate)
	#print ('Moon ra, dec \n')
	mheader			='Moon ra, dec'
	#print (mcoord.ra,mcoord.dec,'\n')
	minfo			= mheader+' '+str(mcoord.ra)+' '+str(mcoord.dec)+'\n'
	mphase			= ephem.Moon(obsdate+' 00:00:00')
	#print ('Moon phase : '+ "%.2f" % mphase.moon_phase)
	mphasestr		='Moon phase : '+ "%.2f" % mphase.moon_phase +'\n'

	#	Angular distance calculation
	def angsep(ra1deg, dec1deg, ra2deg, dec2deg) : 
		ra1rad		= ra1deg*np.pi/180
		dec1rad		= dec1deg*np.pi/180
		ra2rad		= ra2deg*np.pi/180
		dec2rad		= dec2deg*np.pi/180
		cos_a		= np.sin(dec1rad)*np.sin(dec2rad)+(np.cos(dec1rad)*np.cos(dec2rad)*np.cos(ra1rad-ra2rad))
		anglesep	= np.arccos(cos_a)*180/np.pi
		return anglesep

	'''
	targets			= []
	targets.append(rad)
	targets.append(decd)
	targets.append(objname)
	'''
	msep			= angsep(radd,decdd, np.degrees(mcoord.ra), np.degrees(mcoord.dec))

	#sunrise calculation
	observ.horizon	= '-18'
	sunrise			= observ.next_rising(ephem.Sun())
	sunset			= observ.previous_setting(ephem.Sun())

	aaa				= ephem.Date.tuple(sunset)
	#hrr				= int(aaa[3]+obstz+dts+24)
	hrr				= int(aaa[3]+obstz+dts)
	mrr				= aaa[4]
	#print ('sunset : '+str(hrr)+':'+str(mrr))
	sunsetstr		= '-18 deg sunset : '+str(int(hrr))+':'+str(mrr)+'\n'
	sunseti			= hrr + mrr/60. + 0.25

	bbb				= ephem.Date.tuple(sunrise)
	hrr				= bbb[3]+obstz+dts
	mrr				= bbb[4]
	#print ('sunrise : '+str(int(hrr))+':'+str(mrr))
	sunrisestr		= '-18 deg sunrise : '+str(int(hrr))+':'+str(mrr)+'\n'
	sunriseti		= hrr + mrr/60. -0.25

	#f				= open("rts_vis_"+obsdate[0:4]+obsdate[5:7]+obsdate[8:10]+"_"+observatory+".txt",'w')
	f				= open(filename,'w')

	#g				= open("targets.data",'w')

	#header			= '{:25s}	{:12s}	{:10s}	{:5s}	{:5s}	{:5s}	{:5s}	{:1s}'.format('name', 'ra', 'dec', 'rise(LT)', 'transit(LT)', 'set(LT)', 'moon_dist(deg)', 'sort')+'\n'
	header			= 'name ra dec rise(LT) transit(LT) set(LT) moon_dist(deg) sort \n'

	dashline		= '#'+'-'*60+'\n'
	f.write(obsdate)
	f.write('\nobservatory = '+observatory+'\n')
	f.write(sunsetstr)
	f.write(sunrisestr)
	f.write(minfo)
	f.write(mphasestr)
	f.write('alt limit = '+str(altlimit)+'\n')
	f.write('Moon seeperation = '+str(moonseperation)+'\n')
	f.write(dashline)
	f.write(header)

	pobj			= []
	prt				= []
	ptt				= []
	pst				= []

	telescope		= obs.Observer(obslon*u.deg, obslat*u.deg, dts+obstz, obsdate, observatory)

	for n in range(len(objname)):
		ra_hms		= Angle(str(ra[n])+'d').to_string(unit=u.hour, sep=':')[:-2]
		de_hms		= Angle(str(dec[n])+'d').to_string(unit=u.deg, sep=':')[:-3]
		#	35 deg altitute cut
		rtscal		= obs.rts(radd[n], decdd[n], obsdate, obslon, obslat, float(obstz)+dts, limit=altlimit, precision=1440)
		rt			= rtscal[0]
		tt			= rtscal[1]
		st			= rtscal[2]

		if rtscal[0]==None:
			#print (objname[n], ra_hms, de_hms, rtscal[0], rtscal[1], rtscal[2],"%.2f" % msep[n])
			vis=objname[n]+' '+ra_hms+' '+de_hms+ ' '+str(rtscal[0])+' '+ str(rtscal[1])+' '+ str(rtscal[2])+' '+str(int(msep[n]))+ str(prior[n])+'\n'	
			#f.write(vis)
			#print(vis)
		
		elif sunriseti < rtscal[0] < sunseti and sunriseti < rtscal[2] < sunseti and sunriseti < rtscal[1] < sunseti : 
			#print ('It can be seen in daytime!')
			pass
		#	moon seperation = 50 deg cut	
		elif msep[n] < moonseperation :
			#print (objname[n]+' too close to Moon < '+str(moonseperation)+' deg')
			pass
		elif msep[n] > moonsepcut :
			#print (objname[n]+' is close to the Moon by ',str(360-msep[n])+' deg') 		 	
			pass
		else:
			rtp		= "%.2d" % int(rt)+':'+"%.2d" % int((rt-int(rt))*60)
			ttp		= "%.2d" % int(tt)+':'+"%.2d" % int((tt-int(tt))*60)
			stp		= "%.2d" % int(st)+':'+"%.2d" % int((st-int(st))*60)
			vis		= '{:25s}	{:12s}	{:10s}	{:5s}	{:5s}	{:5s}	{:5s}	{:1s}'.format(objname[n], ra_hms, de_hms, rtp, ttp, stp, str(int(msep[n])), str(prior[n]))+'\n'
			f.write(vis)
			#print(vis)
			#targets	= objname[n]+' , '+ra_hms+' hr , '+de_hms+' deg \n'
			#g.write(targets)
			#print (objname[n], ra_hms, de_hms, rtp, ttp, stp, "%.2f" % msep[n])
		
	f.close()
	#g.close()
	#os.system('pluma '+"rts_vis_"+obsdate[0:4]+obsdate[5:7]+obsdate[8:10]+"_loao.txt &")
	#targetfile		="rts_vis_"+obsdate[0:4]+obsdate[5:7]+obsdate[8:10]+".txt"
#------------------------------------------------------------
def ds9regmaker(name, ra, dec, radius=1.5, color='green', dashlist=8, width=2, font='helvetica', fontsize=10, filename='ds9.reg'):
	'''
	name = intbl['Object Name']
	ra = intbl['RA']
	dec = intbl['DEC']
	radius=1.5
	color='green'
	dashlist=8
	width=2
	font='helvetica'
	fontsize=10
	filename='ds9.reg'
	'''

	c = SkyCoord(ra, dec, unit='deg')
	hmsdms = c.to_string('hmsdms')
	hms, dms = [], []

	for hd in hmsdms:
		rahms = hd.split(' ')[0].replace('h', ':').replace('m', ':').replace('s', '')
		dedms = hd.split(' ')[1].replace('d', ':').replace('m', ':').replace('s', '')
		hms.append(rahms)
		dms.append(dedms)

	f = open(filename, 'w')

	head = """# Region file format: DS9 version 4.1\nglobal color={} dashlist={} 3 width={} font="{} {} normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5
	""".format(color, dashlist, width, font, fontsize)

	f.write(head)

	for r,d,n0 in zip(hms, dms, name):
		n = n0.replace('<a HREF="javascript:void(0)">', '').replace('</a>', '')
		reg = """circle({},{},{}") # text={}{}{}
		""".format(r, d, format(radius, ".3f"), '{', n, '}')
		f.write(reg)
	f.close()
#------------------------------------------------------------
def rtsmaker(observatory, headname, save_path, obspath, catpath, start, end, altlimit=30., moonseperation=40., sunlimit='-18', numlimit=100):
	import pytz
	import jdcal
	import ephem
	#from numpy import *
	import numpy as np
	import os, sys
	import string
	import datetime
	import astropy.units as u
	from astropy.io import ascii
	import mskpy.observing as obs
	import astropy.coordinates as coord
	from astropy import units as u
	from astropy.coordinates import SkyCoord

	#------------------------------------------------------------#
	#	INPUT SAMPLE
	#------------------------------------------------------------#
	'''
	observatory		= 'SAO'
	save_path		= './'
	obspath			= "/home/gw/Research/observatory.txt"
	catpath			= 'MS181101ab_Preliminary-all_candidates.txt'
	start			= '2019/04/17'
	end				= '2019/04/19'
	#altitute limit and moon seperation, moon serperation is a little bit close (2~3 deg)
	numlimit		= 100
	altlimit		= 30.
	moonseperation	= 40.
	sunlimit		= '-18'
	'''
	#------------------------------------------------------------#
	#	OBSERVATORY INFO.
	#------------------------------------------------------------#
	obsinfo			= ascii.read(obspath)
	obsname     	= np.copy(obsinfo['name'])
	obsindex		= np.where(obsname == observatory)[0]
	obslat			= (np.copy(obsinfo['latitude(N+)'])[obsindex])[0]
	obslon			= (np.copy(obsinfo['longitude(E+)'])[obsindex])[0]
	obsalt			= (np.copy(obsinfo['altitude'])[obsindex])[0]
	obstz			= (np.copy(obsinfo['timezone'])[obsindex])[0]
	tz				= pytz.timezone(obstz)
	#------------------------------------------------------------#
	observ			= ephem.Observer()
	observ.lat		= str(obslat)
	observ.lon		= str(obslon)
	observ.elevation= obsalt
	observ.horizon	= sunlimit
	#------------------------------------------------------------#
	#objects from catalog file
	tdata			= ascii.read(catpath)
	objname			= tdata['name']
	ra				= tdata['ra']
	dec				= tdata['dec']
	prior			= tdata['sort']
	rank			= tdata['rank']			
	dist			= tdata['dist']

	RA				= coord.Angle(ra, unit = u.deg)
	Dec				= coord.Angle(dec, unit = u.deg)

	radd			= RA.value
	rad				= RA.hour
	decd			= Dec.value
	decdd			= Dec.degree

	#angular distance calculation
	def angsep(ra1deg, dec1deg, ra2deg, dec2deg) : 
		ra1rad		= ra1deg*np.pi/180
		dec1rad		= dec1deg*np.pi/180
		ra2rad		= ra2deg*np.pi/180
		dec2rad		= dec2deg*np.pi/180
		cos_a		= np.sin(dec1rad)*np.sin(dec2rad)+(np.cos(dec1rad)*np.cos(dec2rad)*np.cos(ra1rad-ra2rad))
		anglesep	= np.arccos(cos_a)*180/np.pi
		return anglesep

	#dates to calculate
	fmt				= '%Y/%m/%d'
	startdt			= datetime.datetime.strptime(start, fmt)
	enddt			= datetime.datetime.strptime(end, fmt)
	startmjd		= (jdcal.gcal2jd(startdt.year, startdt.month, startdt.day))[1]
	endmjd			= (jdcal.gcal2jd(enddt.year, enddt.month, enddt.day))[1]

	for i in range(int(endmjd-startmjd+1)):
		onedaymjd   = startmjd+i+1
		oneday      = jdcal.jd2gcal(2400000.5, onedaymjd)
		onedaydt    = datetime.datetime(oneday[0], oneday[1], oneday[2])
		dst         = tz.dst(onedaydt, is_dst=True)
		dst         = dst.seconds/3600

		onedaydt    = datetime.datetime(oneday[0], oneday[1], oneday[2], tzinfo=tz)
		onedayutc   = onedaydt.astimezone(pytz.utc)
		observ.date = onedayutc

		# Moon distance and information
		mcoord		= ephem.Moon()
		mcoord.compute(observ)
		minfo		= 'Moon ra, dec : '+str(mcoord.ra)+' '+str(mcoord.dec)+'\n'
		mphase		= ephem.Moon(observ.date)
		mphasestr	= 'Moon phase   : '+ "%.2f" % mphase.moon_phase +'\n'
		msep		= angsep(radd, decdd, np.degrees(mcoord.ra), np.degrees(mcoord.dec))

		#	SUNSET CALC.
		sunset      = observ.previous_setting(ephem.Sun())
		sunsettu    = ephem.Date.tuple(sunset)
		sunsetdt    = datetime.datetime(sunsettu[0],sunsettu[1],sunsettu[2],sunsettu[3],int(sunsettu[4]),tzinfo=pytz.utc)
		sunsetlocal = sunsetdt.astimezone(tz)
		sunsetstr   = sunlimit+' deg sunset : '+str(sunsetlocal.hour)+':'+str(sunsetlocal.minute)+'\n'
		sunsethour  = sunsetlocal.hour+sunsetlocal.minute/60.+sunsetlocal.second/3600.
		#	SUNRISE CALC.
		sunrise      = observ.next_rising(ephem.Sun())
		sunrisetu    = ephem.Date.tuple(sunrise)
		sunrisedt    = datetime.datetime(sunrisetu[0],sunrisetu[1],sunrisetu[2],sunrisetu[3],int(sunrisetu[4]),tzinfo=pytz.utc)
		sunriselocal = sunrisedt.astimezone(tz)
		sunrisestr   = sunlimit+' deg sunrise : '+str(sunriselocal.hour)+':'+str(sunriselocal.minute)+'\n'
		sunrisehour  = sunriselocal.hour+sunriselocal.minute/60.+sunriselocal.second/3600.
		#print (observatory)
		#print ('Local mid night in UTC : '+str(observ.date))
		#print (minfo,mphasestr,sunsetstr,sunrisestr)

		#	MAKE RESULT FILE
		stryear      = str(oneday[0])
		strmonth     = str(oneday[1])
		strday       = str(oneday[2]-1)
		if int(strmonth) < 10 : strmonth = '0'+strmonth
		if int(strday)   < 10 : strday   = '0'+strday

		f			= open(save_path+'/'+headname+'-'+stryear+strmonth+strday+"-rts_vis-"+observatory+".txt",'w')
		f.write('#\t'+str(observ.date)+' UTC & Day Time Saving +'+str(dst)+'\n')
		f.write('#\tObservatory\t= '+observatory+'\n')
		f.write('#\t'+sunsetstr)
		f.write('#\t'+sunrisestr)
		f.write('#\t'+minfo)
		f.write('#\t'+mphasestr)
		f.write('#\tMoon seperation = '+str(moonseperation)+'\n')
		f.write('#\tAltitude limit = '+str(altlimit)+'\n')
		f.write('#\tRank : the lower rank, the higher priority\n')
		f.write('#------------------------------------------------------- \n')
		f.write('name ra dec rise(LT) transit(LT) set(LT) moon_dist(deg) distance(Mpc) rank\n')

		numcount	= 0
		for n in range(len(rad)):
			#calculate rise transit set time with altitute limit
			param_rts	= dict(	ra=radd[n],
								dec=decdd[n],
								date=onedaydt,
								lon=obslon,
								lat=obslat,
								tz=obstz,
								limit=altlimit,
								precision=1440)
			rtscal		= obs.rts(**param_rts)
			rt			= rtscal[0]
			tt			= rtscal[1]
			st			= rtscal[2]
			if rtscal[0]== None:
				#print (objname[n],ra[n],dec[n], rtscal[0], rtscal[1], rtscal[2],"%.2f" % msep[n])
				pass
			elif sunrisehour < rtscal[0] < sunsethour and sunrisehour < rtscal[2] < sunsethour and sunrisehour < rtscal[1] < sunsethour: 
				#print (objname[n]+' It can be seen in daytime!')
				pass
			elif msep[n] < moonseperation or msep[n] > 360-moonseperation:
				#print (objname[n]+' too close to Moon < '+str(moonseperation)+' deg')
				pass
			else:
				if numcount < numlimit:
					c= SkyCoord(ra=ra[n]*u.degree, dec=dec[n]*u.degree, frame='icrs')
					c_ra= c.ra.hms
					c_dec= c.dec.dms
					
					nra='%02d:%02d:%.3f' %(c_ra[0], abs(c_ra[1]), abs(c_ra[2]))
					ndec='%02d:%02d:%.3f' %(c_dec[0], abs(c_dec[1]), abs(c_dec[2]))
					
					rtp	="%.2d" % int(rt)+':'+"%.2d" % int((rt-int(rt))*60)
					ttp	="%.2d" % int(tt)+':'+"%.2d" % int((tt-int(tt))*60)
					stp	="%.2d" % int(st)+':'+"%.2d" % int((st-int(st))*60)
					vis	='{:8s}	{:12s}	{:12s}	{:5s}	{:5s}	{:5s}	{:3s}		{:4s}	{:4s}'.format(objname[n],str(nra),str(ndec),rtp,ttp,stp,str(int(msep[n])),str(int(dist[n])),str(rank[n]))+'\n'
					f.write(vis)
					#print (objname[n],ra[n],dec[n], rtp,ttp,stp,"%.2f" % msep[n])
					numcount+= 1
				else:
					pass
				'''
				if numcount < numlimit:
					
					rtp="%.2d" % int(rt)+':'+"%.2d" % int((rt-int(rt))*60)
					ttp="%.2d" % int(tt)+':'+"%.2d" % int((tt-int(tt))*60)
					stp="%.2d" % int(st)+':'+"%.2d" % int((st-int(st))*60)
					vis='{:24s}	{:12s}	{:12s}	{:5s}	{:5s}	{:5s}	{:3s}	{:2s}'.format(objname[n],str(ra[n]),str(dec[n]),rtp,ttp,stp,str(int(msep[n])),str(prior[n]))+'\n'
					f.write(vis)
					#print (objname[n],ra[n],dec[n], rtp,ttp,stp,"%.2f" % msep[n])
					numcount+= 1
				else:
					pass
				'''
		f.close()
#-------------------------------------------------------------------------#
def getccdinfo(obs, path_obs):

	'''
	GET CCD INFORMATION (GAIN, PIXEL SCALE, FOV)

	gain, pixscale, fov = getccdinfo(obs, path_obs)
	
	INPUT:
	path_obs = '/home/sonic/Research/table/obs.txt'

	OUTPUT:
	gain, pixscale, fov
	'''
	obstbl = ascii.read(path_obs)
	indx_obs = np.where(obstbl['obs']==obs)

	outdict = dict()
	gain = obstbl[indx_obs]['gain'][0]
	pixscale = obstbl[indx_obs]['pixelscale'][0] *  u.arcsecond / u.pixel
	fov = obstbl[indx_obs]['fov'][0] * u.deg * u.deg
	rdnoise = obstbl[indx_obs]['RDnoise'][0]

	outdict['obs'] = obs
	outdict['gain'] = gain * u.electron / u.second
	outdict['pixelscale'] = pixscale
	outdict['fov'] = fov
	outdict['rdnoise'] = rdnoise

	return outdict
#-------------------------------------------------------------------------#
'''
def wcsremap(inim, tempim, outim='wr.fits'):
	# outim = '{}/wr{}'.format(os.path.dirname(inim), os.path.basename(inim))
	# com = 'wcsremap -template {} -source {} -outIm {}'.format(inim, tempim, outim)
	# com = 'wcsremap -template {} -source {} -outIm {}'.format(tempim, inim, outim)
	com = f'wcsremap -template {tempim} -source {inim} -outIm {outim}'.format(tempim, inim, outim)
	print(f"""INPUT IMAGE\t: {inim}
	TEMP IMAGE\t: {tempim}
	OUTPUT IMAGE\t: {outim}""")
	print(com)
	os.system(com)
'''
#-------------------------------------------------------------------------#
def hotpants(inim, refim, iu=6000000000, il=-100000, tu=6000000000, tl=-100000):

	if os.path.dirname(inim) == '':
		interval = './'
	else:
		interval = '/'

	outim = '{}{}hd{}'.format(os.path.dirname(inim), interval, os.path.basename(inim))
	convfile = '{}{}hc{}'.format(os.path.dirname(inim), interval, os.path.basename(inim))
	com = f'hotpants -c t -n i -iu {iu} -il {il} -tu {tu} -tl {tl} -v 0 -inim {inim} -tmplim {refim} -outim {outim} -oci {convfile}'

	print(com)
	os.system(com)
	return outim
#-------------------------------------------------------------------------#
def epochimcomb(imlist, outim='imcomb.fits', path_save='.'):
	'''
	epochimcomb(imlist, outim='imcomb.fits', path_save='.')

	imlist = glob.glob('Calib*20181229*.fits')
	epochimcomb(imlist)
	'''
	#------------------------------------------------------------
	import numpy as np
	from astropy.nddata import fits_ccddata_reader, fits_ccddata_writer
	# from astropy.nddata import CCDData
	from matplotlib import pyplot as plt  
	from ccdproc import Combiner
	from astropy.time import Time
	from astropy.io import fits
	from imsng import phot
	#------------------------------------------------------------
	#	EXTRACT INFO. FROM THE FIRST IMAGE
	#------------------------------------------------------------
	data0 = fits_ccddata_reader(imlist[0], unit='adu')
	meta0 = data0.meta
	wcs0 = data0.wcs
	part = imlist[0].split('-')
	#------------------------------------------------------------
	#	IMAGE COMBINE
	#------------------------------------------------------------
	comlist = []
	dateobslist = []
	explist = []
	print('{} IMAGE COMBINE START\n'.format(len(imlist)))
	for inim in imlist:
		print(inim)
		hdr = fits.getheader(inim)
		dateobslist.append(Time(hdr['DATE-OBS'], format='isot').jd)
		explist.append(hdr['EXPTIME'])
		comlist.append(fits_ccddata_reader(inim, unit='adu'))

	dateobs = Time(np.mean(dateobslist), format='jd')
	totexp = np.sum(explist)
	try:
		comim = '{}-{}-{}-{}-{}-{}-{}-com.fits'.format(part[0], part[1], part[2], dateobs.isot[0:10].replace('-', ''), dateobs.isot[11:19].replace(':', ''), part[5], int(totexp))
	except:
		print('IMAGE NAME FORMAT IS NOT Calib-... .fits.')
		comim = outim
	c = Combiner(comlist)
	cdata = c.median_combine()
	cdata.meta = meta0
	cdata.wcs = wcs0
	print('OUTPUT IMAGE :\t{}\n'.format(comim))
	fits_ccddata_writer(cdata, path_save+'/'+comim)
	#------------------------------------------------------------
	puthdr(comim, 'TOTEXP', totexp, hdrcomment='Total exposure time in seconds')
	puthdr(comim, 'JD', dateobs.jd, hdrcomment='Center Julian Date at start of exposure')
	puthdr(comim, 'MJD', dateobs.mjd, hdrcomment='Center Modified Julian Date at start of exposure')
	puthdr(comim, 'DATE-OBS', dateobs.isot, hdrcomment='YYYY-MM-DDThh:mm:ss observation start, UT')
	puthdr(comim, 'NCOMBINE', len(imlist), hdrcomment='THE NUMBER OF COMBINED IMAGES')
	for i, inim in enumerate(imlist):
		puthdr(comim, 'COMBINE{}'.format(i+1), inim, hdrcomment='{} COMBINED IMAGE'.format(i+1))
	print('DONE')
	return comim
#------------------------------------------------------------
def combname(imlist):
	import numpy as np
	from astropy.time import Time
	from astropy.io import fits
	#------------------------------------------------------------
	#	EXTRACT INFO. FROM THE FIRST IMAGE
	#------------------------------------------------------------
	part = imlist[0].split('-')
	#------------------------------------------------------------
	#	IMAGE COMBINE
	#------------------------------------------------------------
	comlist = []
	dateobslist = []
	explist = []

	for inim in imlist:
		# print(inim)
		hdr = fits.getheader(inim)
		dateobslist.append(Time(hdr['DATE-OBS'], format='isot').jd)
		explist.append(hdr['EXPTIME'])

	dateobs = Time(np.mean(dateobslist), format='jd')
	totexp = np.sum(explist)

	comim = '{}-{}-{}-{}-{}-{}-{}.com.fits'.format(part[0], part[1], part[2], dateobs.isot[0:10].replace('-', ''), dateobs.isot[11:19].replace(':', ''), part[5], int(totexp))
	return comim, hdr, dateobs, totexp
#------------------------------------------------------------
def swarpcomb(imlist, gain, pixscale, racent=None, decent=None, ximage=None, yimage=None, listname='obj.list', path_save='.', keys_to_get=['OBJECT', 'FILTER']):
	import os, glob
	# import numpy as np
	# from imsng import tool, phot
	'''
	imlist = glob.glob('Calib*.fits')
	path_save = '.'
	path_obs = '/home/sonic/Research/table'
	listname = 'obj.list'
	'''

	# imlist = glob.glob(imkey); imlist.sort()
	f = open(listname, 'w')
	for inim in imlist:
		f.write(inim+'\n')
		# print(inim)
	f.close()

	# comim, hdr, dateobs, totexp = combname(imlist)
	comim, hdr, dateobs, totexp = combname(imlist)
	# part = comim.split('-')

	# gain, pixscale, fov = getccdinfo(part[1], path_obs)
	# gain, pixscale, fov = getccdinfo(part[1], path_obs)

	conf = 'default.swarp'
	os.system('swarp -d > {}/{}'.format(path_save, conf))

	if (racent!=None) & (decent!=None) & (ximage!=None) & (yimage!=None):
		com = f'swarp @{listname} -c {conf} -IMAGEOUT_NAME {comim} -COMBINE_TYPE MEDIAN -RESAMPLE Y -PIXEL_SCALE {pixscale} -GAIN_DEFAULT {gain} -SUBTRACT_BACK Y -CENTER {racent},{decent} -IMAGE_SIZE {ximage},{yimage} -SUBTRACT_BACK N'
	else:
		com = f'swarp @{listname} -c {conf} -IMAGEOUT_NAME {comim} -COMBINE_TYPE MEDIAN -RESAMPLE Y -PIXEL_SCALE {pixscale} -GAIN_DEFAULT {gain} -SUBTRACT_BACK Y'

	print(com)
	os.system(com)

	for key in keys_to_get:
		puthdr(comim, key, hdr[key], hdrcomment=key)

	puthdr(comim, 'TOTEXP', totexp, hdrcomment='Total exposure time in seconds')
	puthdr(comim, 'EXPTIME', totexp, hdrcomment='Total exposure time in seconds')
	puthdr(comim, 'JD', dateobs.jd, hdrcomment='Center Julian Date at start of exposure')
	puthdr(comim, 'MJD', dateobs.mjd, hdrcomment='Center Modified Julian Date at start of exposure')
	puthdr(comim, 'DATE-OBS', dateobs.isot, hdrcomment='YYYY-MM-DDThh:mm:ss observation start, UT')
	puthdr(comim, 'NCOMBINE', len(imlist), hdrcomment='THE NUMBER OF COMBINED IMAGES')
	for i, inim in enumerate(imlist):
		puthdr(comim, 'COMBINE{}'.format(i+1), inim, hdrcomment='{} COMBINED IMAGE'.format(i+1))

	os.system('rm coadd.weight.fits default.swarp obj.list swarp.xml')

	return comim
#------------------------------------------------------------

def swarpalign(imlist, gain, pixscale, racent=None, decent=None, ximage=None, yimage=None, listname='obj.list', outim='align.fits'):
	import os, glob
	# import numpy as np
	# from imsng import tool, phot
	'''
	imlist = glob.glob('Calib*.fits')
	path_save = '.'
	path_obs = '/home/sonic/Research/table'
	listname = 'obj.list'
	'''

	# imlist = glob.glob(imkey); imlist.sort()
	f = open(listname, 'w')
	for inim in imlist:
		f.write(inim+'\n')
		# print(inim)
	f.close()

	conf = 'default.swarp'
	os.system(f'swarp -d > ./{conf}')

	com = f'swarp @{listname} -c {conf} -IMAGEOUT_NAME {outim} -COMBINE_TYPE MEDIAN -RESAMPLE Y -PIXEL_SCALE {pixscale} -GAIN_DEFAULT {gain} -SUBTRACT_BACK Y -CENTER {racent},{decent} -IMAGE_SIZE {ximage},{yimage} -SUBTRACT_BACK N'
	print(com)
	os.system(com)

	os.system('rm coadd.weight.fits default.swarp obj.list swarp.xml')
#------------------------------------------------------------
def trim(inim, position, size, outim='trim.fits'):
	# Load the image and the WCS
	hdu = fits.open(inim)[0]
	wcs = WCS(hdu.header)
	# Make the cutout, including the WCS
	cutout = Cutout2D(hdu.data, position=position, size=size, wcs=wcs)
	# Put the cutout image in the FITS HDU
	hdu.data = cutout.data
	# Update the FITS header with the cutout WCS
	hdu.header.update(cutout.wcs.to_header())
	# Write the cutout to a new FITS file
	hdu.writeto(outim, overwrite=True)
#------------------------------------------------------------
def calc_app(mag, magerr, gwdist0, gwdiststd0, gwdist1, gwdiststd1):
	import numpy as np
	app		= mag+5*np.log10(gwdist1/gwdist0)
	apperr	= np.sqrt( (magerr)**2 + ((5*gwdiststd1)/(np.log(5)*gwdist1))**2 + ((5*gwdiststd0)/(np.log(5)*gwdist0))**2 )
	return app, apperr
#------------------------------------------------------------
def abs2app(mag, magerr, gwdist, gwdiststd):
	import numpy as np
	app		= 5*np.log10(gwdist)-5+mag
	apperr	= 5*gwdiststd/(gwdist*np.log(10))
	return app, apperr
#------------------------------------------------------------
def z2dist(z):
	from astropy import units as u
	from astropy import constants as const
	import numpy as np
	H0 = 70 * u.km / (u.second * u.Mpc)
	c = const.c.to(u.km / u.second)
	d = c*z/H0
	return d
#------------------------------------------------------------
def limitmag(ul0, t0, t):
	import numpy as np
	ul  = ul0 -(-2.5*np.log10(np.sqrt(t/t0)))
	return ul
#------------------------------------------------------------
def exptime4limitmag(ul0, ul, t0):
	import numpy as np
	t = t0*(10.**(2*((ul-ul0)/2.5)))
	return t
#------------------------------------------------------------
def ToO_request(ul0, m0, n, nsigma=5):
	'''
	ul0 : base n sigma depth
	m0 : target magnitude
	n : n*exposure time
	nsigma : ? sigma depth (default 5)

	return depth, target magnitude error
	'''
	import numpy as np
	ul = ul0+2.5*np.log10(np.sqrt(n))
	dul = ul-m0
	mer = 1./(nsigma*(dul*(100**0.2)))
	
	return round(ul, 3), round(mer, 3)
#------------------------------------------------------------
def sqsum(a, b):
	'''
	SQUARE SUM
	USEFUL TO CALC. ERROR
	'''
	return np.sqrt(a**2.+b**2.)
	#------------------------------------------------------------
def puthdr(inim, hdrkey, hdrval, hdrcomment=''):
	from astropy.io import fits
	hdr		=	fits.getheader(inim)
	fits.setval(inim, hdrkey, value=hdrval, comment=hdrcomment)	
	comment     = inim+'\t'+'('+hdrkey+'\t'+str(hdrval)+')'
#------------------------------------------------------------
def gregistering(images_to_align, ref_image):
	import os
	# import sys, glob
	import alipy
	# from multiprocessing import Process, Pool
	# import multiprocessing as mp
	import time
	starttime	= time.time()
	if ref_image == '': ref_image = images_to_align[0]
	identifications = alipy.ident.run(ref_image, images_to_align, visu=False)
	for id in identifications: # list of the same length as images_to_align.
		if id.ok == True: # i.e., if it worked
			print("%20s : %20s, flux ratio %.2f" % (id.ukn.name, id.trans, id.medfluxratio))
		else:
			print("%20s : no transformation found !" % (id.ukn.name))
	outputshape = alipy.align.shape(ref_image)
	for id in identifications:
		if id.ok == True:
			params_align	= dict(	filepath	= id.ukn.filepath,
									uknstarlist	= id.uknmatchstars,
									refstarlist	= id.refmatchstars,
									shape		= outputshape,
									outdir		= os.path.dirname(ref_image),
									makepng		= False)
			alipy.align.irafalign(**params_align)
	deltime = time.time() - starttime
	print('All PROCESS IS DONE.\t('+str(round(deltime, 1))+' sec)')
#-------------------------------------------------------------------------#
def wcsremap(inim, refim, outim, path_com='/data3/wcsremap/wcsremap-1.0.1/wcsremap'):
	import os
	com = f'{path_com} -template {refim} -source {inim} -outim {outim}'
	print(com)
	os.system(com)
	return outim
#-------------------------------------------------------------------------#
def imcombine_routine(images_to_align, ref_image):
	'''
	path_data = '/data3/paek/factory/doao/20201209-1m-IMSNG'
	images_to_align = sorted(glob.glob('/data3/paek/factory/doao/20201209-1m-IMSNG/Calib-DOAO*-R-60.fits'))
	ref_image = '/data3/paek/factory/doao/20201209-1m-IMSNG/Calib-DOAO-NGC6946-20201209-094720-R-60.fits'
	'''
	from pyraf import iraf
	# from imsng import tool
	import glob, os
	from astropy.io import fits
	from astropy.time import Time
	import numpy as np

	# images_to_align.remove(ref_image)
	print('Reference image\t: {}'.format(ref_image))
	print('Input images\t:')
	for inim in images_to_align: print(inim)
	hdr = fits.getheader(ref_image)
	obs = os.path.basename(ref_image).split('-')[1]
	obj = os.path.basename(ref_image).split('-')[2]

	#	Image align
	print('#\tIMAGE REGISTERING')
	# for inim in images_to_align:
		# gregistering(inim, ref_image)
	gregistering(images_to_align, ref_image)
	#	List to combine
	comlist = [ref_image]
	for inim in images_to_align:
		comlist.append(outim_gregistering(inim))

	jdlist = []
	exptimes = []
	print('Summon {}/imcombine.list'.format(os.path.dirname(ref_image)))
	f = open('{}/imcombine.list'.format(os.path.dirname(ref_image)), 'w')
	for i, comin in enumerate(comlist):
		f.write('{}\n'.format(comin))
		hdr_tmp = fits.getheader(comin)
		jdlist.append(float(hdr_tmp['jd']))
		exptimes.append(float(hdr_tmp['exptime']))
		hdr['IMCOMB{}'.format(i)] = comin
	f.close()

	exptime = np.sum(exptimes)
	jd = Time(np.mean(jdlist), format='jd')
	dateobs = jd.isot
	utdate = dateobs.split('T')[0].replace('-', '')
	uttime = dateobs.split('T')[1].replace(':', '')[:6]

	outim = '{}/Calib-{}-{}-{}-{}-{}-{}.com.fits'.format(os.path.dirname(ref_image), obs, obj, utdate, uttime, hdr['filter'], int(exptime))

	param_imcomb = dict(
						input="@{}/imcombine.list".format(os.path.dirname(ref_image)),
						output=outim,
						combine="median",
						project="no",
						reject="none",
						scale="none",
						zero="mode",
						)

	print('#\t{} IMAGE IMCOMBINE'.format(len(comlist)))
	iraf.imcombine(**param_imcomb)
	# puthdr(outim, 'DATE-OBS', dateobs, hdrcomment='YYYY-MM-DDThh:mm:ss observation start, UT')
	# puthdr(outim, 'JD', jd.value, hdrcomment='Julian Date at start of exposure')
	# puthdr(outim, 'EXPTIME', exptime, hdrcomment='Exposure time in seconds')
	for i, comin in enumerate(comlist): puthdr(outim, 'IMCOMB{}'.format(i), os.path.basename(comin), hdrcomment='Combined image {}'.format(i))
	return outim
#-------------------------------------------------------------------------#
def outim_gregistering(inim):
	part = os.path.splitext(inim)
	outim = '{}_gregister{}'.format(part[0], part[1])
	return outim
#-------------------------------------------------------------------------#
def subtraction_routine(inim, refim):
	'''
	obs = 'LOAO'
	path_refim = '/data3/paek/factory/ref_frames/{}'.format(obs)
	inim = '/data3/paek/factory/test/Calib-LOAO-NGC6946-20201213-014607-R-180-com.fits'

	obj = 'NGC6946'
	filte = 'R'
	'''
	# inseeing = fits.getheader(inim)['seeing']
	# refseeing = fits.getheader(refim)['seeing']

	# if inseeing > refseeing:
	# 	images_to_align = [inim]
	# 	ref_image = refim
	# else:
	# 	images_to_align = [refim]
	# 	ref_image = inim
	gregistering([refim], inim)
	#	Registered reference image
	grefim = '{}/{}'.format(os.path.dirname(inim), os.path.basename(outim_gregistering(refim)))
	subim = hotpants(inim, grefim, iu=60000, tu=6000000000, tl=-100000)
	ds9com = 'ds9 {} {} {}&'.format(inim, grefim, subim)
	# os.system(ds9com)
	return subim, ds9com
#-------------------------------------------------------------------------#
def subtraction_routine2(inim, refim):
	'''
	obs = 'LOAO'
	path_refim = '/data3/paek/factory/ref_frames/{}'.format(obs)
	inim = '/data3/paek/factory/test/Calib-LOAO-NGC6946-20201213-014607-R-180-com.fits'

	obj = 'NGC6946'
	filte = 'R'
	'''
	# inseeing = fits.getheader(inim)['seeing']
	# refseeing = fits.getheader(refim)['seeing']

	# if inseeing > refseeing:
	# 	images_to_align = [inim]
	# 	ref_image = refim
	# else:
	# 	images_to_align = [refim]
	# 	ref_image = inim
	outim = refim.replace('.fits', '.wcsremap.fits')
	wcsremap(refim, inim, outim)
	#	Registered reference image
	subim = hotpants(inim, outim, iu=60000, tu=6000000000, tl=-100000)
	ds9com = 'ds9 {} {} {}&'.format(inim, outim, subim)
	# os.system(ds9com)
	return outim, ds9com
#-------------------------------------------------------------------------#
def subtraction_routine3(inim, refim):
	outim = refim.replace('.fits', '.astroalign.fits')
	align_img(target_img=inim, reference_img=refim, outputname=outim,)
	#	Registered reference image
	subim = hotpants(inim, outim, iu=60000, tu=6000000000, tl=-100000)
	ds9com = 'ds9 {} {} {}&'.format(inim, outim, subim)
	# os.system(ds9com)
	return outim, ds9com
#-------------------------------------------------------------------------#
def stampimage(inim, x, y, name='IMSNG_transient', outname='./IMSNG_transient.png'):
	import matplotlib.pyplot as plt
	from astropy.io import fits
	from matplotlib.colors import LogNorm
	from matplotlib.patches import Circle
	from astropy.visualization import ZScaleInterval, LinearStretch
	from astropy.wcs import WCS
	from astropy.visualization import (MinMaxInterval, SqrtStretch, ImageNormalize)
	# from ligo.skymap.plot.marker import reticle
	'''
	PLOT IMAGE AND SHOW DESINATED OBJECTS
	rahms = '06:16:39.2560'
	decdms = '-21:29:59.370'
	name = 'J{}{}'.format(rahms[:11].replace(':','').replace('.',''), decdms[:11].replace(':','').replace('.',''))
	imsngname = 'IMSNG {}'.format(name)
	#	input
	# inim = '/data3/IMSNG/IMSNGgalaxies/NGC2207/LOAO/R/Calib-LOAO-NGC2207-20201216-072752-R-600-com.fits'
	inim = '/data1/Test/Calib-LOAO-NGC2207-20201216-072752-R-600-com.fits'
	outname = '{}/test.png'.format(os.path.dirname(inim))
	txt = imsngname
	# 6:16:39.2560, -21:29:59.370
	# 94.16356667, -21.499825
	# 1448.5771, 1631.7396
	'''

	size = 100
	# x, y = round(1448.5771), round(1631.7396)
	x1, x2 = x-size, x+size
	y1, y2 = y-size, y+size

	#	image information
	data0, hdr = fits.getdata(inim, header=True)
	#	data0[y1:y2, x1:x2]
	data = data0[y1:y2, x1:x2]
	# wcs = WCS(hdr)

	#	plot
	plt.close('all')
	plt.rc('font', family='serif')
	fig = plt.figure()
	x = 1080 / 4 / fig.dpi
	y = 1080 / 4 / fig.dpi
	fig.set_figwidth(x)
	fig.set_figheight(y)
	ax = fig.add_subplot(111)

	norm_zscale	= ImageNormalize(data, interval=ZScaleInterval(), stretch=LinearStretch())
	im = ax.imshow(data, cmap='gray', origin='lower', norm=norm_zscale) 

	#	Marker and text --> not use
	# marker = reticle(
	# 					# inner=50.0, outer=150.0,
	# 					which='lt', angle=180
	# 				)
	# ax.plot(round(x), round(y), markersize=100, markeredgewidth=4, marker=marker, color='yellow')
	# ax.text(x, y+size/10, str(txt), color='gold', fontsize=5)

	ax.set_title(name, fontsize=10)
	ax.grid('both', linestyle='--', color='white', alpha=0.5)
	ax.set_xlabel('x [pix]', fontsize=12)
	ax.set_ylabel('y [pix]', fontsize=12)
	plt.tight_layout()
	plt.minorticks_on()

	fig.savefig(outname, bbox_inches='tight', overwrite=True)
#-------------------------------------------------------------------------#
def dict2table(dictionary, path_save):
	print('Dictionary to table @{}'.format(path_save))
	keys = list(dictionary.keys())
	f = open(path_save, 'w')
	f.write('{}\t{}\n'.format('key', 'value'))
	for key in keys:
		f.write('{}\t{}\n'.format(key, dictionary[key]))
	f.close()
#-------------------------------------------------------------------------#
def imsng_name_correction(inim, alltbl, radius):
	'''
	path_alltarget = '/home/paek/table/alltarget.dat'
	alltbl = ascii.read(path_alltarget)
	inim = '/data3/IMSNG/IMSNGgalaxies/NGC0772/CBNUO/R/Calib-CBNUO-NGC0772-20201210-144734-R-180.fits'
	'''
	import os
	from astropy.coordinates import SkyCoord
	import astropy.io.ascii as ascii
	from astropy import units as u
	from astropy.io import fits
	from astropy.wcs import WCS
	from imsng import calib
	#	center x, y --> ra, dec
	w = WCS(inim)
	hdr = fits.getheader(inim)
	xcent, ycent= hdr['NAXIS1']/2., hdr['NAXIS2']/2.
	ra, dec = w.all_pix2world(xcent, ycent, 1)
	#	matching
	c = SkyCoord(ra, dec, frame='icrs', unit='deg')
	c_all = SkyCoord(alltbl['ra'], alltbl['dec'], unit=(u.hourangle, u.deg))
	indx, sep, _ = c.match_to_catalog_sky(c_all)
	#	object (header, closest matched one)
	obj = hdr['object']
	robj = alltbl['obj'][indx]

	if (obj != robj) & (sep < radius):
		print('Image OBJECT header\t\t:{}'.format(hdr['object']))
		print('Real OBJECT field\t\t:{}'.format(robj))
		print('Separation {} and {}\t:{}'.format(obj, robj, sep.to(u.arcmin)))
		puthdr(inim, 'OBJECT', robj)
		# calib.changehdr(inim, 'OBJECT', robj)
		#	file name change
		# newim = inim.replace(obj, robj)
		# mvcom = 'mv {} {}'.format(inim, newim)
		# print(mvcom)
		# os.system(mvcom)
	else:
		print('Object header is correct.')
		pass
	return robj, sep
#-------------------------------------------------------------------------#
def SE_seeing(inim, obs, path_obs, path_config, seeing_assume, frac=0.68, clean=True):
	# import os
	# import numpy as np
	from util import tool
	# from imsng import phot_tbd
	# from astropy import units as u
	# from astropy.io import fits
	# from astropy.io import ascii
	# import matplotlib.pyplot as plt
	# print('Quick seeing measurement with SE')
	'''
	inim = '/data3/paek/factory/loao/2020_1215/afzobj.NGC2207.20201216.0211.fits'
	path_config = '/home/paek/config'
	obs = 'LOAO'
	path_obs = '/home/paek/table/obs.dat'
	seeing_assume = 3 * u.arcsecond
	frac = 0.68
	'''
	#------------------------------------------------------------
	#	Input
	#------------------------------------------------------------
	hdr = fits.getheader(inim)
	a = hdr['naxis1']/2.
	b = hdr['naxis2']/2.
	#------------------------------------------------------------
	#	CCD information
	obsdict = getccdinfo(obs, path_obs)
	gain = obsdict['gain']
	pixscale = obsdict['pixelscale']
	fov = obsdict['fov']
	# rdnoise = obsdict['readoutnoise']
	#------------------------------------------------------------
	#	OUTPUT NAMES
	fmt0 = '.fits'
	fmt1 = '.fit'
	fmt2 = '.fts'

	if fmt0 in inim:
		cat = '{}/{}'.format(os.path.dirname(inim), os.path.basename(inim).replace(fmt0, '.cat'))
	elif fmt1 in inim:
		cat = '{}/{}'.format(os.path.dirname(inim), os.path.basename(inim).replace(fmt1, '.cat'))
	elif fmt2 in inim:
		cat = '{}/{}'.format(os.path.dirname(inim), os.path.basename(inim).replace(fmt2, '.cat'))

	# cat = '{}/{}'.format(os.path.dirname(inim), os.path.basename(inim).replace('.fits', '.cat'))
	#	SE configurations
	param = '{}/simple.param'.format(path_config)
	conv = '{}/simple.conv'.format(path_config)
	nnw = '{}/simple.nnw'.format(path_config)
	conf = '{}/simple.sex'.format(path_config)
	#	SE parameters
	param_insex = dict(
						#------------------------------
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
						#	PHOTOMETRY
						#------------------------------
						GAIN = str(gain.value),
						PIXEL_SCALE = str(pixscale.value),
						#------------------------------
						#	STAR/GALAXY SEPARATION
						#------------------------------
						SEEING_FWHM = str(seeing_assume.value),
						)

	com = gpphot.sexcom(inim, param_insex)
	os.system(com)
	rawtbl = ascii.read(cat)
	#	Point source selection
	indx_sel = np.where(
						(rawtbl['FLAGS'] == 0) &
						(sqsum((rawtbl['X_IMAGE']-a)/a, (rawtbl['Y_IMAGE']-b)/b) < frac) &
						(rawtbl['CLASS_STAR']>0.9) &
						(rawtbl['FWHM_WORLD']>0.0)
						)
	seltbl = rawtbl[indx_sel]
	#	Seeing in arcsecond/pixel as median value
	seeing = np.median(seltbl['FWHM_WORLD'].to(u.arcsecond))
	peeing = np.median(seltbl['FWHM_IMAGE']) * u.pix
	#	Header update
	try:
		puthdr(inim, hdrkey='SEEING', hdrval=seeing.value, hdrcomment='SEEING [arcsec]')
		puthdr(inim, hdrkey='PEEING', hdrval=peeing.value, hdrcomment='PEEING [pix]')
	except:
		print('try/except: Too low stars to measure seeing. Use 3.0 arcsecond seeing.')
		puthdr(inim, hdrkey='SEEING', hdrval=3.0, hdrcomment='SEEING [arcsec]')
		puthdr(inim, hdrkey='PEEING', hdrval=(3.0*u.arcsecond*pixscale).value, hdrcomment='PEEING [pix]')
	#	Clean output catalog
	if clean == True:
		rmcom = 'rm {}'.format(cat)
		# print(rmcom)
		os.system(rmcom)
	else:
		pass
	return seeing, peeing
#------------------------------------------------------------
def cr_removal(inim, gain, rdnoise):
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
					sigclip=4.0,
					sigfrac=0.3,
					objlim=5.0,
					gain=gain, readnoise=rdnoise, 
					pssl=0.0,
					niter=4,
					sepmed=True,
					cleantype='meanmask',
					fsmode='median',
					psfmodel='gauss',
					psffwhm=hdr['seeing'],
					#	Default
					inmask=None,
					satlevel=65536,
					psfsize=7,
					psfk=None, psfbeta=4.765,
					verbose=False
					)

	time_st = time.time()
	_, crdata = detect_cosmics(**param_cr)
	fits.writeto('{}/cr{}'.format(os.path.dirname(inim), os.path.basename(inim)), crdata, hdr, overwrite=True)
	puthdr('{}/cr{}'.format(os.path.dirname(inim), hdrkey='history', hdrval='Cosmic-rays were removed with the LACosmic {}'.format(time.strftime("%c")), hdrcomment=''))
	time_delta = time.time() - time_st
	print('Remove cosmic-ray for {} [{} sec]'.format(inim, round(time_delta, 3)))
#------------------------------------------------------------
def npstr2str(arr):
	outlist = []
	for i in arr:
		outlist.append(str(i))
	outarr = np.array(outlist)
	# return outarr
	return outlist
#------------------------------------------------------------
def obs_summary(filte, ic_cal_phot, ic_com_phot, path_save):
	'''
	#	Observation summary plots
	#	2020.12.26	Created by 	Gregory S.H. Paek
	filte = 'R'
	ic_cal_phot, ic_com_phot
	path_save = path_data
	'''
	#============================================================
	# import os, glob
	import numpy as np
	import matplotlib.pyplot as plt
	# from astropy.table import Table, vstack
	# from astropy.io import ascii
	# from astropy.io import fits
	from astropy.time import Time
	# from astropy.coordinates import SkyCoord
	# from astropy import units as u
	from astropy.wcs import WCS
	# from astropy import constants as const
	# from imsng import phot, tool
	from matplotlib.gridspec import GridSpec
	# import time
	#============================================================
	#	Input
	#============================================================
	# filte = 'R'
	# path_save = path_data
	#============================================================
	#	USER SETTING
	#============================================================
	#	PATH
	#------------------------------------------------------------
	ic_cl_pht = ic_cal_phot.filter(filter=filte)
	ic_cl_pht.summary.sort('jd')
	cltbl = ic_cl_pht.summary[
								(ic_cl_pht.summary['jd'].mask == False) &
								(ic_cl_pht.summary['ul5_1'].mask == False) &
								(ic_cl_pht.summary['seeing'].mask == False) &
								(ic_cl_pht.summary['skyval'].mask == False) &
								(ic_cl_pht.summary['skysig'].mask == False)
								]

	ic_cm_pht = ic_com_phot.filter(filter=filte)
	ic_cm_pht.summary.sort('jd')
	cmtbl = ic_cm_pht.summary[
								(ic_cm_pht.summary['jd'].mask == False) &
								(ic_cm_pht.summary['ul5_1'].mask == False) &
								(ic_cm_pht.summary['seeing'].mask == False) &
								(ic_cm_pht.summary['skyval'].mask == False) &
								(ic_cm_pht.summary['skysig'].mask == False)
								]
	#------------------------------------------------------------
	# t0 = Time('2020-12-05T12:14:04', format='isot', scale='utc')
	time_t0 = Time(np.min(cltbl['jd']), format='jd')
	t0 = time_t0.jd

	cmtbl['delt'] = cmtbl['jd'] - t0
	cltbl['delt'] = cltbl['jd'] - t0

	delt_com = np.copy(cmtbl['delt']).astype('float64')
	delt_cal = np.copy(cltbl['delt']).astype('float64')

	totalt_obs = 24*(np.max(delt_cal) - np.min(delt_cal))
	#============================================================
	#	PLOT
	#------------------------------------------------------------
	plt.close('all')
	plt.rc('font', family='serif')
	fig = plt.figure()
	x = 1920 / 2 / fig.dpi
	y = 1080 / fig.dpi
	fig.set_figwidth(x)
	fig.set_figheight(y)
	#------------------------------------------------------------
	#	GRID
	#------------------------------------------------------------
	ncols = 3
	nrows = 5
	grid = GridSpec(nrows, ncols,
					left=0.1, bottom=0.15, right=0.94, top=0.94, wspace=3, hspace=0.1)
	ax1 = fig.add_subplot(grid[0:1,		0:ncols])
	ax2 = fig.add_subplot(grid[1:2,		0:ncols])
	ax3 = fig.add_subplot(grid[2:3,		0:ncols])
	ax4 = fig.add_subplot(grid[3:4,		0:ncols])
	ax5 = fig.add_subplot(grid[4:nrows,	0:ncols])
	#------------------------------------------------------------
	param_plot = dict(
						# fmt='o-',
						ms=5,
						marker='o',
						mec='k',
						mfc='None',
						color='silver',
						# alpha=0.75,
						alpha=0.5,
					)
	param_plot_com = dict(
							# fmt='o-',
							ms=10,
							marker='v',
							mec='k',
							mfc='None',
							color='silver',
							alpha=0.75,
						)
	#------------------------------------------------------------
	#	Depth
	#------------------------------------------------------------
	depth = np.copy(cltbl['ul5_1']).astype('float64')
	depth_com = np.copy(cmtbl['ul5_1']).astype('float64')
	ax1.plot(delt_cal, depth, **param_plot)#, label=r'5$\sigma$ depth')
	ax1.axhline(y=np.median(depth), color='dodgerblue', alpha=0.5, linestyle='--', label='Single = {}'.format(round(np.median(depth), 1)))

	ax1.plot(delt_com, depth_com, **param_plot_com)#, label=r'5$\sigma$ depth')
	ax1.axhline(y=np.median(depth_com), color='tomato', linestyle='--', label='Combined = {}'.format(round(np.median(depth_com), 1)))
	#------------------------------------------------------------
	#	Seeing
	#------------------------------------------------------------
	seeing = np.copy(cltbl['seeing']).astype('float64')
	# ax2.plot(delt_cal, seeing, 'o-')#, label=r'5$\sigma$ depth')
	ax2.plot(delt_cal, seeing, **param_plot)#, label=r'5$\sigma$ depth')
	ax2.axhline(y=np.median(seeing), color='dodgerblue', alpha=0.5, linestyle='--', label='Median = {}'.format(round(np.median(seeing), 1)))
	#------------------------------------------------------------
	#	skyval
	#------------------------------------------------------------
	skyval = np.copy(cltbl['skyval']).astype('float64')
	# ax3.plot(delt_cal, skyval, 'o-')#, label=r'5$\sigma$ depth')
	ax3.plot(delt_cal, skyval, **param_plot)#, label=r'5$\sigma$ depth')
	ax3.axhline(y=np.median(skyval), color='dodgerblue', alpha=0.5, linestyle='--', label='Median = {}'.format(int(np.median(skyval))))
	#------------------------------------------------------------
	#	skysig
	#------------------------------------------------------------
	skysig = np.copy(cltbl['skysig']).astype('float64')
	skysig_com = np.copy(cmtbl['skysig']).astype('float64')
	ax4.plot(delt_cal, skysig, **param_plot)#, label=r'5$\sigma$ depth')
	ax4.axhline(y=np.median(skysig), color='dodgerblue', alpha=0.5, linestyle='--', label='Single = {}'.format(int(np.median(skysig))))
	ax4.plot(delt_com, skysig_com, **param_plot_com)#, label=r'5$\sigma$ depth')
	ax4.axhline(y=np.median(skysig_com), color='tomato', linestyle='--', label='Combined = {}'.format(int(np.median(skysig_com))))
	#------------------------------------------------------------
	#	Observation #
	#------------------------------------------------------------
	mins = 1
	step = (1/24/60) * mins
	bins = np.arange(0, np.max(delt_cal)+step, step)
	# ax5.hist(delt_cal, bins=bins, color='grey', alpha=0.5, label='bin step = {} [min]'.format(mins))
	ax5.hist(delt_cal, bins=bins, color='grey', alpha=0.5, label='single:{}, com:{}'.format(len(cltbl), len(cmtbl)))
	#============================================================
	#	USER SETTING
	#============================================================
	ax1.legend(loc='best', ncol=2, prop={'size': 12}, framealpha=1.0)
	ax2.legend(loc='best', ncol=2, prop={'size': 12}, framealpha=1.0)
	ax3.legend(loc='best', ncol=2, prop={'size': 12}, framealpha=1.0)
	ax4.legend(loc='best', ncol=2, prop={'size': 12}, framealpha=1.0)
	ax5.legend(loc='best', ncol=2, prop={'size': 12}, framealpha=1.0)
	#------------------------------------------------------------
	ax1.grid(linestyle='--', color='grey', alpha=0.5)
	ax2.grid(linestyle='--', color='grey', alpha=0.5)
	ax3.grid(linestyle='--', color='grey', alpha=0.5)
	ax4.grid(linestyle='--', color='grey', alpha=0.5)
	ax5.grid(linestyle='--', color='grey', alpha=0.5)
	#------------------------------------------------------------
	ax1.minorticks_on()
	ax2.minorticks_on()
	ax3.minorticks_on()
	ax4.minorticks_on()
	ax5.minorticks_on()
	#------------------------------------------------------------
	params_notick = dict(
							axis='x',
							which='both',
							bottom=False,
							top=False,
							labelbottom=False,
							direction='in',
							labelsize=10,
						)
	ax1.tick_params(**params_notick)
	ax2.tick_params(**params_notick)
	ax3.tick_params(**params_notick)
	ax4.tick_params(**params_notick)
	# ax5.tick_params(**params_notick)
	#------------------------------------------------------------
	#	SETTING
	#------------------------------------------------------------
	ax1.set_title('(Total obs. time in {}-band) = {} hours'.format(filte, round(totalt_obs, 2)), fontsize=20)
	left, right = ax1.set_xlim()
	lo, up = ax1.set_ylim()
	ax5.set_xlim([left, right])
	ax1.set_ylim([up, lo])
	ax1.set_ylabel(r'5$\sigma$ depth [AB]', fontsize=12)
	ax2.set_ylabel('Seeing [arcsec]', fontsize=12)
	ax3.set_ylabel('Sky bkg', fontsize=12)
	ax4.set_ylabel('Sky sig', fontsize=12)
	ax5.set_ylabel('Observed #', fontsize=12)
	ax5.set_xlabel('Days after {}'.format(time_t0.isot), fontsize=20)
	#------------------------------------------------------------
	fig.savefig('{}/obs.summary.{}.png'.format(path_save, filte), bbox_inches='tight', dpi=500, overwrite=True)
	fig.savefig('{}/obs.summary.{}.pdf'.format(path_save, filte), bbox_inches='tight', overwrite=True)
#------------------------------------------------------------
def sendmail(path_machine, path_address, path_contents, subject='Hello World!', path_attach_txt='', path_attach_png=''):
	'''
	obs = 'LOAO'
	path_data = '/data3/paek/factory/test'
	path_machine = '/home/paek/table/mail.machine.dat'
	path_address = '/home/paek/table/mail.recivers.txt'
	subject = '[IMSNG] {} '.format(obs)
	path_contents = '{}/helloworld.txt'.format(path_data)
	path_attach_png = None
	path_attach_txt = None
	'''
	#============================================================
	#   Mail
	#	20.12.28	Created by Gregory S.H. Paek
	#============================================================
	import os
	import glob
	import time
	import codecs
	import smtplib
	from astropy.io import ascii
	from email.mime.base import MIMEBase
	from email.mime.text import MIMEText
	from email.mime.image import MIMEImage
	from email.mime.multipart import MIMEMultipart
	from email.header import Header 
	#------------------------------------------------------------
	sttime = time.time()
	comment	= '#\tSend IMSNG pipeline status'
	print(comment)
	#------------------------------------------------------------
	#	Path
	#------------------------------------------------------------
	txtlist = sorted(glob.glob(path_attach_txt))
	pnglist = sorted(glob.glob(path_attach_png))
	#	Table
	reciver = ascii.read(path_address)
	machine = ascii.read(path_machine)
	#	Write mails ...
	# subject = '[IMSNG] {} '.format(obs)
	contents = codecs.open('{}'.format(path_contents), 'rb', 'utf-8')
	#	Mailer
	fromID = machine['id'][0].item()
	fromPW = machine['pw'][0].item()
	#	Recivers
	toIDs = ''
	for address in reciver['address']: toIDs += address+','
	toIDs = toIDs[:-1]
	# toIDs = "gregorypaek94@gmail.com"
	ccIDs = 'gundam_psh@naver.com'
	#------------------------------------------------------------
	#	Mail object
	#------------------------------------------------------------
	msg = MIMEMultipart()
	msg['Subject'] = Header(s=subject, charset="utf-8")
	msg['From'] = fromID
	msg['To'] = toIDs
	if ccIDs != None: msg['Cc'] = ccIDs
	msg.attach(MIMEText(contents.read()))
	#	ATTACH TEXT FILE ON MAIL
	if len(txtlist) != 0:
		for txt in txtlist:
			part = MIMEBase("application", "octet-stream")
			part.set_payload(open(txt, 'rb').read())
			part.add_header('Content-Disposition', 'attachment; filename="%s"'% os.path.basename(txt))
			msg.attach(part)
	else:
		print('No attached text files.')
		pass
	#	ATTACH PNG FILE ON MAIL
	if len(pnglist) != 0:
		for png in pnglist:
			fp = open(png, 'rb')
			img = MIMEImage(fp.read())
			fp.close()
			img.add_header('Content-Disposition', 'attachment', filename=os.path.basename(png))
			msg.attach(img)
	else:
		print('No attached png files.')
	pass
	#	ACCESS TO GMAIL & SEND MAIL
	smtp_gmail = smtplib.SMTP_SSL('smtp.gmail.com', 465)
	smtp_gmail.login(fromID, fromPW)
	smtp_gmail.sendmail(msg["From"], msg["To"].split(",") + msg["Cc"].split(","), msg.as_string())
	smtp_gmail.quit()

	deltime = time.time() - sttime
	comment = 'Send mails to {}. [{} sec]'.format(toIDs, round(deltime))
	print(comment)
#------------------------------------------------------------
import requests
def slack_bot(token, channel, text):
	response = requests.post("https://slack.com/api/chat.postMessage",
		headers={"Authorization": "Bearer "+token},
		data={"channel": channel,"text": text}
	)
	print(response)
#------------------------------------------------------------


#%%
from astropy.io import fits
from astropy import wcs
from ccdproc import CCDData
import os
import astroalign as aa
#%%
def align_img(target_img, reference_img, outputname,):
    """
	Author: HH Choi
	Modified by Gregory S.H. Paek

    parameters
    ----------
    1. target_img : str
            absolute path of a target image 
    2. reference_img : str
            absolute path of a reference image
            The image as a reference frame. All images will be aligned on the basis of this image
    3. prefix = str
            prefix of the aligned iamge
      
    returns
    ----------
    1. outputname : str
            absolute path of the aligned image
    notes
    ----------
    usages
	if __name__ == '__main__':
		tgtimg = '/Users/hhchoi1022/Desktop/Calib-KCT_STX16803-NGC1566-20220325-013537-g-120.fits'
		refimg = '/Users/hhchoi1022/Desktop/Ref-KCT_STX16803-NGC1566-g-4440.com.fits'
		align_img(target_img= tgtimg, reference_img= refimg)
    ----------
    """
    tgt_data, tgt_hdr = fits.getdata(target_img, header = True)
    ref_data, ref_hdr = fits.getdata(reference_img, header = True)
    ref_wcs = wcs.WCS(ref_hdr)
    wcs_hdr = ref_wcs.to_header()
    tgt_hdr.update(wcs_hdr)
    tgt_data = tgt_data.byteswap().newbyteorder()
    ref_data = ref_data.byteswap().newbyteorder()
    
    aligned_data, footprint = aa.register(tgt_data, ref_data, fill_value = 0, detection_sigma = 3, max_control_points = 30)
    aligned_tgt = CCDData(aligned_data, unit = 'adu', header = tgt_hdr)
    # outputname = f'{os.path.dirname(target_img)}/{prefix}{os.path.basename(target_img)}'
    fits.writeto(outputname, aligned_tgt.data, aligned_tgt.header, overwrite = True)
    return outputname


#------------------------------------------------------------

from astropy.wcs import WCS

def get_wcs_coordinates(filename):
    # fits 파일 읽기
    hdulist = fits.open(filename)

    # WCS 정보 가져오기
    wcs = WCS(hdulist[0].header)

    # 이미지 데이터 가져오기
    data = hdulist[0].data

    # 중심 좌표 계산
    center_x = data.shape[1] // 2
    center_y = data.shape[0] // 2
    center_ra, center_dec = wcs.all_pix2world(center_x, center_y, 0)

    # 좌상단 꼭짓점 좌표 계산
    top_left_x = 0
    top_left_y = 0
    top_left_ra, top_left_dec = wcs.all_pix2world(top_left_x, top_left_y, 0)

    # 우상단 꼭짓점 좌표 계산
    top_right_x = data.shape[1] - 1
    top_right_y = 0
    top_right_ra, top_right_dec = wcs.all_pix2world(top_right_x, top_right_y, 0)

    # 좌하단 꼭짓점 좌표 계산
    bottom_left_x = 0
    bottom_left_y = data.shape[0] - 1
    bottom_left_ra, bottom_left_dec = wcs.all_pix2world(bottom_left_x, bottom_left_y, 0)

    # 우하단 꼭짓점 좌표 계산
    bottom_right_x = data.shape[1] - 1
    bottom_right_y = data.shape[0] - 1
    bottom_right_ra, bottom_right_dec = wcs.all_pix2world(bottom_right_x, bottom_right_y, 0)

    # 꼭짓점 좌표를 1차원 배열로 변환
    vertices = np.array([(top_left_ra, top_left_dec),
                         (top_right_ra, top_right_dec),
                         (bottom_right_ra, bottom_right_dec),
                         (bottom_left_ra, bottom_left_dec)])

    # fits 파일 닫기
    hdulist.close()

    # 결과 반환
    return (center_ra, center_dec), vertices
