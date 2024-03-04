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
def ds9regmaker(filename, intbl, name, racol, deccol):
	import os,sys
	import string
	from astropy.io import ascii 
	import numpy as np
	import math
	'''
	racol		= 'ALPHA_J2000'
	deccol		= 'DELTA_J2000'
	name		= 'NUMBER'
	intbl		= ascii.read(filename)
	'''

	ra		=	intbl[racol]
	dec		=	intbl[deccol]
	starname=	intbl[name]            # intbl['NUMBER']

	radius	= """ 5" """
	color	= "green"
	f		= open(filename+'.reg', 'w')

	head1	= "# Region file format: DS9 version 4.1\n"
	head2	= """global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n"""
	head3	= "fk5\n"

	f.write(head1)
	f.write(head2)
	f.write(head3)

	for n in range(len(ra)):
		body="circle("+str(ra[n])+","+str(dec[n])+","+radius+") # color="+color+" text={"+str(starname[n])+"}\n"	
		f.write(body)
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
	obs = 'SAO'
	path_obs = '/home/sonic/Research/table'

	OUTPUT:
	gain, pixscale, fov
	'''
	obstbl = ascii.read(path_obs+'/obs.txt')
	indx_obs = np.where(obstbl['obs']==obs)
	gain = obstbl[indx_obs]['gain'][0]
	pixscale = obstbl[indx_obs]['pixelscale'][0]
	fov = obstbl[indx_obs]['fov'][0]
	rdnoise = obstbl[indx_obs]['RDnoise'][0]
	return gain, pixscale, rdnoise
#-------------------------------------------------------------------------#
def wcsremap(inim, tempim):
	outim   = 'wr'+tempim
	# com     = 'wcsremap -template '+inim+' -source '+tempim+' -outIm '+outim
	com     = 'wcsremap -template '+tempim+' -source '+inim+' -outIm '+outim	
	print(com)
	os.system(com)
	return outim
#-------------------------------------------------------------------------#
'''
def wcsremap(inim, template):
	templ   = inim
	source  = template
	outim   = 'Ref'+source
	com     = 'wcsremap -template '+templ+' -source '+source+' -outIm '+outim
	print(com)
	os.system(com)
	return outim
'''
#-------------------------------------------------------------------------#
def hotpants(inim, refim, iu=60000, tu=60000, tl=-10000):
	starttime = time.time()
	if os.path.dirname(inim) == '':
		interval = './'
	else:
		interval = '/'
	outim = os.path.dirname(inim)+interval+'hd'+os.path.basename(inim)
	convfile = os.path.dirname(inim)+interval+'hc'+os.path.basename(inim)

	com = 'hotpants -c t -n i -iu {} -tu {} -tl {} -v 0 -inim {} -tmplim {} -outim {} -oci {}'.format(iu, tu, tl, inim, refim, outim, convfile)

	# com = 'hotpants -c t -n i -iu 60000 -tu 60000 -tl -10000 -v 0 -inim ' + \
		# inim+' -tmplim '+refim+' -outim '+outfile+' -oci '+convfile
	os.system(com)
	deltime = time.time() - starttime
	print('All PROCESS IS DONE.\t('+str(round(deltime, 1))+' sec)')
	return outim
#-------------------------------------------------------------------------#
def epochimcomb(imlist, outim='imcomb.fits', path_save='.'):
	'''
	epochimcomb(imlist, outim='imcomb.fits', path_save='.')

	imlist = glob.glob('Calib*20181229*.fits')
	epochimcomb(imlist)
	'''
	#------------------------------------------------------------
	import os, glob
	import numpy as np
	from astropy.nddata import CCDData, fits_ccddata_reader, fits_ccddata_writer
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
	phot.puthdr(comim, 'TOTEXP', totexp, hdrcomment='Total exposure time in seconds')
	phot.puthdr(comim, 'JD', dateobs.jd, hdrcomment='Center Julian Date at start of exposure')
	phot.puthdr(comim, 'MJD', dateobs.mjd, hdrcomment='Center Modified Julian Date at start of exposure')
	phot.puthdr(comim, 'DATE-OBS', dateobs.isot, hdrcomment='YYYY-MM-DDThh:mm:ss observation start, UT')
	phot.puthdr(comim, 'NCOMBINE', len(imlist), hdrcomment='THE NUMBER OF COMBINED IMAGES')
	for i, inim in enumerate(imlist):
		phot.puthdr(comim, 'COMBINE{}'.format(i+1), inim, hdrcomment='{} COMBINED IMAGE'.format(i+1))
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

	comim = '{}-{}-{}-{}-{}-{}-{}-com.fits'.format(part[0], part[1], part[2], dateobs.isot[0:10].replace('-', ''), dateobs.isot[11:19].replace(':', ''), part[5], int(totexp))
	return comim, hdr, dateobs, totexp
#------------------------------------------------------------
def swarpcomb(imlist, listname='obj.list', path_save='.', path_obs = '/home/sonic/Research/table'):
	import os, glob
	import numpy as np
	from imsng import tool, phot
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

	comim, hdr, dateobs, totexp = tool.combname(imlist)
	part = comim.split('-')

	gain, pixscale, fov = tool.getccdinfo(part[1], path_obs)

	conf = 'default.swarp'
	os.system('swarp -d > {}/{}'.format(path_save, conf))

	com = 'swarp @{} -c {} -IMAGEOUT_NAME {} -COMBINE_TYPE MEDIAN -RESAMPLE N -PIXEL_SCALE {} -GAIN_DEFAULT {} -SUBTRACT_BACK Y'.format(listname, conf, comim, pixscale, gain)

	print(com)
	os.system(com)

	phot.puthdr(comim, 'OBJECT', hdr['OBJECT'], hdrcomment='OBJECT')
	phot.puthdr(comim, 'TOTEXP', totexp, hdrcomment='Total exposure time in seconds')
	phot.puthdr(comim, 'JD', dateobs.jd, hdrcomment='Center Julian Date at start of exposure')
	phot.puthdr(comim, 'MJD', dateobs.mjd, hdrcomment='Center Modified Julian Date at start of exposure')
	phot.puthdr(comim, 'DATE-OBS', dateobs.isot, hdrcomment='YYYY-MM-DDThh:mm:ss observation start, UT')
	phot.puthdr(comim, 'NCOMBINE', len(imlist), hdrcomment='THE NUMBER OF COMBINED IMAGES')
	for i, inim in enumerate(imlist):
		phot.puthdr(comim, 'COMBINE{}'.format(i+1), inim, hdrcomment='{} COMBINED IMAGE'.format(i+1))

	os.system('rm coadd.weight.fits default.swarp obj.list swarp.xml')

	return comim
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