#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#============================================================
#	Module
#------------------------------------------------------------
import os
import sys
import glob
import time
#	IMSNGpy modules
sys.path.append('/home/paek/imsngpy')
from misc import *
#	Astropy
from astropy.io import ascii
from astropy.table import Table, vstack
#============================================================
#	Function
#------------------------------------------------------------
def get_size(start_path = '.'):
	total_size = 0
	for dirpath, dirnames, filenames in os.walk(start_path):
		for f in filenames:
			fp = os.path.join(dirpath, f)
			#	skip if it is symbolic link
			if not os.path.islink(fp):
				total_size += os.path.getsize(fp)
	return total_size
#------------------------------------------------------------
#	Path
#------------------------------------------------------------
path_obsdata = '/data6/obsdata'
path_table = '/home/paek/table'
path_log = '/home/paek/log'
path_gppy = '/home/paek/imsngpy/imsngpy'
path_check_rasa36 = '/home/paek/qsopy/monitor/classify_rasa36.py'
path_check_lsgt = '/home/paek/qsopy/monitor/classify_lsgt.py'
path_preprocess = '/home/paek/qsopy/reduction/gregorycalib_2021.py'
#	Slack
keytbl = ascii.read(f'{path_table}/keys.dat')
OAuth_Token = keytbl['key'][keytbl['name']=='slack'].item()
#------------------------------------------------------------
#	Data information
#------------------------------------------------------------
obsdict = dict(
	#	LOAO
	loao=dict(
		path_base='/data6/obsdata/LOAO',
		path_new='',
		log=f'{path_log}/loao.log',
		size=0,	#	[bytes]
		core=1,	#	4
	),
	#	DOAO
	doao=dict(
		path_base='/data6/obsdata/DOAO',
		path_new='',
		log=f'{path_log}/doao.log',
		size=0,	#	[bytes]
		core=1,	#	4
	),	
	#	SOAO
	soao=dict(
		path_base='/data6/obsdata/SOAO',
		path_new='',
		log=f'{path_log}/soao.log',
		size=0,	#	[bytes]
		core=1,	#	4
	),	
	#	CBNUO
	cbnuo=dict(
		path_base='/data6/obsdata/CBNUO',	#	./2021_0101
		path_new='',
		log=f'{path_log}/cbnuo.log',
		size=0,	#	[bytes]
		core=1,	#	4
	),	
	#	KHAO
	khao=dict(
		path_base='/data6/obsdata/KHAO',	#	./2021_0101
		path_new='',
		log=f'{path_log}/khao.log',
		size=0,	#	[bytes]
		core=2,	#	4
	),	
	#	MDFTS
	mdfts=dict(
		path_base='/data6/obsdata/MDFTS',	#	./2021_0101
		path_new='',
		log=f'{path_log}/mdfts.log',
		size=0,	#	[bytes]
		core=2,	#	4
	),	
	#	KCT_STX16803
	kct_stx16803=dict(
		path_base='/data6/obsdata/KCT_STX16803',
		path_new='',
		log=f'{path_log}/kct_stx16803.log',
		size=0,	#	[bytes]
		core=1,	#	4
	),	
	#	RASA36
	rasa36=dict(
		path_base='/data6/obsdata/RASA36',
		path_new='',
		log=f'{path_log}/rasa36.log',
		size=0,	#	[bytes]
		core=1,	#	4
	),	
)
#------------------------------------------------------------
obslist = ['LOAO', 'DOAO', 'SOAO', 'CBNUO', 'KHAO', 'KCT_STX16803', 'RASA36', 'LSGT']
print('OBSERVATOR LIST :', end='')
print(obslist)
obs = input('obs:').upper()
# obs = 'LOAO'
delay = 10
ncore = input('# of cores (i.e. 8):')
'''
print(f"Wrong input in variable 'sphere' (sphere={sphere})")
print('Process all obs. data')
obslist = ['loao', 'doao', 'soao', 'cbnuo',]+['kct_stx16803', 'rasa36']
'''
#============================================================
#	Main body
#------------------------------------------------------------
print(f"{'='*60}\n\n[gpwatch/o_o] Watching new data for {obs} with {ncore} cores \n\n{'='*60}")
st = time.time()
while True:
	try:
		#	Time
		et = time.time()
		delt = int(et - st)
		h = delt // (60*60)
		m = delt // 60 
		s = delt % 60
		timer = '{:02d}:{:02d}:{:02d}'.format(h, m, s)
		print(timer, end="\r")
		log = f"{path_log}/{obs.lower()}.log"
		path_base = f"{path_obsdata}/{obs}"
		#
		logtbl = ascii.read(log)
		dirlist = os.listdir(path_base)
		#	
		for f in dirlist:
			path_new = f"{path_base}/{f}"
			if (path_new not in logtbl['date']) & (f"{path_new}/" not in logtbl['date']) & (os.path.isdir(path_new)):
				print()
				#------------------------------------------------------------
				#	Slack message
				#------------------------------------------------------------	
				channel = '#pipeline'
				text = f'[gpwatch/{obs}] Detected New {os.path.basename(path_new)} Data'
				param_slack = dict(
					token = OAuth_Token,
					channel = channel,
					text = text,
				)
				slack_bot(**param_slack)
				#	
				print(text)
				init_size = get_size(path_new)
				while True:
					if obs.upper() == 'LOAO':
						time.sleep(int(delay*6))
					else:
						time.sleep(int(delay*2))
					now_size = get_size(path_new)
					if init_size != now_size:
						print(f'Still uploading {os.path.basename(path_new)} : {init_size} --> {now_size}')
						init_size = now_size
					else:
						#------------------------------------------------------------
						#	RASA36 exception
						#------------------------------------------------------------
						# if (obs.upper() == 'RASA36'):
						if (obs.upper() == 'SKIP'):
							com = f'python {path_check_rasa36} {path_new}'
							print(com)
							os.system(com)

							if len(dirlist) == len(os.listdir(path_base)):
								com = f"python {path_preprocess} {obs} {ncore}"
								print(com)
								os.system(com)
							else:
								break
						#------------------------------------------------------------
						#	LSGT exception
						#------------------------------------------------------------
						elif (obs.upper() == 'LSGT'):
							com = f'python {path_check_lsgt} {path_new}'
							print(com)
							os.system(com)

							if (len(dirlist) == len(os.listdir(path_base))) & (len(glob.glob(f"{path_new}/*.fit"))>0):
								com = f"python {path_preprocess} {obs} {ncore}"
								print(com)
								os.system(com)
								st = time.time()
							else:
								#	Write log --> skip no data folder
								logtbl = ascii.read(log)
								logtbl_ = Table()
								logtbl_['date'] = [path_new]
								logtbl = vstack([logtbl, logtbl_])
								logtbl.write(f'{path_log}/{obs.lower()}.log', format='ascii.tab', overwrite=True)
								break
				
						else:
							#	Run python code
							com = f"python {path_preprocess} {obs} {ncore}"
							print(com)
							os.system(com)

							print(f"[gpwatch/{obs}] Process for {os.path.basename(path_new)} is done.")
							print(f"{'='*60}\n\n[gpwatch/o_o] Watching new data for {obs} with {ncore} cores \n\n{'='*60}")
							break
	except Exception as e:
		print(e)
		#------------------------------------------------------------
		#	Slack message
		#------------------------------------------------------------	
		channel = '#pipeline'
		text = f'[gpwatch/{obs}] Error\n{e}'
		param_slack = dict(
			token = OAuth_Token,
			channel = channel,
			text = text,
		)
		slack_bot(**param_slack)
		print(text)
	time.sleep(1)

