#============================================================
#	python script for gregister using Alipy code
#	iraf-imcombine.py for combine "gregistered" images
#	Usage : 
#	python gregister-script.py str1 ref_image
#	OR JUST RUN AND WRITE INPUT & REF IMAGE
#	2017.12.22	Changsu Choi
#	2019.05.08	GREGORY S.H. PAEK
#============================================================
import os, sys, glob
import alipy
from multiprocessing import Process, Pool
import multiprocessing as mp
import time
#============================================================
def gregistering(images_to_align, ref_image):
	starttime	= time.time()
	if ref_image == '': ref_image = images_to_align[0]
	identifications = alipy.ident.run(ref_image, images_to_align, visu=False)
	for id in identifications: # list of the same length as images_to_align.
		if id.ok == True: # i.e., if it worked
			print "%20s : %20s, flux ratio %.2f" % (id.ukn.name, id.trans, id.medfluxratio)
		else:
			print "%20s : no transformation found !" % (id.ukn.name)
	outputshape = alipy.align.shape(ref_image)
	for id in identifications:
		if id.ok == True:
			params_align	= dict(	filepath	= id.ukn.filepath,
									uknstarlist	= id.uknmatchstars,
									refstarlist	= id.refmatchstars,
									shape		= alipy.align.shape(ref_image),
									outdir		= './',
									makepng		= False)
			alipy.align.irafalign(**params_align)
	deltime		= time.time() - starttime
	print('All PROCESS IS DONE.\t('+str(round(deltime, 1))+' sec)')
#------------------------------------------------------------
#	INPUT
#------------------------------------------------------------
os.system('ls -C *.fits')
try     :
	str1 = sys.argv[1]
except  : str1        = raw_input('Registering Target (i.e. ref.fits) : ')
if str1     == '' : str1 = 'ref.fits'
images_to_align = sorted(glob.glob(str1))
n=len(images_to_align)

try     :
	print(sys.argv[1])
	ref_image	= images_to_align[0]
except  : ref_image        = raw_input('Reference Image (Calib*.fits)      : ')
if ref_image    == '' : ref_image = 'Calib*.fits'

print str(n),'images will be aligned'
for i in images_to_align : print i
#------------------------------------------------------------
#	MULTI PROCESSING
#------------------------------------------------------------
if __name__ == '__main__':
	jobs	= []
	p			= mp.Process(target=gregistering, args=(images_to_align, ref_image))
	jobs.append(p)
	p.start()