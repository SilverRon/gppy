# gppy
The gpPy is an automatic pipeline to handle the optical/NIR images from different IMSNG/GECKO facilities from data reduction, astrometry, image stacking, calibration, image subtraction and, transient search.

# Author
Gregory S.H. Paek (백승학)

## Version log
- 2023.09.01: version 1.0
- 2018: development start

# INDEX
- What is IMSNG/GECKO?
- b
- c

# 1. What is IMSNG/GECKO?
## 1.1. Intensive Monitoring Survey of Nearby Galaxies (IMSNG; [M. Im et al. (2019)](http://koreascience.or.kr/article/JAKO201912262463280.page))
The IMSNG project is monitoring close and star-forming galaxies with a high cadence ($\rm <1\:day$) to catch the early emission of SNe.

## 1.2. Gratitational-wave Electromagnetic-wave Counterpart in Korea Observatory (GECKO; [M. Im et al. (2019)](http://yokohamagrb2019.wikidot.com/proceedings))
The GECKO project is aiming to find KN, the optical/NIR counterpart of GW, with the network of 1-2m class telescopes in the world.

## 1.3. Facilties
Both projects share the same facilities. They consist of more than 10 telescopes described in below:
- Korea
	- SAO 1m
	- SOAO 0.61m (~2020)
	- DOAO 1m
	- KHAO 0.4m
	- MDFTS 0.76m
	- MAAO 0.7m
	- DNSM 1.0m
- United States
	- LOAO 1.0m
	- McDonald Observatory 0.25m, 0.8m, 2.1m (~2020)
- Australia
	- LSGT 0.43m
- Chile
	- KCT 0.36m
	- RASA36 0.36m
	- 7-Dimensional Telescope (7DT) - TBD

# 2. Requirements
- This pipeline requires the version of Python > 3.X.
- Recommanded version of Python == 3.11

# 3. Installation
```
$ git clone https://github.com/SilverRon/gppy
```

# 4. Structure and Usage
`gpPy` consists of two parts, one for the start to monitor the upload

`gpwatch.py`: 
```
cd ~/gppy

# python gpwatch.py [*observatory name] [**number of cpu to use]
# *require either capital letter or lower case
# **require cpu >= 1

python gpwatch.py loao 2
```

`IMSNG_routine.py`:

`GECKO_routine.py`:

# 5. Features 
## 5.1. asdf
1. Monitor the uploaded data in certain directories
2. Check the change of file size
3. If the change is stopped in bytes unit, run the main script (`IMSNG_routine.py` or `GECKO_routine.py`)
4. Header Correction
5. Data Reduction
6. Astrometry (`Astrometry.net`)
7. Removal of Cosmic-ray (`LAcosmic`)
8. Convert filename to convertional format
9. Photometry for single frame
	- Reference catalog for the zeropoint calibration
		- PansTARRS DR1
		- SDSS DR?
		- APASS
		- 2MASS
10. Image align & combine
	- Alignment: `gregistering`
	- Combine: `imcombine` task in `pyraf` or `SWarp`
11. Photometry for combined frame
12. Image Subtraction
13. Photometry for subtracted frame

*- Note that this pipeline is designed to handle only IMSNG/GECKO facilities, so some functionalities are limited to our system or structure. If you want to apply it to your system, Feel free to let me know, and we can discuss together.*

## 5.X. Calibration and Photometry
### 5.X.X. keywords
- `mag_auto`
- `mag_aper`
- `mag_aper_1`
- `mag_aper_2`
- `mag_aper_3`

## 5.2. Output
```
Calib*.fits      # reduced & calibrated single frame
Calib*.cat
Calib*com.fits   # reduced & calibrated combined frame
hcCalib*com.fits 
hdCalib*com.fits
*.phot.cat
```

# 6. Future Update
- WCS correction (https://github.com/evertrol/sippv)
- Handle binning
- ...

# 7. License and Copyright
Free!

# 8. Contact
- Gregory S.H. Paek (gregorypaek94___at___gmail.com) @Seoul National University

# 9. Reference
## 9.1. Project
1. IMSNG
2. GECKO
## 9.2. Software
asdf
## 9.3. 
asdf

# 10. Acknowledgments
Thank to our IMSNG team members, Dr. Changsu Choi, Dr. Gu Lim, and Dr. Sophia Kim.
Especially Dr. Changsu Choi, who made techincal foundations in the beggining of the IMSNG project, inspired and motivated me to develop this pipeline.
