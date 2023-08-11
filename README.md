# gpPy
The gpPy is an automatic pipeline to handle the optical/NIR images from different IMSNG/GECKO facilities from data reduction, astrometry, image stacking, calibration, image subtraction and, transient search using multithreads of CPU with `multiprocessing`.

Either IMSNG or GECKO aims to find transients in the optical/NIR images to catch an early emission from them and follow up them rapidly.

Purpose of the gpPy is to reduce the workload from a daily data process & analysis from a massive number of images by automation from data reduction to transient search. Feature I want to address is that gpPy converts inhomogeneous raw data from different telesocopes and cameras to homogeneouse output to organize them in our internal data base.

We have tested gpPy since 2020 in our internal system. After numerous debugging processes has been performed, it is now working stably, dealing with *infinte* unexpected situations on the data.

# Author
Gregory S.H. Paek (백승학)

## Version log
- 2023.09.01: version 0.1
- 2018: development start

# INDEX
- What is IMSNG/GECKO?
- TBD
- ...

# 1. What is IMSNG/GECKO?
## 1.1. Intensive Monitoring Survey of Nearby Galaxies (IMSNG; [M. Im et al. (2019)](http://koreascience.or.kr/article/JAKO201912262463280.page))
The IMSNG project is monitoring close and star-forming galaxies with a high cadence ($\rm <1\:day$) to catch the early emission of supernovae (SNe).

## 1.2. Gratitational-wave Electromagnetic-wave Counterpart in Korea Observatory (GECKO; [M. Im et al. (2019)](http://yokohamagrb2019.wikidot.com/proceedings))
The GECKO project is aiming to find kilonovae (KNe), the optical/NIR counterpart of GW, with the network of 1-2m class telescopes in the world.

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
We recommand the following requirements about the version:
- This pipeline requires the version of Python == 3.11.3

**We highly recommend to use the state-of-the-art `Python >= 3.11` to maximize the computing speed.**

## 2.1. Python Library
- `numpy == 1.23.5`
- `scipy == 1.10.1`
- `matplotlib == 3.7.1`
- `astropy == 5.1`
- `astroscrappy == 1.1.0`
- `ccdproc == 2.4.0`
- `requests == 2.28.1`: To utilize slack API
- *`alipy`
<!-- - `PyRAF >= X` -->

## 2.2. External Software
- `SourceEXtractor == 2.19.5`: `sex`
- Astrometry.net: `solve-field`
- `HOTPANTS == 5.1.11`: `hotpants`
- `SWarp == 2.38.0` (if necessary): `swarp`

*`alipy` for Python 3 needs extra actions to be performed. Please note the section X.X to install the `alipy` for the Python 3.

# 3. Installation
```
$ git clone https://github.com/SilverRon/gppy
```

# 4. Structure and Usage
`gpPy` consists of four parts.

## 4.1. `gpwatch.py`
```
cd ~/gppy

# python gpwatch.py [*observatory name] [**number of cpu to use]
# *require either capital letter or lower case
# **require cpu >= 1

python gpwatch.py loao 2
```
## 4.2. Main
### 4.2.1. `IMSNG_routine.py`
### 4.2.2. `GECKO_routine.py`
### 4.2.3. `7DT_routine.py` (TBD)

## 4.3. Photometry
They are located in `./phot/`.

### 4.3.1. `gregoryphot_mp_2021.py`
TBD

### 4.3.2. `gregoryphot_sub_2021.py`
TBD

## 4.4. Transient Search
### 4.4.1. `gregoryfind_bulk_mp_2021.py`

# 5. Features 
## 5.1. Structure
1. Monitor the uploaded data in certain directories
2. Check the change of file size - *slack alert!*
3. If the change is stopped in bytes unit, run the main script (`IMSNG_routine.py` or `GECKO_routine.py`)
4. Header Correction
5. Data Reduction
	5.1. Bias correction
	5.2. Dark correction
	5.3. Flat Fielding
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
14. Filtering transient candidates among detections in the subtracted frame
15. Generate stamp images for the filtered transient candidates
16. Move output files to designated directory
17. Make the final log file - *slack alert!*

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
- GPU-accelerated data reduction and `SourceEXtractor` for a large size of images
- WCS correction (https://github.com/evertrol/sippv)
- Handle different binning size
- ...

# 7. License and Copyright
TBD

# 8. Contact
- Gregory S.H. Paek (gregorypaek94___at___gmail.com) @Seoul National University

# 9. Reference
1. IMSNG
2. GECKO

# 10. Acknowledgments
Thanks to our IMSNG team members, Prof. Myungshin Im, Dr. Changsu Choi, Dr. Seo-won Chang, Dr. Gu Lim, and Dr. Sophia Kim.
Especially, special thanks to Dr. Changsu Choi, who made techincal foundations in the beggining of the IMSNG project, inspired and motivated me to develop this pipeline.
