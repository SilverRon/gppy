# gpPy

[![DOI](https://zenodo.org/badge/631772578.svg)](https://zenodo.org/badge/latestdoi/631772578)

__The gpPy is an automatic pipeline to handle the optical/NIR images from different IMSNG/GECKO facilities from data reduction, astrometry, image stacking, calibration, image subtraction and, transient search using multithreads of CPU with `multiprocessing`.__

Either IMSNG or GECKO aims to find transients in the optical/NIR images to catch an early emission from them and follow up them rapidly.

Purpose of the gpPy is to reduce the workload from a daily data process & analysis from a massive number of images by automation from data reduction to transient search. Feature I want to address is that gpPy converts inhomogeneous raw data from different telesocopes and cameras to homogeneouse output to organize them in our internal data base.

We have tested gpPy since 2020 in our internal system. After numerous debugging processes has been performed, it is now working stably, dealing with *infinte* unexpected situations on the data.

## What is IMSNG/GECKO?

### Intensive Monitoring Survey of Nearby Galaxies (IMSNG; [M. Im et al. (2019)](http://koreascience.or.kr/article/JAKO201912262463280.page))
The IMSNG project is monitoring close and star-forming galaxies with a high cadence ($\rm <1\:day$) to catch the early emission of supernovae (SNe).

### Gratitational-wave Electromagnetic-wave Counterpart in Korea Observatory (GECKO; [M. Im et al. (2019)](http://yokohamagrb2019.wikidot.com/proceedings))
The GECKO project is aiming to find kilonovae (KNe), the optical/NIR counterpart of GW, with the network of 1-2m class telescopes in the world.

### Facilties
Both projects share the same facilities. They consist of more than 10 telescopes, and the `gpPy` can handle the data described in below:

|Facility|Location|Description|
|:---:|:---:|:---:|
|SAO        |Korea       |TBD|
|SOAO       |Korea       |TBD|
|DOAO       |Korea       |TBD|
|KHAO       |Korea       |TBD|
|MDFTS      |Korea       |TBD|
|MAAO       |Korea       |TBD|
|DNSM       |Korea       |TBD|
|CBNUO      |Korea       |TBD|
|LOAO       |USA         |TBD|
|McDonald   |USA         |TBD|
|LSGT       |Austrailia  |TBD|
|KCT        |Austrailia  |TBD|
|RASA36     |Chile       |TBD|
|7DT        |Chile       |TBD|
|KMTNet_SSO |Austrailia  |TBD|
|KMTNet_SAAO|South Africa|TBD|
|KMTNet_CTIO|Chile       |TBD|

# Installation
```
$ git clone https://github.com/SilverRon/gppy
```

## Requirements
We recommand the following requirements about the version:
- This pipeline requires the version of Python == 3.11.3

**We highly recommend to use the state-of-the-art `Python >= 3.11` to maximize the computing speed.**

### Python Library
- `numpy == 1.23.5`
- `scipy == 1.10.1`
- `matplotlib == 3.7.1`
- `astropy == 5.1`
- `astroscrappy == 1.1.0`
- `ccdproc == 2.4.0`
- `requests == 2.28.1`: To utilize slack API
- *`alipy`
<!-- - `PyRAF >= X` -->

*`alipy` for Python 3 needs extra actions to be performed. Please follow below steps (a little bit messy though...):

### External Software
- `SourceEXtractor == 2.19.5`: `sex`
- Astrometry.net: `solve-field`
- `HOTPANTS == 5.1.11`: `hotpants`
- `SWarp == 2.38.0` (if necessary): `swarp`


1. Unzip the `alipy.tar.gz` in `./sample_code` directory
2. Install `alipy`
	```
	cd alipy
	# python setup.py install is not working properly
	pip install .
	```
3. Install `asciidata` in `./sample_code/alipy/asciidata` directory
	```
	cd asciidata/
	pip instsall .
	```
4. Edit source code of `alipy` by editing align.py in `{path_to_site-pakages}/alipy/align.py`
	
	a. Change importing libraries
	```
	# line 7
	## Before
	import alipy

	## After
	import astropy.io.fits as pyfits
	```

	b. Change the `writer` variable
	```
	## Before

	# Step 1, we write the geomap input.
	...
	geomap = open(geomapinpath, "wb")
	writer = csv.writer(geomap, delimiter="\t")
	...
	
	## After
	...
	geomap = open(geomapinpath, "w")
	writer = csv.writer(geomap, delimiter="\t", lineterminator="\n", quotechar="'", quoting=csv.QUOTE_NONE)
	...
	```

<!-- # Structure and Usage
`gpPy` consists of four parts.  -->

# Structure 
1. `gpwatch.py`
	- Executing this script alone automates the following processes:
	
	1. Monitor the uploaded data in certain directories

	2. Check the change of file size - *slack alert!*

	3. If the change is stopped in bytes unit, run the main script 

2. `IMSNG_routine.py` (or `GECKO_routine.py`)
	- This script can be executed automatically via gpwatch.py or used independently to process a single date's data folder for a specific observatory. Processes marked with an asterisk (*) leverage multi-processing capabilities based on the number of threads specified at runtime.

	1. Header Correction
		- Because `gpPy` handles 10+ telescope data with various configurations, this step is the most important to make homogeneous outputs
	2. Making Master frames
		- All master frames are saved in the directory where the user set. If there is no master frames in the same day, it will use the one taked in the nearest date.
		- `bn_median`, which is twcie faster than the `np.median`. 
		1. Master Bias (Zero)
		2. Master Dark
		3. Master Flat
	5. Data Reduction - Bias & Dark Correction and Flat Fielding
		- If you set to correct the fringe pattern in the image (e.g. LOAO z and y-bands), the defringe process will be done. This process requires the master defringe frame and the scale factors of dark and bright areas of the master defringe frame.

	6. *Astrometry (`solve-field` commands from `Astrometry.net`)
	7. *Removal of Cosmic-ray (`LAcosmic`): Computational cost is high.
	8. Convert filename to convertional format
		- Format: `Calib-{obs}-{object or field}-{date}-{time}-{filter}-{exposure time}.fits`
			- The `gpPy` refers and modifies the name format of iTelescope
			- example: Calib-LOAO-M101-20240223-130150-R-60.fits
			- `obs`: Observatory or Telescope name
			- `object or field`: Observed object or field, corresponding the header keyword `OBJECT`
			- `date` & `time`: Observed UTC date and time, corresponding the header keyword `DATE-OBS`
			- `exposure time`: Exposure time, corresponding the header keyword `EXPTIME`
	9. Run photometry code for reduced single frames
		- If all process goes well, there are the same number of `Calib*fits` with the raw OBJECT frames. Photometry will be applied on the `Calib*fits`.

3. `phot/gregoryphot_mp_2021.py`: Photometry for single frames
	- This photometry code requires `gphot.config` file to set the important paramter. If there is no configuration file in the directory where the code runs, it automatically default one (`config/gphot.config`). 
	- Source detection and Photometry (`source-extractor`)
		- To make more accurate photometric measurement and shape definition of sources, `source-extractor` runs twice. 
			1. First run uses `config/growthcurve.*` configurations, and the detection threshold is higher and generate less number of columns than second run. The goal of the first run is to draw a signal-to-noise curve to define the most optimized aperture size where the SNR is the highest, and to measure preliminary FWHM (seeing) of detected sources for input parameter of second run.
			2. Second run uses `config/gregoryphot.*` configurations.  
	- Below is the reference catalog for the zeropoint calibration. Considering the survey coverage and its accuracy, PanSTARR DR1 is highly recommended for the calibration of telescopes in the Northen hemisphere and APASS is recommended for the calibraiton of telescopes in the Southern hemisphere. 
		- PanSTARRS DR1 (_grizy_)
		- SDSS DR12 (_ugriz_)
		- APASS (_BVgri_)
		- 2MASS (_JHK_)
	- The `gpPy` calibrate all magnitudes to the AB magnitude. If the calibration of Johnson-Cousine filter system (_UBVRI_) is required, the conversion equations to use the SDSS or PS1 magnitude are applied (Blaton+07). 

4. Image align & combine
	- Alignment: `gregistering`
	- Combine: `imcombine` task in `pyraf` or `SWarp`
5. Photometry for combined frame
6. Image Subtraction: `HOTPANTS` (Becker+)
7. Photometry for subtracted frame: `phot/gregoryphot_sub_2021.py`
8. Transient Search: `phot/gregoryfind_bulk_mp_2021.py`
	1. Filter transient candidates among detections in the subtracted frame
	2. Generate stamp images for the filtered transient candidates
9. Move output files to designated directory
10. Make the final log file - *slack alert!*

*- Note that this pipeline is designed to handle only IMSNG/GECKO facilities, so some functionalities are limited to our system or structure. If you want to apply it to your system, Feel free to let me know, and we can discuss together.*

<!-- ## Calibration
1. Master Bias
2. Master Dark
3. Master Flat -->

## Photometry
The `gpPy` uses the differential photometry with the reference sources in the images

### Reference catalogs and Conversion of filter systems
All magnitudes that `gpPy` calculates are AB magnitude systems. Some reference catalogs (e.g. APASS) has Vega system (_BV_). Not only that, filter system of telescopes (e.g. LOAO, ...) follows Johnson-Cousine (_BVRI_). 

### Aperture Keywords
- `mag_auto`
- `mag_aper`: Best aperture diameter assuming the gaussian profile (SEEING x1.3462)
- `mag_aper_1`: Set the optimized aperture size based on the SNR curve
- `mag_aper_2`: (Aperture Diameter) = (SEEING x2)
- `mag_aper_3`: (Aperture Diameter) = (SEEING x2)
- `mag_aper_4`: Fixed 3" diameter aperture 
- `mag_aper_5`: Fixed 5" diameter aperture

## Output
```
# reduced & calibrated single frame
Calib*.fits
Calib*.cat

# reduced & calibrated combined frame
Calib*com.fits
hcCalib*com.fits 
hdCalib*com.fits
*.phot.cat
```

# Usage
## Monitoring the new data
You can find the `gpwatch.py`. It has a infinite `while` loof to monitor the newly uploaded data in the specific directory where you direct. 
```
cd ~/gppy

# python gpwatch.py [*observatory name] [**number of cpu to use]
# *require either capital letter or lower case
# **require cpu >= 1

# Example for the LOAO observatory data
python gpwatch.py loao 2
```

Or, you can manually type the observatory name to process like this:

```
> python gpwatch.py 
OBSERVATOR LIST :['LOAO', 'DOAO', 'SOAO', 'CBNUO', 'KHAO', 'KCT_STX16803', 'RASA36', 'LSGT']
obs:loao
```

If there is the new data, it starts to measure the size of the folder (e.g. `2024_0229`) to recognize whether the data is complete its upload. If there is no change in the total data size in the new data directory, then it runs the main pipeline to process a single data folder (`IMSNG_routine.py`).
```
> python gpwatch.py 
OBSERVATOR LIST :['LOAO', 'DOAO', 'SOAO', 'CBNUO', 'KHAO', 'KCT_STX16803', 'RASA36', 'LSGT']
obs:loao
============================================================

[`gpwatch`/o_o] Watching new data for LOAO with 1 cores 

============================================================
00:00:00
<Response [200]>
[`gpwatch`/LOAO] Detected New 2024_0229 Data
python ./IMSNG_Routine.py LOAO 1
```

## How to implement the data configuration of new observatory
The primary feature of gpPy is its expandability, allowing it to handle various datasets from different observatories efficiently. Therefore, implementing a new configuration is straightforward and can be done in a matter of minutes. Follow these steps to get started:

1. Basic Information: First, gather the basic information about the camera and telescope. The most crucial parameters to note are as follows:
	```
	obs             ccd         gain    RDnoise     dark    pixelscale  fov
	LOAO            E2V         2.68    4.84        0.0     0.794       28.1
	...
	```
	The location of `obs.dat` in current version is located here: `config/obs.dat`. 

2. Observational Data Structure: Next, organize the observational data using the following folder tree structure:
	```
	{observatory}/{date}
	├── image.000.fits
	├── image.001.fits
	├── image.002.fits
	...
	```
3. Log File: Finally, create a log file with a specific filename format ({observatory}.log) that contains the following information:
	```
	date
	/data6/obsdata/LOAO/2024_0101
	/data6/obsdata/LOAO/2024_0102
	/data6/obsdata/LOAO/2024_0103
	...
	```
<!-- path_obs = '/home/paek/table/obs.dat' -->
<!-- path_changehdr = '/home/paek/table/changehdr.dat' -->

<!-- ## Processing the data
There are two pipelines depending on the types of dat, one is for the IMSNG, the other for the GECKO. The main difference between them is that GECKO pipeline has much simpler steps to make faster process for the real-time identification of transient, but they are fundamentally identical.

### `IMSNG_routine.py`
This pipeli

### `GECKO_routine.py`

## 3. Photometry
They are located in `./phot/`.

### `gregoryphot_mp_2021.py`
TBD

### `gregoryphot_sub_2021.py`
TBD

## 4. Transient Search
### `gregoryfind_bulk_mp_2021.py` -->

# Future Update
<!-- - GPU-accelerated data reduction and `SourceEXtractor` for a large size of images -->
- WCS correction ([Github](https://github.com/evertrol/sippv))
- Handle different binning size
- ...

# License and Copyright
TBD

# Contact
|Name|Institution|Email|ORCID|
|:---:|:---:|:---:|:---:|
|Gregory S.H. Paek|*SNU|gregorypaek94_at_g_mail|[ORCID](https://orcid.org/my-orcid?orcid=0000-0002-6639-6533)|

*SNU: Seoul National Univeristy

<!-- # 9. Reference
1. IMSNG
2. GECKO -->

# Acknowledgments
Thanks to our IMSNG team members, Prof. Myungshin Im, Dr. Changsu Choi, Dr. Seo-won Chang, Dr. Gu Lim, and Dr. Sophia Kim.
Especially, special thanks to Dr. Changsu Choi, who made techincal foundations in the beggining of the IMSNG project, inspired and motivated me to develop this pipeline.

# Version log
- 2024.03.01: v0.2 documentation 
- 2023.09.01: v0.1
- 2018: development start
