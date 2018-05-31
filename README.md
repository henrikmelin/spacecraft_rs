# Spacecraft Remote Sensing 

This python3 project calculates the viewing geometry for cameras and spectrometers mounted on spacecraft orbiting planets in our solar system. This geometry can then be used to project observations onto coordinate grids, such as longitude and latitude. 

## Disclaimer 

This program is free software. It is open source software. It is not space agency software. It is not endorsed by any space agency or any spacecraft mission. It may not work as expected. It may be riddled with bugs. It may produce inconsistent and/or incorrect results. Use at your own risk. 

I'm always interested in hearing thoughts and comments - feel free to contact me. 

## Requirements 

### 1. Python dependencies 

One of many beautiful aspect of open source software is that we can drink from wells that we did not dig. Here, we are using a number of existing modules: 
* [numpy](http://www.numpy.org/)
* [spiceypy](https://github.com/AndrewAnnex/SpiceyPy)


### 2. NAIF pointing kernels

The software will not function without the spacecraft pointing information provided by the [NASA Navigation and Ancillary Information Facility](https://naif.jpl.nasa.gov/naif/) (NAIF) kernels (see the website for a wealth of information and tutorials). These are a collection of files describing the geometry and the pointing of a spacecraft relative to some solar system body for a specific time. Any remote sensing instruments, such as cameras and spectrometers, are defined relative to the spacecraft, with each having their individual kernel. You can get these kernels from the [NAIF FTP site](ftp://naif.jpl.nasa.gov/pub/naif/), for a whole host of spacecraft missions. 

For example, to download (or to update) all the NAIF kernels for the Juno mission to Jupiter,  you can use the following command (but do note that the entire collection is very large, over 64 Gb of data):

```
wget -nH --cut-dirs=4 -m ftp://naif.jpl.nasa.gov/pub/naif/JUNO/
```
Once this behemoth is downloaded, you can use `juno_kernel_loader.py` to automatically load the required set of kernels for a particular time. 

```
import juno_pointing_kernels

# Say that the wget command above downloaded to /path/to/kernels/JUNO/, then:
time = '2017-01-02 12:34:56'
dir  = '/path/to/kernels/JUNO/'
juno_pointing_kernels.load_kernels(dir, timestr = time)
```

Alternatively you can download individual kernels for any specific interval, although that requires some tedious parsing of the `ck/` and `spk/` filenames. 

### 3. Mission data

Since spacecraft missions are funded by public resources, the data they produced are as a rule available free of charge to the public after a certain time period. Both NASA and ESA have archives with vast amounts of data on them, provided in a variety of file-formats. 

## An Example - Juno JIRAM

The Jovian Infrared Auroral Mapper ([JIRAM](http://www.iaps.inaf.it/solarsystem/jiram/about-jiram/)) is the near-infrared imager and spectrograph onboard the Juno spacecraft, currently in orbit about Jupiter. It also has a spectroscopic mode. It has two imaging modes, one at 3.5 microns capturing auroral emissions, and one at 5 microns, capturing thermal emissions from the deep troposphere. In this example, we will use an observation at 5 microns and assign a set of geometries to each of the pixels that make up the image. 

The JUNO JIRAM data is stored on the Planetary Data System (PDS), and data from the Juno perijoves (closest approach) can be downloaded from [here](e.g. ftp://pds-atmospheres.nmsu.edu/PDS/data/jnojir_2003/). 

In this example we will assign geometry to each pixel of the following image from the seventh perijove:

```
https://pds-atmospheres.nmsu.edu/PDS/data/jnojir_2007/DATA/JIR_IMG_RDR_2017191T083944_V01.IMG
```

We will also need the associated label (LBL) file:

```
https://pds-atmospheres.nmsu.edu/PDS/data/jnojir_2007/DATA/JIR_IMG_RDR_2017191T083944_V01.LBL
```

This gorgeous image looks like this:

![Jupiter by JIRAM](https://github.com/henrikmelin/spacecraft_rs/blob/master/resources/JIR_IMG_RDR_2017191T083944_V01.IMG.png)

Our aim here is to add grid lines of longitude and latitude, shown in the example below: 

```
#!/usr/bin/python
import spacecraft_rs as srs
import matplotlib.pyplot as plt

# Read the data and the label file
file    = '/path/to/data/JIR_IMG_RDR_2017191T083944_V01.IMG'
lbl     = srs.read_label(file.replace('.IMG', '.LBL'))
im      = srs.read_image(file, lbl, 'IMAGE')
timestr = str(lbl['START_TIME'])

# Setup the kernel-loader and load the appropriate ones. This is based purely on time.
dir     = '/path/to/kernels/JUNO/'
kl      = srs.juno_kernel_loader(dir, timestr = timestr)

# Generate the geometry
jiram   = srs.juno_jiram_geometry()
jiram.set_mode(lbl['INSTRUMENT_MODE_ID'])
ret     = jiram.sc_pos(timestr)
fov     = jiram.get_corners()

# This method extracts the parameters using nice-names
lat = jiram.extract_param('lat', pixel=1)
lon = jiram.extract_param('lon', pixel=1)

# Plot our result
plt.figure(figsize=(8, 3))
plt.imshow(im, cmap='hot')
plt.contour(lat, levels=range(-90, 90, 10), colors='white', linestyles='solid')
plt.contour(lon, levels=range(10, 360, 10), colors='white', linestyle='solid')
plt.tight_layout()
plt.show()
```

Which produces the result: 

![Jupiter by JIRAM, with a grid!](https://github.com/henrikmelin/spacecraft_rs/blob/master/resources/JIR_IMG_RDR_2017191T083944_V01_grid.png)










