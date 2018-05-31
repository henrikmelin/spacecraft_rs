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