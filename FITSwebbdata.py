# -*- coding: utf-8 -*-
"""
Created on Sun Apr  2 11:52:42 2023

@author: mwate
"""

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.visualization import make_lupton_rgb
from PIL import Image
from skimage import exposure


# Open the FITS files and create objects containing the image data in a numpy array
dataf470n = fits.open('jw02733-o001_t001_nircam_f444w-f470n_i2d.fits') # f470 long wavelength data
dataf444w = fits.open('jw02733-o001_t001_nircam_f405n-f444w_i2d.fits') # f444w long wavelength data
dataf187n = fits.open('jw02733-o001_t001_nircam_clear-f187n_i2d.fits') # f187n short wavelength data
dataf356w = fits.open('jw02733-o001_t001_nircam_clear-f356w_i2d.fits') # f356w medium wavelength data
dataf090w = fits.open('jw02733-o001_t001_nircam_clear-f090w_i2d.fits') # f090w short wavelength data
image_dataf444w = dataf444w[1].data # Retrieve image data from index of FITS files
image_dataf356w = dataf356w[1].data
image_dataf187n = dataf187n[1].data
image_dataf090w = dataf090w[1].data
image_dataf470n = dataf470n[1].data

# resize image data and save as array
def resize(image): # define image resizing function
    image_resized = Image.fromarray(image)
    image_resized = image_resized.resize((2348, 2355))
    image_resized = np.asarray(image_resized)
    return image_resized.transpose()    
    
# pass image data into resizing function
image_dataf187n_resized = resize(image_dataf187n)
image_dataf356w_resized = resize(image_dataf356w)
image_dataf090w_resized = resize(image_dataf090w)


# normalization (0-255)
def norm(image): # define image normalizing function
    image_norm = ((image - image.min())/(image.max() - image.min()))*255
    return image_norm

# pass image data into normalization function
image_dataf187n_norm = norm(image_dataf187n_resized)
image_dataf356w_norm = norm(image_dataf356w_resized)
image_dataf444w_norm = norm(image_dataf444w)
image_dataf090w_norm = norm(image_dataf090w_resized)
image_dataf470n_norm = norm(image_dataf470n)

#%%Create false color composites
# plt.close('all')
# composite1 = make_lupton_rgb(image_dataf470n_norm, image_dataf356w_norm, image_dataf187n_norm, minimum=[10,5,8], stretch=2)
# plt.imshow(composite1)

# Stack images
stack = image_dataf470n_norm + image_dataf356w_norm + image_dataf187n_norm
stack = norm(stack)
plt.figure()
plt.imshow(stack, vmin=19.8, vmax=20.09)

stack_eq = exposure.equalize_hist(stack)
plt.figure()
plt.imshow(stack_eq, cmap='magma', vmin=0.66, vmax=1)
plt.colorbar()
#%%
stackflat = stack_eq.flatten()
plt.figure()
plt.hist(stackflat, bins='auto')

# plt.figure()
# cp1flat = composite1.flatten()
# plt.hist(cp1flat, bins='auto')

# plt.figure()
# stackflat = stack.flatten()
# plt.hist(stackflat, bins='auto')

# composite2 = make_lupton_rgb(image_dataf444w_norm, image_dataf444w_norm, image_dataf356w_norm, Q=10, stretch=0.05)
# plt.figure()
# plt.imshow(composite2)

# composite3 = make_lupton_rgb(image_dataf356w_norm, image_dataf444w_norm, image_dataf187n_norm, stretch=0.05)
# plt.figure()
# plt.imshow(composite3)

# composite4 = make_lupton_rgb(image_dataf444w_norm, image_dataf356w_norm, image_dataf090w_norm, minimum=5, stretch=90)
# plt.figure()
# plt.imshow(composite4)

