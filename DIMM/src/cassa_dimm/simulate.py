import numpy as np
from astropy.io import fits
from astropy.modeling.functional_models import Gaussian2D
import os

def generate_dimm_fits_cube(filename="dimm_simulated_cube.fits", frames=500):
    print(f"Simulating {frames} frames of raw DIMM FITS data...")
    
    # Sensor ROI dimensions (e.g., 64x128 pixels to fit both spots tightly)
    y_size, x_size = 64, 128
    
    # Base sub-pixel coordinates for the two spots
    # Simulating a prism that splits the beam into left and right spots
    x1_base, y1_base = 32.5, 32.0 
    x2_base, y2_base = 96.5, 32.0 
    
    # Telescope / Star characteristics
    fwhm = 3.5  # Base size of the star on the sensor
    stddev = fwhm / 2.355
    amplitude = 15000.0  # Peak ADU of the stars
    sky_bg = 300.0       # Sky background ADU
    read_noise = 10.0    # Camera read noise
    
    # Simulate atmospheric seeing variance (the "wobble")
    # This dictates how much the spots move relative to each other
    seeing_jitter_std = 1.2 # pixels
    
    # Pre-allocate the 3D numpy array for the FITS cube
    datacube = np.zeros((frames, y_size, x_size), dtype=np.float32)
    
    # Create the base coordinate grid
    y, x = np.mgrid[0:y_size, 0:x_size]
    
    for i in range(frames):
        # Add random atmospheric jitter to the coordinates
        x1 = x1_base + np.random.normal(0, seeing_jitter_std)
        y1 = y1_base + np.random.normal(0, seeing_jitter_std)
        x2 = x2_base + np.random.normal(0, seeing_jitter_std)
        y2 = y2_base + np.random.normal(0, seeing_jitter_std)
        
        # Model the two stars
        star1 = Gaussian2D(amplitude=amplitude, x_mean=x1, y_mean=y1, x_stddev=stddev, y_stddev=stddev)
        star2 = Gaussian2D(amplitude=amplitude, x_mean=x2, y_mean=y2, x_stddev=stddev, y_stddev=stddev)
        
        clean_frame = star1(x, y) + star2(x, y) + sky_bg
        
        # Apply Poisson noise (photon shot noise) and Gaussian read noise
        noisy_frame = np.random.poisson(clean_frame) + np.random.normal(0, read_noise, clean_frame.shape)
        
        datacube[i, :, :] = noisy_frame.astype(np.float32)

    # Save as a single multi-extension FITS file
    hdu = fits.PrimaryHDU(datacube)
    
    # Add essential metadata
    hdu.header['EXPTIME'] = 0.005  # 5 millisecond exposure
    hdu.header['INSTRUME'] = 'SIM_DIMM_8INCH'
    hdu.header['FRAMES'] = frames
    
    if os.path.exists(filename):
        os.remove(filename)
    hdu.writeto(filename)
    print(f"-> Successfully saved high-speed ROI datacube to {filename}")

if __name__ == "__main__":
    generate_dimm_fits_cube()