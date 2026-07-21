import numpy as np
from astropy.io import fits
from scipy.ndimage import median_filter
from photutils.centroids import centroid_com
import warnings

# Ignore benign photutils warnings about background levels
warnings.filterwarnings('ignore')

def calculate_dimm_from_fits(filename="dimm_simulated_cube.fits"):
    print(f"Loading FITS datacube: {filename}")
    
    # Load the 3D Datacube [Frames, Y, X]
    with fits.open(filename) as hdul:
        datacube = hdul[0].data
    
    num_frames, y_size, x_size = datacube.shape
    print(f"Loaded {num_frames} frames of shape {x_size}x{y_size}.")

    # --- HARDWARE SPECIFICATIONS FOR 8-INCH NEWTONIAN ---
    # You must calibrate these to your exact setup!
    pixel_size_um = 3.75        # e.g., ZWO ASI120/224 pixel size
    focal_length_mm = 1000.0    # 8-inch 200p f/5 focal length
    D = 0.06                    # Mask hole diameter in meters (60mm)
    d = 0.15                    # Distance between mask hole centers in meters (150mm)
    wavelength = 500e-9         # Standard reference wavelength (500nm)
    
    # Calculate Plate Scale
    plate_scale = 206.265 * (pixel_size_um / focal_length_mm) # arcsec/pixel
    arcsec_to_rad = 4.8481368e-6
    scale_rad = plate_scale * arcsec_to_rad

    # Arrays to store the differential distances
    dx_list = []
    dy_list = []

    print("Centroiding stars across all frames...")
    # To prevent the centroiding algorithm from getting confused, 
    # we split the frame in half: Left (Spot 1) and Right (Spot 2)
    midpoint_x = x_size // 2

    for i in range(num_frames):
        raw_frame = datacube[i]
        
        # --- NEW STEP: DESPIKING ---
        # Apply a 3x3 spatial median filter. 
        # This replaces every pixel with the median of its 8 neighbors.
        # Single "hot pixels" are instantly erased, but the multi-pixel star remains untouched.
        clean_frame = median_filter(raw_frame, size=3)
        
        # Sub-divide the cleaned image
        left_half = clean_frame[:, :midpoint_x]
        right_half = clean_frame[:, midpoint_x:]
        
        # Subtract local background before centroiding for sub-pixel accuracy
        left_bg_sub = np.maximum(left_half - np.median(left_half), 0)
        right_bg_sub = np.maximum(right_half - np.median(right_half), 0)
        
        # Calculate Center of Mass (COM) for both spots
        x1, y1 = centroid_com(left_bg_sub)
        x2_rel, y2 = centroid_com(right_bg_sub)
        
        # Adjust x2 back to global coordinates
        x2 = x2_rel + midpoint_x
        
        # Calculate differential distances
        dx_list.append(x1 - x2)
        dy_list.append(y1 - y2)

    # Convert pixel differences to Radians
    dx_rad = np.array(dx_list) * scale_rad
    dy_rad = np.array(dy_list) * scale_rad

    # Calculate Variances (sigma squared)
    var_l = np.var(dx_rad) # Longitudinal (parallel to baseline)
    var_t = np.var(dy_rad) # Transverse (perpendicular to baseline)

    # --- SARAZIN & RODDIER (1990) DIMM EQUATIONS ---
    # The coefficients Kl and Kt depend slightly on the ratio of D/d.
    # For D=60mm and d=150mm, the standard approximations hold well.
    Kl = 0.358 
    Kt = 0.242 
    
    # Calculate Fried Parameter (r0) in meters
    r0_l = (Kl * (wavelength**2) / (var_l * (d**(1/3))))**(3/5)
    r0_t = (Kt * (wavelength**2) / (var_t * (d**(1/3))))**(3/5)
    
    r0_avg = (r0_l + r0_t) / 2.0
    
    # Calculate Astronomical Seeing (FWHM) in arcseconds
    seeing_rad = 0.98 * (wavelength / r0_avg)
    seeing_arcsec = seeing_rad / arcsec_to_rad

    print("\n--- DIMM ANALYSIS RESULTS ---")
    print(f"Plate Scale:          {plate_scale:.3f} arcsec/pixel")
    print(f"Longitudinal r0:      {r0_l * 100:.2f} cm")
    print(f"Transverse r0:        {r0_t * 100:.2f} cm")
    print(f"Average Fried (r0):   {r0_avg * 100:.2f} cm")
    print(f"ESTIMATED SEEING:     {seeing_arcsec:.2f} arcseconds")
    print("-----------------------------")

if __name__ == "__main__":
    calculate_dimm_from_fits()