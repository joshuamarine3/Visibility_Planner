import numpy as np 
import photutils.aperture 
from dataclasses import dataclass 
import astropy 
import os 
import pandas as pd 
import matplotlib.pyplot as plt 

# from . import aperture 
# from . import calibration 
# from . import lightcurve 
from . import misc 
# from . import plate_solve 
# from . import plotting 
# from . import target_selection 




# Calculate the read noise by taking pairs of bias frames, subtracting them, and finding the standard deviation of the difference image 
def calc_read_noise(bias_folder): 

    read_noise_list = [] 
    fullpaths = [bias_folder/filename for filename in os.listdir(bias_folder) if "aster" not in filename]
    for i in range(len(fullpaths)-1): 

        file1 = misc.load_fits(fullpaths[i]) 
        file2 = misc.load_fits(fullpaths[i+1]) 
        
        data1 = file1["data"].astype(float)
        data2 = file2["data"].astype(float)
        diff = data1 - data2

        read_noise = np.std(diff) / np.sqrt(2) 
        read_noise_list.append(read_noise)
    
    read_noise = np.median(read_noise_list)
    return read_noise





# Class that holds information about a specific run of calc_aperture 
@dataclass 
class ApertureStarResult: 
    star: misc.Star 
    aperture: photutils.aperture.CircularAperture   # Aperture plotting object 
    annulus: photutils.aperture.CircularAnnulus     # Annulus plotting object 
    total_counts: float | int                       # Total counts in aperture 
    source_counts: float | int                      # Counts from just the source in the aperture 
    sky_counts: float | int                         # Counts from the background (or sky) in the aperture 
    n_pix: float | int                              # Number of pixels in the aperture 
    read_noise: float | int                         # Read noise 
    SNR: float | int                                # Signal-to-noise-ratio 





@dataclass 
class ApertureFileResult: 
    target: ApertureStarResult 
    norm: ApertureStarResult 
    test: ApertureStarResult 
    obs_time: float 
    filter: float 





def calc_aperture_one_star(star: misc.Star, file, read_noise): 

    # Load fits file 
    header = file["header"]
    data = file["data"]  
    wcs = astropy.wcs.WCS(header) 

    # Source position on image 
    target_skycoord = astropy.coordinates.SkyCoord(
        ra=star.position[0]*astropy.units.deg, 
        dec=star.position[1]*astropy.units.deg, frame='icrs')
    target_x_y = tuple([float(i) for i in wcs.world_to_pixel(target_skycoord) ])

    # Aperture 
    aperture = photutils.aperture.CircularAperture(target_x_y, r=star.r_aperture) 
    aperture_phot_table = photutils.aperture.aperture_photometry(data, aperture)
    aperture_counts = aperture_phot_table['aperture_sum'][0]
    
    # Annulus 
    r_annulus_in = star.r_annulus_in_relative*star.r_aperture
    r_annulus_out = 1.75*star.r_annulus_in_relative*star.r_aperture 
    annulus = photutils.aperture.CircularAnnulus(target_x_y, r_in=r_annulus_in, r_out=r_annulus_out)  
    annulus_phot_table = photutils.aperture.aperture_photometry(data, annulus)
    annulus_counts = annulus_phot_table['aperture_sum'][0]

    # Subtract background from original flux
    bkg_per_pixel = annulus_counts / annulus.area
    bkg_in_aperture = bkg_per_pixel * aperture.area
    corrected_counts = aperture_counts - bkg_in_aperture 

    # Calculate signal-to-noise-ratio (SNR) 
    S_star = corrected_counts 
    npix_times_S_sky = bkg_in_aperture
    npix_times_R_squared = aperture.area * read_noise**2 
    SNR = S_star / np.sqrt(S_star + npix_times_S_sky + npix_times_R_squared)

    # Store results in a custom dataclass 
    results = ApertureStarResult( 
        star = star, 
        aperture = aperture, 
        annulus = annulus, 
        total_counts = aperture_counts, 
        source_counts = corrected_counts, 
        sky_counts = aperture_counts - corrected_counts, 
        n_pix = aperture.area, 
        read_noise = read_noise, 
        SNR = SNR
    )
    return results 





def calc_aperture_one_file(file, read_noise, target_star: misc.Star, norm_star: misc.Star, test_star: misc.Star):

    results_target = calc_aperture_one_star(target_star, file, read_noise) 
    results_norm = calc_aperture_one_star(norm_star, file, read_noise) 
    results_test = calc_aperture_one_star(test_star, file, read_noise) 

    obs_time = pd.to_datetime(file["header"]["DATE-OBS"]) 
    filter = file["header"]["FILTER"]

    results = ApertureFileResult(
        target = results_target, 
        norm = results_norm, 
        test = results_test, 
        obs_time = obs_time, 
        filter = filter 
    )
    return results 





def calc_aperture_folder(bias_folder, calibrated_folder, target_star, norm_star, test_star): 

    read_noise = calc_read_noise(bias_folder) 
    
    aperture_folder_results = [] 

    fullpaths = [calibrated_folder/filename for filename in os.listdir(calibrated_folder)] 
    for fullpath in fullpaths: 
    
        file = misc.load_fits(fullpath)
    
        aperture_file_result = calc_aperture_one_file(file, read_noise, target_star, norm_star, test_star)
        aperture_folder_results.append(aperture_file_result)

    return aperture_folder_results 




