import numpy as np 
import os 
import astropy 
import matplotlib.pyplot as plt 
import scipy

# from . import aperture 
# from . import calibration 
# from . import lightcurve 
from . import misc 
# from . import plate_solve 
from . import plotting 
# from . import target_selection 





# Median combine all files in a folder 
def create_master_file(folder, master_filename_nofilter, filter="", master_bias=None, normalize=False, overwrite=False): 
    print("\n")
    fullpaths = [folder / filename for filename in os.listdir(folder) if filter in filename]

    # Add filter to filename 
    if filter != "": 
        base, ext = master_filename_nofilter.rsplit(".fits", 1) # Split from the right, only once 
        master_filename = base + f"_{filter}.fits" + ext  # 'ext' will be whatever comes after, usually empty
    else: 
        master_filename = master_filename_nofilter 
    
    # Escape if master file already calcualted 
    master_fullpath = folder / master_filename 
    if master_fullpath.exists() and overwrite==False: 
        print(f"{master_fullpath} already exists") 
        return 
    
    # Warn the user if overwriting a file 
    if master_fullpath.exists() and overwrite==True: 
        print(f"Overwriting {master_fullpath}")
    
    print(f"Combining {len(fullpaths)} files")
    
    # Loop over all files in folder 
    image_list = [] 
    for fullpath in fullpaths: 

        image = misc.load_fits(fullpath).data 

        # Subtract master bias from each frame as its loaded (if master_bias is provided)
        if master_bias is not None: 
            image = image - master_bias.data 

        image_list.append(image)

    # Median combine 
    master_image = np.median(image_list, axis=0)

    # Normalize 
    if normalize == True: 
        master_image = master_image / np.max(master_image) 

    # Copy header from the first file, add a comment 
    header = plotting.load_fits(fullpaths[0]).header 
    header['COMMENT'] = f"Master image, created by median combining {len(fullpaths)} images"

    # Save master image to new fits file with updated header 
    hdu = astropy.io.fits.PrimaryHDU(master_image, header=header)
    hdu.writeto(master_fullpath, overwrite=True) 
    print(f"Saving file to: {master_fullpath}")





# Calibrate a single image 
def calibrate_one_file(raw_fullpath, calibrated_fullpath, master_bias, master_dark, master_flat, overwrite=False): 

    # Escape if file already calcualted 
    if calibrated_fullpath.exists() and overwrite==False: 
        print(f"{calibrated_fullpath} already exists")
        return  
    
    # Warn the user if overwriting a file 
    if calibrated_fullpath.exists() and overwrite==True: 
        print(f"Overwriting {calibrated_fullpath}")
    
    # Load fits files of raw image 
    raw_file = plotting.load_fits(raw_fullpath)

    # Calculate calibrated image 
    calibrated_image = ( raw_file.data - master_bias.data - master_dark.data ) / master_flat.data  

    # Save calibrated image to new fits file
    hdu = astropy.io.fits.PrimaryHDU(calibrated_image, header=raw_file.header) 
    hdu.writeto(calibrated_fullpath, overwrite=True) 
    print(f"Saving calibrated image to: {calibrated_fullpath}") 





# Calibrate all images in a folder with a specific filter 
def calibrate_all_files_one_filter(raw_folder, calibrated_folder, filter, bias_folder, dark_folder, flat_folder): 

    # Load bias, dark, and flat FITS files that apply to all files in this folder 
    master_bias = plotting.load_fits(plotting.find_fullpath(bias_folder, ["master"]))
    master_dark = plotting.load_fits(plotting.find_fullpath(dark_folder, ["master", filter]))
    master_flat = plotting.load_fits(plotting.find_fullpath(flat_folder, ["master", filter]))

    raw_fullpaths = [raw_folder / filename for filename in os.listdir(raw_folder) if filter in filename] 
    for raw_fullpath in raw_fullpaths: 
        
        print("\n")
        calibrated_fullpath = calibrated_folder / raw_fullpath.name.replace("Raw", "Calibrated") 
        calibrate_one_file(raw_fullpath, calibrated_fullpath, master_bias, master_dark, master_flat)





# Remove glare from a calibrated image (used when flats didn't work)
def remove_glare_one_file(fullpath, show_image=False, save=True): 
    file1 = plotting.load_fits(fullpath)
    print("\n")     
    print(f"Removing glare from: {fullpath}")
    background_blocky = scipy.ndimage.grey_opening(file1.data, size=(20,20))
    background_smooth = scipy.ndimage.gaussian_filter(background_blocky, sigma=60)  
    corrected_image = file1.data - background_smooth 
    
    if save==True: 
        hdu = astropy.io.fits.PrimaryHDU(corrected_image, header=file1.header)
        hdu.writeto(fullpath, overwrite=True) 
        print(f"Overwriting {fullpath}")

    if show_image==True: 

        plotting.plot_image(background_smooth) 
        plt.title(f"Approximated glare \n{fullpath}", fontsize=20)
        
        plotting.plot_image(file1.data)
        plt.title(f"Original image \n{fullpath}", fontsize=20)

        plotting.plot_image(corrected_image)
        plt.title(f"Glare removed \n{fullpath}", fontsize=20)




# Remove glare from all images in a folder 
def remove_glare_folder(folder): 
    fullpaths = [folder/filename for filename in os.listdir(folder)]
    for fullpath in fullpaths: 
        remove_glare_one_file(fullpath) 

