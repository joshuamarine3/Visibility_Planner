import scipy 
import astropy 
import twirl 
import pandas as pd 
import os 
import matplotlib.pyplot as plt 

# from . import aperture 
# from . import calibration 
# from . import lightcurve 
# from . import misc 
# from . import plate_solve 
from . import plotting 
# from . import target_selection 





# Identify stars in an image 
def find_peaks(image, num=15): 
    smoothed_image = scipy.ndimage.median_filter(image, size=3)
    stars_pixel_coords = twirl.find_peaks(smoothed_image, threshold=4.0)[0:num]
    return stars_pixel_coords 





# Query nearby stars to a given RA, Dec position from GAIA
def find_nearby_GAIA_stars(center, num=25): 
    pixel = 0.63*astropy.units.arcsec  # known pixel scale
    fov = 3072 * pixel.to(astropy.units.deg)
    all_GAIA_stars = twirl.gaia_radecs(center, 1.2*fov)
    return all_GAIA_stars[0:num]





# Plate solve one image 
def plate_solve_one_file(fullpath, stars_sky_coords): 

    # Load fits file 
    file = plotting.load_fits(fullpath)
    image = file.data 
    header = file.header 

    # Locate stars in image and nearby stars according to GAIA 
    stars_pixel_coords = find_peaks(image)
    wcs = twirl.compute_wcs(stars_pixel_coords, stars_sky_coords) 

    # Add WCS info to header and overwrite current file 
    wcs_header = wcs.to_header()
    header.update(wcs_header) 
    astropy.io.fits.writeto(fullpath, image, header, overwrite=True) 
    return wcs 





# Plate solve an entire folder 
def plate_solve_folder(folder, image_center, df_save_fullpath=None): 

    fullpaths = [folder/filename for filename in os.listdir(folder)]
    stars_sky_coords = find_nearby_GAIA_stars(image_center)

    # Re-order fullpaths so you do L1, B1, V1, L2, B2, V2, etc 
    order = {"L": 0, "B": 1, "V": 2}
    fullpaths = sorted(
        fullpaths,
        key=lambda f: (
            int(f.stem.split("-")[1][0:4]),      
            order[f.stem[-1]]                    
        )
    )

    x = [] 
    y = [] 

    for fullpath in fullpaths: 
        
        wcs = plate_solve_one_file(fullpath, stars_sky_coords) 

        # Add x/y coordinates of image center to csv 
        skycoord = astropy.coordinates.SkyCoord(
            ra=image_center[0]*astropy.units.deg, 
            dec=image_center[1]*astropy.units.deg, frame='icrs') 
        position = tuple([float(i) for i in wcs.world_to_pixel(skycoord) ])
        x.append(position[0])
        y.append(position[1])

    if df_save_fullpath is not None: 
        data = {
            "pixel_x": x,
            "pixel_y": y,
            "fullpaths": fullpaths}
        df = pd.DataFrame(data)
        df.to_csv(df_save_fullpath, index=False)





def plot_plate_solve(fullpath, image_center, aperture_stars): 

    file = plotting.load_fits(fullpath)
    wcs = astropy.wcs.WCS(file.header)

    plotting.plot_image(file.data) 
    plt.title(f"{fullpath}", fontsize=20)

    # Peaks in image (RED)
    stars_pixel_coords = find_peaks(file.data)  
    plt.scatter(
        stars_pixel_coords[:,0], stars_pixel_coords[:,1], 
        marker="o", facecolors="none", edgecolors="red", s=50, lw=2, label="Stars identified by looking for peaks in image") 

    # Stars targeted for aperture photometry (GOLD) 
    for ra_dec_pair in aperture_stars: 
        skycoord = astropy.coordinates.SkyCoord(
            ra=ra_dec_pair[0]*astropy.units.deg, 
            dec=ra_dec_pair[1]*astropy.units.deg, frame='icrs') 
        pixel_coord = tuple([float(i) for i in wcs.world_to_pixel(skycoord)]) 
        label = None
        if ra_dec_pair == aperture_stars[0]: 
            label="Stars targeted for aperture photometry"
        plt.scatter([
            pixel_coord[0]], [pixel_coord[1]], 
            marker="o", facecolors="none", edgecolors="gold", s=200, lw=2, label=label) 

    # GAIA stars (GREEN) 
    stars_sky_coords = find_nearby_GAIA_stars(image_center) 
    for i in range(len(stars_sky_coords)): 
        ra_dec_pair = stars_sky_coords[i]  
        skycoord = astropy.coordinates.SkyCoord(
            ra=ra_dec_pair[0]*astropy.units.deg, 
            dec=ra_dec_pair[1]*astropy.units.deg, frame='icrs') 
        pixel_coord = tuple([float(i) for i in wcs.world_to_pixel(skycoord)]) 
        label = None 
        if i == 0: 
            label = "GAIA star positions according to plate solve"
        plt.scatter([
            pixel_coord[0]], [pixel_coord[1]],  
            marker="o", facecolors="none", edgecolors="limegreen", s=500, lw=2, label=label) 

    plt.legend() 
    plt.xlim((0, 3072))
    plt.ylim((0, 2047))
    



