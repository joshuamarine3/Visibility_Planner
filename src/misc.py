import numpy as np 
from dataclasses import dataclass, field
import astropy 
import os 
from typing import Optional, Tuple, Dict  

# from . import aperture 
# from . import calibration 
# from . import lightcurve 
# from . import misc 
# from . import plate_solve 
# from . import plotting 
# from . import target_selection 



@dataclass(frozen=True)
class Filter:
    name: str
    plot_color: str





@dataclass(frozen=True)
class Star:
    name: str
    label: str
    color: str 
    id: str 
    position: Tuple[float, float]
    r_aperture: float
    r_annulus_in_relative: float = 4.0
    magnitudes: Dict[Filter, float] = field(
        default_factory=dict,
        compare=False,
        hash=False
    )





# Find the path to a file that contains a set of strings in its name 
def find_fullpath(folder, subset_strs_list):    
    fullpaths = [
        folder / filename
        for filename in os.listdir(folder)
        if all(subset_str in filename for subset_str in subset_strs_list)] 
    if fullpaths: 
        return fullpaths[0] 
    return None 





# Load a FITS file (include the data as a 2d array, the header, and the full path to the file that was loaded)
def load_fits(fullpath):
    print(f"Loading {fullpath}") 
    data= astropy.io.fits.getdata(fullpath, ext=0)  
    header = astropy.io.fits.open(fullpath)[0].header 
    file = {"data": data, "header": header, "fullpath": fullpath}
    return file





# Convert from RA and dec coordinates on the sky to pixel coordinates on a particular platesolved image 
def radec_to_pixelcoords(file, radec): 

    # Load fits file 
    header = file["header"]
    wcs = astropy.wcs.WCS(header) 

    # RA, dec -> position on image 
    target_skycoord = astropy.coordinates.SkyCoord(
        ra=radec[0]*astropy.units.deg, 
        dec=radec[1]*astropy.units.deg, frame='icrs')
    target_x_y = tuple([float(i) for i in wcs.world_to_pixel(target_skycoord) ]) 
    return target_x_y 





# Convert from a flux ratio to a magnitude difference
def fluxratio_to_magdiff(flux_ratio): 
    delta_m = 2.5*np.log10(flux_ratio) 
    return delta_m 

# Convert from a magnitude difference to a flux ratio 
def magdiff_to_fluxratio(magdiff): 
    flux_ratio = 10**(-magdiff/2.5) 
    return flux_ratio 







