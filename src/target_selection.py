import astropy 
import astroplan 
from dataclasses import dataclass 
import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.dates as mdates
import datetime as dt 
import pytz 
import re 

# from . import aperture 
# from . import calibration 
# from . import lightcurve 
# from . import misc 
# from . import plate_solve 
# from . import plotting 
# from . import target_selection 







# Get Observer object used to calculate when targets are visible
def get_observer(lat=41.66, lon=-91.53, height=200, tz="US/Central"):
    location = astropy.coordinates.EarthLocation(
        lat=lat*astropy.units.deg,
        lon=lon*astropy.units.deg,
        height=height*astropy.units.m)
    return astroplan.Observer(location=location, timezone=tz)





# Get FixedTarget object corresponding to some (RA, Dec) position on the sky 
def get_target(dataframe_row=None, coord_deg_tuple=None, coord_str=None): 

    # Pass only one input. Raise an error if none or more than one are provided. 
    if sum(x is not None for x in [dataframe_row, coord_deg_tuple, coord_str]) != 1: 
        raise ValueError("Provide one and only one input to get_target()") 
    
    # If input is a row of the dataframe 
    if dataframe_row is not None: 
        target_coord = astropy.coordinates.SkyCoord(dataframe_row["Coords"], unit=(astropy.units.hourangle, astropy.units.deg)) 

    # If input is a tuple of the form (ra, dec), both in degrees 
    if coord_deg_tuple is not None: 
        target_coord = astropy.coordinates.SkyCoord(
        ra=coord_deg_tuple[0]*astropy.units.deg, 
        dec=coord_deg_tuple[1]*astropy.units.deg) 

    # If input is a string of the form "21 28 24.56 +46 40 30.8" 
    if coord_str is not None: 
        target_coord = astropy.coordinates.SkyCoord(coord_str, unit=(astropy.units.hourangle, astropy.units.deg))

    target = astroplan.FixedTarget(coord=target_coord) 
    return target 





# Class that holds all visibility information, such as the duration the target is visible and its altitude vs time
@dataclass 
class Visibility: 
    observer: astroplan.Observer            # Observer object used to calc this visibility 
    target: astroplan.FixedTarget           # Target object used to calc this visibilty 
    date: astropy.time.core.Time            # Night of observation (observing can continue past midnight onto the next day)
    min_alt: float | int                    # Minimum altitude for a target to be counted as visible (degrees)

    rise_time: astropy.time.core.Time       # Time the target rises (target altitude > minimum altitude) (observing can start if the sun has set)
    set_time: astropy.time.core.Time        # Time the target sets (target altitude < minimum altitude) (observing must end)
    
    sunset: astropy.time.core.Time          # Time the sun sets (observing can start if the target is visible) 
    sunrise: astropy.time.core.Time         # Time the sun rises (observing must end)
    
    start: astropy.time.core.Time           # Time observing can start (target must be above horizon and it must be night)
    end: astropy.time.core.Time             # Time observing must end (either target has set or sun has risen)
    duration: dt.timedelta                  # Duration that the target is visible 
    
    time_arr: np.array                      # Array of linearly spaced points in time (used in altitude vs time plot)
    alt_arr: np.array                       # Array of altitudes (heights of target above horizon in deggres)





# Calculate all aspects of visibility 
def calc_visibility(observer, target, date, min_alt=30):
    
    rise_time = observer.target_rise_time(date, target, which='next', horizon=min_alt*astropy.units.deg)
    set_time  = observer.target_set_time(date, target, which='next', horizon=min_alt*astropy.units.deg)

    sunset = observer.sun_set_time(date, which='next')
    sunrise = observer.sun_rise_time(date, which='next')

    time_arr = np.linspace(sunset, sunrise)
    alt_arr = observer.altaz(time_arr, target).alt 

    # Currently visible 
    if set_time<rise_time: 

        # Sets before sunset 
        if set_time<sunset: 

            if rise_time<sunset: 
                start=sunset 
                end=sunrise 
            if sunset<rise_time<sunrise: 
                start=rise_time 
                end=sunrise 
            if rise_time>sunrise: 
                start=sunrise 
                end=sunrise 
            
        # Sets during the night 
        if sunset<set_time<sunrise: 
            start=sunset 
            end=set_time 
        
        # Sets after sunrise 
        if set_time>sunrise: 
            start=sunset 
            end=sunrise 

    # Currently not visible 
    if set_time>rise_time: 

        # Rises before sunset 
        if rise_time<sunset: 

            if set_time<sunset: 
                start=sunset 
                end=sunset 
            if sunset<set_time<sunrise: 
                start=sunset 
                end=set_time 
            if set_time>sunrise: 
                start=sunset 
                end=sunrise 

        # Rises during the night 
        if sunset<rise_time<sunrise: 

            if set_time<sunrise: 
                start=rise_time 
                end=set_time 
            if set_time>sunrise: 
                start=rise_time
                end=sunrise
        
        # Rises after sunrise 
        if rise_time>sunrise: 
            start = sunrise  
            end = sunrise 

    try: 
        duration = end.to_datetime() - start.to_datetime() 

    # Error: If always up or never up 
    except UnboundLocalError:  

        # Always up 
        if np.nanmax(alt_arr.deg)>min_alt: 
            start = sunset 
            end = sunrise  

        # Never up 
        if np.nanmax(alt_arr.deg)<min_alt: 
            start = sunset 
            end = sunset 
    
    duration = end.to_datetime() - start.to_datetime() 

    visibility = Visibility(
        observer = observer, 
        target = target, 
        date = date, 
        min_alt = min_alt, 

        rise_time = rise_time, 
        set_time = set_time, 
        
        sunset = sunset, 
        sunrise = sunrise, 
        
        start = start, 
        end = end, 
        duration = duration, 
        
        time_arr = time_arr, 
        alt_arr = alt_arr, 
    )
    return visibility 
    




# Plot altitude vs time using the visibility information calculated with calc_visibility() 
def plot_target_altitude(visibility, target_row=None): 

    central = pytz.timezone("US/Central")

    # Create strings to use in title and legend label 
    date_str = f"{visibility.sunset.to_datetime(timezone=central).date()}" 
    if target_row is None: 
        title_str = f"Visibility of target from Iowa City on {date_str}"
        legend_label = "Target"
    if target_row is not None: 
        if len(target_row["Name"])>12: 
            target_str=f"{target_row['Name']:.12s}..."
        else: 
            target_str=f"{target_row['Name']}"
        title_str = f"Visibility of {target_str} from Iowa City on {date_str}"
        period_str = f"period={target_row['period_days']*24:.2f} hours"
        mag_diff_str = f"min flux={target_row['relative_flux_min']*100:.0f}%" 
        legend_label = f"{target_row['Name']} ({period_str}, {mag_diff_str})"

    plt.figure(figsize=(12, 6))
    plt.plot(visibility.time_arr.to_datetime(timezone=central), visibility.alt_arr, lw=3, label=legend_label) 

    # Add vertical lines at sunset/sunrise; shaded region to represent the visibility window; horizontal line at minimum altitude 
    plt.axhline(visibility.min_alt, ls="dotted", color="red", label=f"Minimum altitude = {visibility.min_alt} deg") 
    plt.axvspan(visibility.start.to_datetime(timezone=central), visibility.end.to_datetime(timezone=central), label=f"Visibility window ({visibility.duration.total_seconds()/3600:.2f} hours)", color="limegreen", alpha=0.1)

    datetime_obj = visibility.sunset.to_datetime(timezone=central)
    label = f"Sunset ({datetime_obj.time().strftime('%H:%M:%S')})"
    plt.axvline(datetime_obj, color="cornflowerblue", label=label)

    datetime_obj = visibility.sunrise.to_datetime(timezone=central)
    label = f"Sunrise ({datetime_obj.time().strftime('%H:%M:%S')})"
    plt.axvline(datetime_obj, color="gold", label=label)

    plt.xlabel("Local Time")
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M', tz=central))
    plt.ylabel("Altitude [deg]")
    plt.ylim((0, 90))
    plt.title(title_str) 
    plt.grid(alpha=0.5) 
    plt.legend() 





# Add a magnitude difference column to the dataframe 
def extract_magdiff(mag_str): 

    paren_match = re.search(r'\(([-+]?\d*\.\d+|\d+)\)', mag_str)

    # Case 1: Check if there's a " - " indicating two numbers to subtract
    if " - " in mag_str:
        numbers = re.findall(r'[-+]?\d*\.\d+|\d+', mag_str)
        if len(numbers) >= 2:
            magdiff = np.abs(float(numbers[0]) - float(numbers[1])) 
        else:
            magdiff = None  # Not enough numbers to subtract
        
    # Case 2: Check if there's a number inside parentheses
    elif paren_match:
        magdiff = float(paren_match.group(1))

    # Otherwise, no difference info found
    else: 
        magdiff = None
    
    return magdiff 





# Add a magnitude difference column to the dataframe 
def extract_magmax(mag_str):

    paren_match = re.search(r'\(([-+]?\d*\.\d+|\d+)\)', mag_str)

    # Case 1: Check if there's a " - " indicating two numbers to subtract
    if " - " in mag_str:
        numbers = re.findall(r'[-+]?\d*\.\d+|\d+', mag_str)
        numbers_float = [float(n) for n in numbers]
        magmax = np.max(numbers_float) 

    # Case 2: Check if there's a number inside parentheses
    elif paren_match:
        try: 
            magmax = float(mag_str.split(" ")[0].strip()) 
        except: 
            magmax = None
    
    else: 
        magmax = None 

    return magmax 




# Add declination in degrees column to the dataframe 
def extract_dec_degrees(s):
    match = re.search(r'([+-]\d+)', s)
    if match:
        return int(match.group(1))
    else:
        return None

