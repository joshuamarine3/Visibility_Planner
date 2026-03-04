import numpy as np 
import pandas as pd 
import scipy 

# from . import aperture 
# from . import calibration 
# from . import lightcurve 
from . import misc 
# from . import plate_solve 
# from . import plotting 
# from . import target_selection 





class LightcurveData: 
    def __init__(self, aperture_folder_results, filter_list, target_star, norm_star, test_star): 
        
        self.aperture_folder_results = aperture_folder_results 

        # Define available values 
        self.filters = filter_list 
        self.stars = [target_star, norm_star, test_star] 
        self.grids = ["original", "interped"] 
        
        # Initialize a dict inisde a dict inside a dict 
        self.data = {} 
        for f in self.filters: 
            self.data[f] = {} 
            for s in self.stars: 
                self.data[f][s] = {} 
                for g in self.grids: 
                    self.data[f][s][g] = {} 
        


        # Interpolated times (same for all filters) 
        times_interped = self.retrieve_times() 

        # Loop over all 3 filters 
        for filter in self.filters: 

            # Stuff that only depends on filter (will be used for all 3 stars)
            norm_star_mag = norm_star.magnitudes[filter] 
            times = self.retrieve_times(filter) 
            for star in self.stars: 
                self.set(filter=filter, star=star, grid="original", key="obs_time", value=times)
                self.set(filter=filter, star=star, grid="interped", key="obs_time", value=times_interped)
            
            # Raw counts -> smooth -> interpolate 
            for star in self.stars: 
                if star == target_star: 
                    smoothing_window=15
                else: 
                    smoothing_window=30 
                counts = self.retrieve_counts(filter, star)
                counts_smoothed = scipy.ndimage.uniform_filter1d(counts, size=smoothing_window, mode="nearest") 
                counts_smoothed_interped = np.interp(times_interped, times, counts_smoothed) 
                self.set(filter=filter, star=star, grid="original", key="counts", value=counts)
                self.set(filter=filter, star=star, grid="original", key="counts_smoothed", value=counts_smoothed)
                self.set(filter=filter, star=star, grid="interped", key="counts_smoothed", value=counts_smoothed_interped)

            # Relative flux and magnitude by comparing to norm star 
            for star in self.stars: 

                # Include magnitude for both original grid and interpolated/smoothed grid 
                for (grid, counts_key) in [("interped", "counts_smoothed"), ("original", "counts")]: 

                    target_counts = self.get(filter=filter, star=star, grid=grid, key=counts_key) 
                    norm_counts = self.get(filter=filter, star=norm_star, grid=grid, key="counts_smoothed")  
                    relative_flux_target = target_counts / norm_counts  
                    mag_diff_target = misc.fluxratio_to_magdiff(relative_flux_target)
                    mag_target = norm_star_mag - mag_diff_target 
                    self.set(filter=filter, star=star, grid=grid, key="mag", value=mag_target)




        
    # Raise an error if you try to access an invalid value 
    def check_valid(self, filter, star, grid): 
        if filter not in self.filters: 
            raise NameError(f"Filter '{filter}' not recognized. Must be one of {self.filters}") 
        if star not in self.stars: 
            raise NameError(f"Star '{star}' not recognized. Must be one of {self.stars}") 
        if grid not in self.grids: 
            raise NameError(f"Grid '{grid}' not recognized. Must be one of {self.grids}") 



    # Set a value (used when initiallizing the object) 
    def set(self, filter, star, grid, key, value): 
        self.check_valid(filter, star, grid) 
        self.data[filter][star][grid][key] = value 
        


    # Return a value (makes syntax for calling data more intuitive) 
    def get(self, filter, star, grid, key): 
        self.check_valid(filter, star, grid) 
        return self.data[filter][star][grid].get(key, None) 
    


    
    # Helper: Retrieve times array from results list 
    def retrieve_times(self, filter: misc.Filter = None): 
        if filter is not None: 
            times = [
                r.obs_time
                for r in self.aperture_folder_results
                if r.filter == filter.name
            ]
        else: 
            times = [
                r.obs_time
                for r in self.aperture_folder_results
            ] 
        times = pd.to_datetime(times)
        return times 



    # Helper: Retrieve counts array from results list 
    def retrieve_counts(self, filter: misc.Filter, star: misc.Star): 
        counts = np.array([
            getattr(r, star.id).source_counts
            for r in self.aperture_folder_results
            if r.filter == filter.name
        ]) 
        return counts 






