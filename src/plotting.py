import numpy as np
from typing import Optional, Dict, Any, Tuple
import pytz 

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec
import matplotlib.patheffects as path_effects

from . import aperture 
# from . import calibration 
from . import lightcurve 
from . import misc 
# from . import plate_solve 
# from . import plotting 
# from . import target_selection 





# Helper: compute proper imshow limits for integer pixel indexing
def _imshow_xlim_ylim_from_center(center, window, shape):
    """Return (x0,x1), (y0,y1) in data coordinates suitable for imshow axes.
    center: (xpix, ypix) in pixel coordinates (0..N-1)
    window: side length in pixels (float or int)
    shape: (ny, nx) image shape
    """
    nx = shape[1]
    ny = shape[0]
    cx, cy = float(center[0]), float(center[1])
    half = float(window) / 2.0
    # imshow pixel centers go from 0 to N-1; imshow extent default is [-0.5, N-0.5]
    x0 = cx - half - 0.5
    x1 = cx + half - 0.5
    y0 = cy - half - 0.5
    y1 = cy + half - 0.5
    # Respect image boundaries
    x0 = max(x0, -0.5)
    x1 = min(x1, nx - 0.5)
    y0 = max(y0, -0.5)
    y1 = min(y1, ny - 0.5)
    return (x0, x1), (y0, y1)




def plot_image_ax(
        image: np.ndarray,
        ax: Optional[plt.Axes] = None,
        fig: Optional[plt.Figure] = None,
        vrange_std: Optional[tuple]=None,
        vrange_absval: Optional[tuple]=None,
        center_position: Optional[Tuple[float,float]]=None,
        window_size: int = 200,
        cmap='gray', 
        obs_time = None, 
    ) -> Dict[str, Any]:
    """Plot image on ax. Returns dict with fig, ax, im, cbar (or None)."""
    
    if ax is None:
        fig = plt.figure(figsize=(9, 6))
        ax = fig.add_axes([0, 0, 1, 1])  # full figure
    else:
        fig = ax.figure 

    
    if vrange_std is None and vrange_absval is None:
        vrange_std = (-1, 4)
    if vrange_std is not None and vrange_absval is not None:
        raise ValueError("Provide either vrange_std or vrange_absval, not both")

    mean = np.mean(image)
    std = np.std(image)
    if vrange_std is not None:
        vmin = mean + min(vrange_std) * std
        vmax = mean + max(vrange_std) * std
    else:
        vmin, vmax = min(vrange_absval), max(vrange_absval)

    im = ax.imshow(image, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)

    # cbar = None
    # if show_colorbar:
    #     divider = make_axes_locatable(ax)
    #     cax = divider.append_axes("right", size="4%", pad=0.05)
    #     colorbar_kwargs = colorbar_kwargs or {}
    #     cbar = fig.colorbar(im, cax=cax, **colorbar_kwargs)

    # handle center cropping robustly using pixel coords
    if center_position is not None:
        (x0, x1), (y0, y1) = _imshow_xlim_ylim_from_center(center_position, window_size, image.shape)
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        # add a little margin so we don't clip text drawn outside axis
        ax.margins(x=0)
    else:
        # ensure axis covers whole image in pixel coordinates
        ax.set_xlim(-0.5, image.shape[1]-0.5)
        ax.set_ylim(-0.5, image.shape[0]-0.5)

        # Set ticks 
    ax.set_xticks([t for t in ax.get_xticks() if t>ax.get_xlim()[0]+1 and t<ax.get_xlim()[1]-1])
    ax.set_yticks([t for t in ax.get_yticks() if t>ax.get_ylim()[0]+1 and t<ax.get_ylim()[1]-1])

    # Set tick parameters 
    for label in ax.get_yticklabels():
        label.set_rotation(90)
        label.set_ha('center')   # horizontal alignment
        label.set_va('center')   # vertical alignment
    ax.tick_params(direction='in', length=8, width=1, color="white", pad=-15, labelcolor="white", labelsize=9) 

    for tick in ax.get_xticklabels() + ax.get_yticklabels():
        tick.set_path_effects([
            path_effects.Stroke(linewidth=0.8, foreground='tomato'),
            path_effects.Normal()
        ])

    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

    # Add observation time to plot  
    if obs_time is not None: 
        t0 = obs_time 
        t1 = t0.tz_localize("UTC")
        t2 = t1.tz_convert("US/Central")
        title_str = t2.strftime("%H:%M:%S")
        text = ax.text(
            0.92, 0.96, title_str, color='white',
            ha='center', va='center', transform=ax.transAxes, fontsize=20)
        text.set_path_effects([
                path_effects.Stroke(linewidth=1.5, foreground='lime'),
                path_effects.Normal()
            ])
    
    return dict(fig=fig, ax=ax, im=im)#, cbar=cbar)





def draw_apertures_on_ax(
        ax: plt.Axes, 
        apertures: list[aperture.ApertureStarResult] = [], 
        label_style="verbose", # "verbose", "minimal", or None 
    ):
    
    artists = []
    num_labels_placed = 0 
    for atr in apertures:

        # plot aperture + annulus objects (photutils .plot returns artists)
        a_art = atr.aperture.plot(ax=ax, color=atr.star.color, lw=2)
        ann_art = atr.annulus.plot(ax=ax, color=atr.star.color, lw=2)
        artists.extend([a_art, ann_art])

        lines = [
            f"{atr.star.name}",
            f"{atr.star.label}",
            f"Source counts = {atr.source_counts/1000:.3g}k",
            f"Total counts = {atr.total_counts/1000:.3g}k",
            f"Source / total = {atr.source_counts/atr.total_counts*100:.3g}%",
            f"SNR = {atr.SNR:.3g}",
        ]

        if label_style == "minimal":
            info_text = "\n".join(lines[:2])   # first two lines only
        elif label_style == "verbose":
            info_text = "\n".join(lines)       # all lines

        if label_style is not None: 
            label_positions = [0.98, 0.80, 0.62]
            ax.text(
                0.05, label_positions[num_labels_placed], 
                info_text, transform=ax.transAxes,
                fontsize=9, verticalalignment="top",
                bbox=dict(boxstyle="round", edgecolor=atr.star.color, lw=4, facecolor="white", alpha=0.7)
            )
            num_labels_placed += 1 

    return artists





def plot_lightcurve_ax(
        lc: lightcurve.LightcurveData, 
        star: misc.Star, 
        filters = None, 
        ax: Optional[plt.Axes]=None, 
        fig: Optional[plt.Figure]=None,
        aperture_file_result=None, 
        show_legend=True, 
        include_smoothed_line=True, 
        include_title=True  

    ):
    if ax is None:
        fig, ax = plt.subplots(figsize=(8,6)) if fig is None else (fig, fig.add_subplot(111))
    else:
        fig = ax.figure

    if filters is None:
        filters = lc.filters

    artists = []
    ax.invert_yaxis()
    for filt in filters: 

        x1 = lc.get(filter=filt, star=star, grid="original", key="obs_time")
        y1 = lc.get(filter=filt, star=star, grid="original", key="mag")
        ax.scatter(x1, y1, s=10, color=filt.plot_color, label=f"{filt.name} filter", )
        
        if include_smoothed_line==True: 
            x2 = lc.get(filter=filt, star=star, grid="interped", key="obs_time")
            y2 = lc.get(filter=filt, star=star, grid="interped", key="mag")
            l, = ax.plot(x2, y2, lw=2, color=filt.plot_color, alpha=0.5, zorder=0.1)
            artists.append(l)

    if include_title: 
        ax.set_title(f"Lightcurve of {star.label} ({star.name})")
    ax.set_xlabel("Time of observation")
    ax.set_ylabel("Magnitude")
    if show_legend:
        ax.legend()
    ax.grid(alpha=0.5)
    central = pytz.timezone("US/Central")
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M', tz=central)) 

    # draw current-time line if requested
    current_line = None
    if aperture_file_result is not None:
        kwargs = dict(ls='--', color=star.color, zorder=0, lw=2)
        current_line = ax.axvline(aperture_file_result.obs_time, **kwargs)

    return dict(fig=fig, ax=ax, artists=artists, current_line=current_line)






def combined_frame_plot(
        file, 
        aperture_for_zoom: aperture.ApertureStarResult, 
        aperture_file_result: aperture.ApertureFileResult, 
        lc: lightcurve.LightcurveData, 
    ):
        
    fig = plt.figure(figsize=(9 + 0.1 + 6, 6 + 0.1 + 2.1 + 0.5))

    # 2 rows, 3 columns
    gs = GridSpec(
        nrows=4,
        ncols=4, 
        figure=fig,
        height_ratios=[6, 0.1, 2.1, 0.5], 
        width_ratios=[0.7, 8.3, 0.1, 6]
    )
    ax1 = fig.add_subplot(gs[0, 0:2])  
    ax2 = fig.add_subplot(gs[0, 3])    
    ax3 = fig.add_subplot(gs[2, 1:])

    # Full image
    res_full = plot_image_ax(file['data'], ax=ax1, vrange_std=(-1,4), obs_time=aperture_file_result.obs_time) 
    draw_apertures_on_ax(ax1, apertures=[aperture_for_zoom], label_style="minimal")

    # Zoomed image (if star provided)
    center = misc.radec_to_pixelcoords(file, aperture_for_zoom.star.position)
    res_zoom = plot_image_ax(file['data'], ax=ax2, center_position=center, window_size=200, vrange_std=(-1,4))
    draw_apertures_on_ax(ax2, apertures=[aperture_for_zoom], label_style=None)
    
    # Lightcurve
    filters = [filter for filter in lc.filters if filter.name[0] == f'{file["fullpath"].name[-5]}']
    res_lc = plot_lightcurve_ax(
        ax=ax3, 
        lc=lc, 
        star=aperture_for_zoom.star, 
        filters=filters, 
        aperture_file_result=aperture_file_result, 
        include_smoothed_line=False, 
        include_title=False 
    )

    fig.subplots_adjust(left=0, right=1, bottom=0, top=1, hspace=0, wspace=0)
    return dict(fig=fig, axes=[ax1, ax2, ax3], image_full=res_full, image_zoom=res_zoom, lc=res_lc)





# Raw counts (all 3 stars at one filter)
def plot_raw_counts(lc: lightcurve.LightcurveData, filter: misc.Filter, include_smoothed_line=True): 
    fig, ax = plt.subplots(figsize=(8, 6))

    for star in lc.stars: 
        x = lc.get(filter=filter, star=star, grid="original", key="obs_time") 
        y1 = lc.get(filter=filter, star=star, grid="original", key="counts") 
        y2 = lc.get(filter=filter, star=star, grid="original", key="counts_smoothed") 
        plt.scatter(x, y1, s=10, color=star.color, label=f"{star.label} ({star.name})") 
        if include_smoothed_line == True: 
            plt.plot(x, y2, lw=2, color=star.color, alpha=0.5)

    for spine in plt.gca().spines.values():
        spine.set_edgecolor(filter.plot_color)
        spine.set_linewidth(2)

    plt.title(f"Raw counts in {filter.name} filter")
    plt.xlabel("Time (Central time zone)")  
    plt.ylabel("Counts") 
    plt.legend() 
    plt.grid(alpha=0.5) 

    central = pytz.timezone("US/Central")
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M', tz=central)) 





# Normalize and convert to magnitude (all 3 stars at one filter)
def plot_normalized_magnitudes(lc: lightcurve.LightcurveData, filter: misc.Filter, include_smoothed_line=True): 
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.invert_yaxis() 

    for star in lc.stars: 
        
        x1 = lc.get(filter=filter, star=star, grid="original", key="obs_time") 
        y1 = lc.get(filter=filter, star=star, grid="original", key="mag") 
        plt.scatter(x1, y1, s=10, color=star.color, label=f"{star.label} ({star.name})") 

        if include_smoothed_line == True: 
            x2 = lc.get(filter=filter, star=star, grid="interped", key="obs_time") 
            y2 = lc.get(filter=filter, star=star, grid="interped", key="mag") 
            plt.plot(x2, y2, lw=2, color=star.color, alpha=0.5)

    for spine in plt.gca().spines.values():
        spine.set_edgecolor(filter.plot_color)
        spine.set_linewidth(2)

    plt.title(f"Normalized magnitudes in {filter.name} filter") 
    plt.xlabel("Time (Central time zone)")  
    plt.ylabel("Magnitude") 
    plt.legend() 
    plt.grid(alpha=0.5) 

    central = pytz.timezone("US/Central")
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M', tz=central)) 





def plot_color(lc: lightcurve.LightcurveData, star: misc.Star, B_filter: misc.Filter, V_filter: misc.Filter, include_smoothed_line=True): 

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.invert_yaxis() 

    if include_smoothed_line==True: 
        times = lc.get(filter=B_filter, star=star, grid="interped", key="obs_time") 
        B = lc.get(filter=B_filter, star=star, grid="interped", key="mag") 
        V = lc.get(filter=V_filter, star=star, grid="interped", key="mag") 
        B_minus_V = B-V
        plt.plot(times, B_minus_V, lw=2, color="black", alpha=0.5)

    times = lc.get(filter=B_filter, star=star, grid="original", key="obs_time") 
    B = lc.get(filter=B_filter, star=star, grid="original", key="mag") 
    V = lc.get(filter=V_filter, star=star, grid="original", key="mag") 
    B_minus_V = B-V
    plt.scatter(times, B_minus_V, s=20, color="gold", ec="black")

    plt.title(f"Color variation of {star.label} ({star.name})") 
    plt.xlabel("Time (Central time zone)")  
    plt.ylabel("B-V (color index) ") 
    # plt.legend() 
    plt.grid(alpha=0.5) 

    central = pytz.timezone("US/Central")
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M', tz=central)) 

    # Add blue-to-red gradient in the background of the plot 
    gradient = np.linspace(0, 1, 256).reshape(-1, 1) 
    current_xlim = ax.get_xlim() 
    current_ylim = ax.get_ylim()
    plt.imshow(
        gradient,
        extent=[current_xlim[0], current_xlim[1], current_ylim[0], current_ylim[1]],
        origin="lower",
        aspect="auto",
        cmap="coolwarm_r",
        alpha=0.25 
    ) 






