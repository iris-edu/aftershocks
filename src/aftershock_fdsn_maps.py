#!/usr/bin/env python

import sys
import os

import glob
from datetime import datetime, timedelta

import math
from PIL import Image

import getopt

import warnings

import matplotlib.pyplot
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec

from mpl_toolkits.basemap import Basemap

import subprocess

from obspy.imaging.beachball import beach, mt2plane, aux_plane
from obspy.geodetics import degrees2kilometers

import matplotlib
from matplotlib.colors import LinearSegmentedColormap

matplotlib.use('Agg')
from matplotlib import ticker
import pickle

# Ignore user warnings.
warnings.filterwarnings("ignore")

# Import the aftershocks parameters and libraries.
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
lib_dir = os.path.join(parent_dir, 'lib')
param_dir = os.path.join(parent_dir, 'param')

sys.path.append(param_dir)
sys.path.append(lib_dir)

import aftershock_fdsn_maps_param as param
import aftershocks_fdsn_lib as dp_lib
import event_fdsn_lib as event_lib

"""
    Description:

    This is the Python code behind the IRIS DMC's Aftershocks Data product. It is capable of producing plots and 
    animations that are part of the IRIS DMC's Aftershocks Data product (http://ds.iris.edu/spud/aftershock).

    The code can be configured via its parameter file "aftershock_maps_param.par" or via the command line arguments. 
    Currently parameters are optimized for use with Mercator map projection, IRIS event services and IRIS GCMT data 
    product. Change in projection, resolution or data provider will require parameter tuning or code editing.

    Copyright and License:

    This software Copyright (c) 2018 IRIS (Incorporated Research Institutions for Seismology).

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or (at
    your option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/.

    History:
        2020-09-16 Manoch: V.2020.260 R2.1 public release
        2020-08-22 Manoch: V.2020.236 FDSN support
        2020-08-01 Manoch: V.2020.214 R2 release.
        2014-12-17 Alex: R1, development and initial release.

"""

# Script info.
script_version = 'V.2020.266'
script = sys.argv[0]
script = os.path.basename(script)


def event_density_range(base_map, lat_list, lon_list, bin_width_km):
    """Calculate event density and its range.
    """

    if lon_list:
        start_lon = math.floor(min(events['longitude']))
        end_lon = math.ceil(max(events['longitude']))
        start_lat = math.floor(min(events['latitude']))
        end_lat = math.ceil(max(events['latitude']))

        # form the bins
        ll_x, ll_y = base_map(start_lon, start_lat)
        ur_x, ur_y = base_map(end_lon, end_lat)
        x_km = (ur_x - ll_x) / 1000.0
        y_km = (ur_y - ll_y) / 1000.0

        # compute appropriate bins to aggregate data
        # nx is number of bins in x-axis, i.e. longitude
        # ny is number of bins in y-axis, i.e. latitude
        nx = math.ceil(x_km / bin_width_km)
        ny = math.ceil(y_km / bin_width_km)
        lon_bins = np.linspace(start_lon, end_lon, nx)
        lat_bins = np.linspace(start_lat, end_lat, ny)

        # aggregate the number of earthquakes in each bin, we will only use the density
        density, lat_edges, lon_edges = np.histogram2d(lat_list, lon_list, [lat_bins, lon_bins])
        density_min = min(map(min, density))
        density_max = max(map(max, density))

        return density, density_min, density_max

    else:
        return None, None, None


def usage():
    """The usage message.
    """
    new_line = '\n'
    print(f'{new_line}{new_line}{script} ({script_version}):')
    print(f'{new_line}This is the Python code behind the IRIS DMC\'s Aftershocks Data product. It is capable of '
          f'{new_line}producing plots and animations that are part of the IRIS DMC\'s Aftershocks Data product '
          f'{new_line}(http://ds.iris.edu/spud/aftershock).', flush=True)

    print(f'{new_line}{new_line}This is the Python code behind the IRIS DMC\'s Aftershocks Data product. '
          f'It is capable of producing plots and {new_line}animations that are part of the IRIS DMC\'s Aftershocks '
          f'Data product (http://ds.iris.edu/spud/aftershock).{new_line}{new_line}'

          f'The code can be configured via its parameter file "aftershock_maps_param.par" or via the command {new_line}'
          f'line arguments. Currently parameters are optimized for use with the Mercator map projection, NEIC/USGS FDSN'
          f'{new_line}event services and GCMT event catalog. Changes in projection, resolution or data provider may '
          f'require parameter tuning {new_line}and/or code update.'
          f'{new_line}{new_line}command line options:{new_line}'
          f'\t-h --help\t\tthis message{new_line}'
          f'\t-v --verbose\t\trun in verbose mode{new_line}'
          f'\t-e --eid [event ID]\tprocess this FDSN event ID as the mainshock{new_line}'
          f'\t-b --before [integer]\tnumber of days before the event days to display {new_line}'
          f'\t-a --after [integer]\tnumber of days after the event days to display {new_line}'
          f'\t-m --minmag [float]\tminimum event magnitude to consider{new_line}'
          f'\t-l --label [string]\tinclude this label in the file name{new_line}'
          f'\t-p --plot [1,...]\tplots to create. Use one or more (comma separated) indices from the list of possible '
          f'plots (see the "map_description" variable in the param file)." {new_line}')
    print('\nNote: Options 2 and 3 require GCMT solution.\n')
    print(f'\t\tindex\t\tdescription')
    print(f'\t\t{5 * "="}\t\t{11 * "="}')
    _desc = param.map_description
    _ind = sorted(_desc.keys())
    for _i in _ind:
        print(f'\t\t{_i}\t\t{_desc[_i]}')

    print(f'\n\t-r --refid [event ID]\tprocess this FDSN event ID as a reference event{new_line}'
          f'\t-s --scale [float] ETOPO uses scale to reduce map resolution for maps that are zoomed '
          f'in. {new_line}The default scale is automatically set between 1.2 '
          f'and 2.6. The larger scale values require more memory{new_line}'
          f'\t-x --xdate\tif -x is present, x-axis will have date labels rather than day labels{new_line}'
          f'\t-T --title [double-quoted]\tuse this title for  plots instead of the default title'
          f'{new_line}'
          f'{new_line}NOTE: either eid or refid should be provided', flush=True)
    print(f'{new_line}Example:'
          f'{new_line}{new_line}From: https://earthquake.usgs.gov/fdsnws/event/1/query?format=text&starttime='
          f'2020-07-22&endtime=2020-07-23&minmag=6.8&nodata=404'
          f'{new_line}{new_line} we obtain the event ID and then run:{new_line}'
          f'   {script} -v -e us7000asvb'
          f'{new_line} to create aftershock plots for event us7000asvb', flush=True)
    print('\n\n', flush=True)


# Get parameters from the param file.
verbose = param.verbose

basemap_resolution = param.animation_basemap_resolution
animation_basemap_boundaries = param.animation_basemap_boundaries

etopo_scale = 1.0

rate_bar_color = param.rate_bar_color
rate_bar_axes_color = param.rate_bar_axes_color
rate_bar_alpha = param.rate_bar_alpha

elevation_bar_color = param.elevation_bar_color
elevation_bar_axes_color = param.elevation_bar_axes_color
elevation_bar_alpha = param.elevation_bar_alpha

image_dir = dp_lib.mkdir(param.image_dir)
video_dir = dp_lib.mkdir(param.video_dir)
scratch_dir = dp_lib.mkdir(param.scratch_dir)
log_dir = dp_lib.mkdir(param.log_dir)
data_dir= dp_lib.mkdir(param.data_dir)

spud_service_gcmt_url = param.spud_service_gcmt_url

use_gcmt_for_mainshock = param.use_gcmt_for_mainshock

catalog_begin = param.catalog_begin

days_before = param.days_before
days_before_animation = param.days_before_animation
days_after = param.days_after
request_min_mag = param.request_min_mag

plot_label_shift = param.plot_label_shift

focal_mechanism_by_mt = param.focal_mechanism_by_mt

aftershocks_depth_radius = param.aftershocks_depth_radius

plot_label_location = param.plot_label_location

tag = param.map_tag

label_font_size = param.label_font_size
title_font_size = param.title_font_size

topo_file = param.topo_file

marker_edge_width = param.marker_edge_width

figure_size = param.figure_size

key_mag_small = param.key_mag_small
key_mag_large = param.key_mag_large

plot_unassociated_gcmts = param.plot_unassociated_gcmts

map_grid = param.map_grid

fade_duration = param.fade_duration

mainshock_beach_ball_color = param.mainshock_beach_ball_color

mag_scale = param.mag_scale

heatmap_relief_alpha = param.heatmap_relief_alpha
d_km = param.heatmap_bin_width_km

catalog_dc = param.catalog_dc
gcmt_dc = param.gcmt_dc

location_map_relief_alpha = param.location_map_relief_alpha

scale = None
_fig = None
_ax1 = None
_m = None

title_text = None
reference_mag = None

# Should xaxis be date format? (None or date)
date_xaxis = False

map_index = param.map_index
plots = param.map_list

moment_items = dict()
nodal_plane = dict()
tensor = dict()

line_break = '\n'

norm_min = None
norm_max = None
mainshock = None

gcmt_url = None
fdsn_url = None

reference_event_id = None
mainshock_id = None

em_dash = u'\u2014'

gcmt_tensor = dict()
gcmt_events = list()

log_file = sys.stdout
log_to_screen = param.log_to_screen
if not log_to_screen:
    if log_dir is None:
        log_file = open(sys.stdout, 'a')
        print(f'[WARN] Output to the log file is OFF')
    else:
        log_file_name = os.path.join(param.log_dir, script.replace('.py', ''))
        log_file_name = f"{log_file_name}_{datetime.now().strftime('%Y-%m-%d')}"
        log_file_name = '.'.join([log_file_name, 'log'])
        print(f'[INFO] Output is going to the log file {log_file_name}')
        log_file = open(log_file_name, 'a')
        sys.stdout = log_file


def set_date_xlabels(ax):
    """Get the x labels as number of days and return the corresponding dates."""
    # Get the current locations and labels.
    labels = [item.get_text() for item in ax.get_xticklabels()]
    # Set text labels and properties.
    d_label = list()
    for d in labels:
        d_label.append((mainshock['datetime'] + float(d) * 86400).strftime("%Y-%m-%d"))

    return d_label


def magnitude_vs_days(magnitude, map_tag='aftershocks'):
    """
    Plot of Magnitude vs Days after the mainshock.
    Plot beach balls when possible.
    """

    # Clean up first.
    plot_file_name = f'{mainshock_id}_{map_tag}_{map_index[1]}.png'
    plot_file_name = os.path.join(image_dir, plot_file_name)
    if os.path.exists(plot_file_name):
        os.remove(plot_file_name)

    # Set up the figure.
    fig = matplotlib.pyplot.figure(figsize=figure_size, dpi=param.dpi)

    print(f'[INFO] Starting figure magnitude_vs_days M {magnitude} {map_tag}', flush=True, file=log_file)

    # The main plot [left, bottom, width, height].
    ax1 = fig.add_axes([0.08, 0.15, 0.82, 0.6])
    ax1.spines['top'].set_visible(False)
    plt.ylabel('Magnitude')
    if not date_xaxis:
        plt.xlabel(f'Days from {event_type.lower()}')

    ax1.grid(b=True, which='major', color=map_grid['color'], alpha=map_grid['alpha'])

    # The ratio of y to x axes with 1.6 empirical corrections to keep beach balls as circles.
    ax1.set_xlim(- dp_days_before, dp_days_after)
    ax1.set_ylim(magnitude, max_mag)
    ratio = 1.6 * (ax1.get_ylim()[1] - ax1.get_ylim()[0]) / (ax1.get_xlim()[1] - ax1.get_xlim()[0])

    # Must rescale x axes due to the beach balls.
    ax1.set_xlim(ax1.get_xlim()[0] * ratio, ax1.get_xlim()[1] * ratio)

    # Production date time stamp.
    prod_date = dp_lib.version_timestamp(script_version)
    ax1.text(0.0, -0.1, prod_date, horizontalalignment='left', fontsize=5,
             verticalalignment='top', transform=ax1.transAxes)

    # The plot title goes on top.
    plt.suptitle(title, fontsize=title_font_size)

    # Plot the Large and smaller Mag circles for reference.
    circle_size = (key_mag_small['size'] / aftershock_scale_factor) ** 2.5
    plt.plot(key_mag_small['x'] + plot_label_shift['ur']['x'], key_mag_small['y'] +
             plot_label_shift['ur']['y'], 'o', markersize=circle_size, transform=ax1.transAxes,
             fillstyle='none', color='k')
    plt.text(key_mag_small['text_x'] + plot_label_shift['ur']['x'], key_mag_small['text_y'] +
             plot_label_shift['ur']['y'],
             f'M{key_mag_small["size"]}',
             fontsize=label_font_size, color='k', va='center', transform=ax1.transAxes)

    circle_size = (key_mag_large['size'] / aftershock_scale_factor) ** 2.5
    plt.plot(key_mag_large['x'] + plot_label_shift['ur']['x'], key_mag_large['y'] + plot_label_shift['ur']['y'],
             'o', markersize=circle_size, transform=ax1.transAxes,
             fillstyle='none', color='k')
    plt.text(key_mag_large['text_x'] + plot_label_shift['ur']['x'], key_mag_large['text_y'] +
             plot_label_shift['ur']['y'],
             f'M{key_mag_large["size"]}, MT\'s use GCMT\ndepth & magnitude',
             fontsize=label_font_size, color='k', va='center', transform=ax1.transAxes)

    # The main event's beach ball.
    main_color = dp_lib.get_rgba(plt, min_event_depth, max_event_depth, mainshock['depth'])
    if mainshock['mt'] is not None:
        b = beach(mainshock['mt'], xy=(0, mainshock['magnitude']),
                  width=((mainshock['magnitude'] * mag_scale[map_index[1]] / (4.6 * aftershock_scale_factor)) ** 2.5),
                  linewidth=0.5, alpha=1,
                  facecolor=main_color, edgecolor='k', axes=ax1)
        b.set_zorder(10000)
        ax1.add_collection(b)
    else:
        size = dp_lib.event_marker_size(mainshock['magnitude'], aftershock_scale_factor)
        plt.plot(0, mainshock['magnitude'], 'o', markersize=size, zorder=10000,
                 markerfacecolor=main_color, markeredgecolor='k',
                 markeredgewidth=marker_edge_width)

    # Plot the aftershocks.
    hourly_rate = dict()

    # Other events.
    hourly_rate['0.0'] = 1
    warn_issued = False
    for _index, _elon in enumerate(events['longitude']):

        hour = f'{math.floor(events["time_diff_sec"][_index] / 3600.0):.1f}'

        if events['magnitude'][_index] >= min_mag:
            if hour not in hourly_rate.keys():
                hourly_rate[hour] = 1
            else:
                hourly_rate[hour] += 1

        size = dp_lib.event_marker_size(events['magnitude'][_index], aftershock_scale_factor)
        x = ratio * events['time_diff_sec'][_index] / 86400
        y = events['magnitude'][_index]

        # GCMT event.
        not_gcmt = True
        if events['mt'][_index] is not None:
            try:
                b = beach(events['mt'][_index],
                          xy=(x, y),
                          width=((events['magnitude'][_index] *
                                  mag_scale[map_index[1]] / (4.6 * aftershock_scale_factor)) ** 2.5),
                          linewidth=0.5, alpha=1,
                          facecolor=dp_lib.get_rgba(plt, min_event_depth, max_event_depth, events['depth'][_index]),
                          edgecolor='k', axes=ax1)
                b.set_zorder(5000)
                ax1.add_collection(b)
                not_gcmt = False
            except Exception as ex:
                if not warn_issued:
                    print(f'[WARN] GCMP beach ball failed for event {_index}')
                    warn_issued = True
                not_gcmt = True
        # Regular events.
        if not_gcmt:
            plt.plot(x, y, 'o', markersize=size,
                     markerfacecolor=dp_lib.get_rgba(plt, min_event_depth, max_event_depth, events['depth'][_index]),
                     markeredgecolor='k', zorder=500,
                     markeredgewidth=marker_edge_width)

    # Fudge to fix the X labeling
    ax1.get_xaxis().set_ticks([])
    xtick_locations, xtick_labels = dp_lib.days_ticks(ratio, dp_days_before, dp_days_after,
                                                      plot_days_tick_interval)
    plt.xticks(np.array(xtick_locations), np.array(xtick_labels))

    # Should the axis have date label?
    if date_xaxis:
        ax1.set_xticklabels(set_date_xlabels(ax1), fontsize=6, ha='right', rotation=20)

    # Plot the colorbar.
    ax = fig.add_axes([0.915, 0.15, 0.01, 0.6])
    cmap, norm = dp_lib.get_rgba(plt, min_event_depth, max_event_depth, None)
    cb = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='vertical')
    cb.ax.invert_yaxis()
    tick_locator = ticker.MaxNLocator(nbins=5)
    cb.locator = tick_locator
    cb.update_ticks()
    cb.set_label("depth (km)")

    # Shade times before the mainshock.
    ax1.axvspan(2.0 * ax1.get_xlim()[0], 0.0, alpha=0.05, color='k')

    # Plot the hourly rate bar chart using the same rati that we used for the main plot.
    hours_large = list()
    values_large = list()
    for key in sorted(hourly_rate):
        this_day = float(key) / 24.0
        if - dp_days_before <= this_day <= dp_days_after:
            hours_large.append(ratio * this_day)
            values_large.append(hourly_rate[key])

    # The hourly rate plot [left, bottom, width, height].
    ax1t = fig.add_axes([0.08, 0.75, 0.82, 0.15], sharex=ax1)
    ax1t.spines['right'].set_visible(False)
    ax1t.spines['top'].set_visible(False)
    ax1t.spines['left'].set_visible(False)
    ax1t.spines['bottom'].set_visible(False)
    plt.setp(ax1t.get_xticklabels(), visible=False)

    ax1t.bar(hours_large, values_large, color=rate_bar_color, width=0.05, alpha=rate_bar_alpha)

    # Various labels and legends.
    _text = f'M ≥ {magnitude} events {date_range_label} within {int(max_radius_km)} km radius and ' \
        f'within {aftershocks_depth_radius} km depth'

    ax1t.set_title(_text, fontsize=label_font_size)
    ax1t.set_ylabel(f'Event rate\n(per hour)', fontsize=label_font_size, color=rate_bar_axes_color)
    ax1t.tick_params(axis='y', colors=rate_bar_axes_color, labelsize=label_font_size, labelcolor=rate_bar_axes_color)

    # Plot the logo.
    if os.path.isfile(param.logo_image):
        im = Image.open(param.logo_image)
        im.thumbnail((param.logo_width, param.logo_height), Image.ANTIALIAS)  # resizes image in-place
        fig.figimage(im, param.logo_x, param.logo_y)

    # Save the figure to a .png file
    matplotlib.pyplot.savefig(plot_file_name)

    plt.close('all')
    print(f'[INFO] Done with figure {map_tag}: {plot_file_name} ', flush=True, file=log_file)


def strike_vs_day(distance_along_strike, map_tag='aftershocks'):
    """

     Group 2: Make figure of Dist along Strikes vs Days from mainshock.
     Plot beac hballs when possible.

    """
    warn_issued = False
    # Two strike angles are possible.
    for strike_index, dist in enumerate(distance_along_strike):
        # Clean up first!
        plot_file_name = f'{mainshock_id}_{map_tag}_{map_index[2]}_{strike_index + 1}.png'
        plot_file_name = os.path.join(image_dir, plot_file_name)
        if os.path.exists(plot_file_name):
            os.remove(plot_file_name)

        # Set up the figure.
        fig = matplotlib.pyplot.figure(figsize=figure_size, dpi=param.dpi)

        # Set up the main subplot (ax1).
        gs1 = gridspec.GridSpec(3, 3)
        gs1.update(left=0.12, right=0.9, bottom=0.15, top=0.9, wspace=0.01)
        ax1 = plt.subplot(gs1[:, :])

        # Production date time stamp.
        production_date = dp_lib.version_timestamp(script_version)
        ax1.text(0.0, -0.1, production_date, horizontalalignment='left', fontsize=5,
                 verticalalignment='top', transform=ax1.transAxes)

        # Axis labels, limits and grids.
        plt.ylabel('Distance Along Strike (km)')

        if not date_xaxis:
            plt.xlabel(f'Days from {event_type.lower()}')

        ax1.set_xlim(-1 * dp_days_before, dp_days_after, plot_days_tick_interval)
        ax1.set_ylim(1.02 * min(distance_along_strike_km), 1.02 * max(distance_along_strike_km))
        ax1.grid(b=True, which='major', color=map_grid['color'], alpha=map_grid['alpha'])

        # Ratio of y to x axes with 1.6 empirical corrections to keep beach balls as circles.
        ratio = 1.6 * (ax1.get_ylim()[1] - ax1.get_ylim()[0]) / (ax1.get_xlim()[1] - ax1.get_xlim()[0])

        # Plot the Large and smaller Mag circles for reference.
        this_size = (key_mag_small['size'] / aftershock_scale_factor) ** 2.5
        plt.plot(key_mag_small['x'], key_mag_small['y'], 'o', markersize=this_size, transform=ax1.transAxes,
                 fillstyle='none', color='k', clip_on=False)
        plt.text(key_mag_small['text_x'], key_mag_small['text_y'],
                 f'M{key_mag_small["size"]}',
                 fontsize=label_font_size, color='k', va='center', transform=ax1.transAxes, clip_on=False)

        this_size = (key_mag_large['size'] / aftershock_scale_factor) ** 2.5
        plt.plot(key_mag_large['x'], key_mag_large['y'], 'o', markersize=this_size, transform=ax1.transAxes,
                 fillstyle='none', color='k', clip_on=False)
        plt.text(key_mag_large['text_x'], key_mag_large['text_y'],
                 f'M{key_mag_large["size"]}, MT\'s use GCMT\nlocation, depth & magnitude',
                 fontsize=label_font_size, color='k', va='center', transform=ax1.transAxes, clip_on=False)

        # Main event's beach ball.
        if mainshock['mt'] is not None:
            b = beach(mainshock['mt'], xy=(0.0, 0.0),
                      width=((mainshock['magnitude'] * mag_scale[map_index[2]] /
                              (4.8 * aftershock_scale_factor)) ** 2.5),
                      linewidth=0.5, alpha=1,
                      facecolor=mainshock_color, edgecolor='k', axes=ax1)
            b.set_zorder(10000)
            ax1.add_collection(b)
        else:
            size = dp_lib.event_marker_size(mainshock['magnitude'], aftershock_scale_factor)
            plt.plot(0.0, 0.0, 'o', markersize=size, zorder=10000,
                     markerfacecolor=mainshock_color, markeredgecolor='k',
                     markeredgewidth=marker_edge_width)

        # Plot the aftershocks.
        for lon_index, elon in enumerate(events['longitude']):
            size = dp_lib.event_marker_size(events['magnitude'][lon_index], aftershock_scale_factor)

            # GCMT event.
            not_gcmt = True
            if events['mt'][lon_index] is not None:
                try:
                    b = beach(events['mt'][lon_index], xy=(ratio * events['time_diff_sec'][lon_index] / 86400.0,
                                                        events['distance_along_strike'][strike_index][lon_index]),
                              width=((events['magnitude'][lon_index] * mag_scale[map_index[2]]
                                      / (4.6 * aftershock_scale_factor)) ** 2.5),
                              linewidth=0.5, alpha=1,
                              facecolor=aftershock_color[lon_index], edgecolor='k', axes=ax1)
                    b.set_zorder(5000)
                    ax1.add_collection(b)
                    not_gcmt = False
                except:
                    if not warn_issued:
                        print(f'[WARN] GCMP beach ball failed for event {lon_index}')
                        warn_issued = True
                    not_gcmt = True
            # Regular events.
            if not_gcmt:
                plt.plot(ratio * events['time_diff_sec'][lon_index] / 86400.0,
                         events['distance_along_strike'][strike_index][lon_index], 'o',
                         markersize=size, markeredgecolor='k', markeredgewidth=marker_edge_width,
                         markerfacecolor=aftershock_color[lon_index], zorder=500)

        # Set the plot title.
        _text = f'M ≥ {min_mag} events {date_range_label} within {int(max_radius_km)} km radius and ' \
            f'within {aftershocks_depth_radius} km depth along strike {int(azimuth[strike_index])} deg'
        plt.suptitle(title, fontsize=title_font_size)
        ax1.set_title(_text, fontsize=label_font_size)

        # Plot directions, N/S/E/W on the plot.
        coord1, coord2 = dp_lib.get_directions(azimuth[strike_index])
        plt.text(0.01, 0.97, coord1,
                 fontsize=label_font_size,
                 color='k', transform=ax1.transAxes)
        plt.text(0.01, 0.01, coord2,
                 fontsize=label_font_size,
                 color='k', transform=ax1.transAxes)

        # Fudge to fix the Xaxis labeling
        ax1.get_xaxis().set_ticks([])
        xtick_locations, xtick_labels = dp_lib.days_ticks(ratio, dp_days_before, dp_days_after,
                                                          plot_days_tick_interval)
        plt.xticks(np.array(xtick_locations), np.array(xtick_labels))

        # Should the axis have date label?
        if date_xaxis:
            ax1.set_xticklabels(set_date_xlabels(ax1), fontsize=6, ha='right', rotation=20)

        # Plot the colorbar.
        ax = fig.add_axes([0.915, 0.15, 0.01, 0.75])
        cmap, norm = dp_lib.get_rgba(plt, min_event_depth, max_event_depth, None)
        cb = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='vertical')
        cb.ax.invert_yaxis()
        tick_locator = ticker.MaxNLocator(nbins=5)
        cb.locator = tick_locator
        cb.update_ticks()

        cb.set_label("depth (km)")

        # Shade times before the mainshock.
        ax1.axvspan(2.0 * ax1.get_xlim()[0], 0.0, alpha=0.05, color='k')

        # Plot the logo.
        if os.path.isfile(param.logo_image):
            im = Image.open(param.logo_image)
            im.thumbnail((param.logo_width, param.logo_height), Image.ANTIALIAS)  # resizes image in-place
            fig.figimage(im, param.logo_x, param.logo_y)

        # Print the figure to a .png file
        matplotlib.pyplot.savefig(plot_file_name)

        plt.close('all')

        print(f'[INFO] Done with color figure {map_tag}: {plot_file_name} ', flush=True, file=log_file)


def strike_vs_depth(distance_along_strike, map_tag='aftershocks'):
    """

     Group 3: Make figure of Dist along Strikes vs depth .
     Plot beachballs when possible.

    """
    # Two strike angles are possible.
    for strike_index, dist in enumerate(distance_along_strike):
        # Clean up first.
        plot_file_name = f'{mainshock_id}_{map_tag}_{map_index[3]}_{strike_index + 1}.png'
        plot_file_name = os.path.join(image_dir, plot_file_name)
        if os.path.exists(plot_file_name):
            os.remove(plot_file_name)

        # Here we are saying that we should not do small magnitudes (index 0).

        # Set up the figure.
        fig = matplotlib.pyplot.figure(figsize=figure_size, dpi=param.dpi)

        # The main plot.
        # Set up the axes [left, bottom, width, height]
        ax1 = fig.add_axes([0.08, 0.15, 0.82, 0.6], clip_on=False)
        ax1.spines['top'].set_visible(False)
        plt.ylabel('Depth (km)')
        plt.xlabel('Distance Along Strike (km)')
        plt.gca().invert_yaxis()
        ax1.grid(b=True, which='major', color=map_grid['color'], alpha=map_grid['alpha'])

        # Ratio of y to x axes with 1.6 empirical corrections to keep beach balls as circles (we do not plot beach balls
        # here but keep the code consistent with the rest).
        ax1.set_ylim(max_event_depth * 1.2, 0)
        ax1.set_xlim(1.02 * min(distance_along_strike_km), 1.02 * max(distance_along_strike_km))

        # Production date time stamp.
        production_date = dp_lib.version_timestamp(script_version)
        ax1.text(0.0, -0.1, production_date, horizontalalignment='left', fontsize=5,
                 verticalalignment='top', transform=ax1.transAxes)

        # The plot title goes on top.
        plt.suptitle(title, fontsize=title_font_size)

        # Plot the Large and smaller Mag circles for reference.
        this_size = (key_mag_small['size'] / aftershock_scale_factor) ** 2.5
        plt.plot(key_mag_small['x'], key_mag_small['y'], 'o', markersize=this_size, transform=ax1.transAxes,
                 fillstyle='none', color='k')
        plt.text(key_mag_small['text_x'], key_mag_small['text_y'],
                 f'M{key_mag_small["size"]}',
                 fontsize=label_font_size, color='k', va='center', transform=ax1.transAxes)

        this_size = (key_mag_large['size'] / aftershock_scale_factor) ** 2.5
        plt.plot(key_mag_large['x'], key_mag_large['y'], 'o', markersize=this_size, transform=ax1.transAxes,
                 fillstyle='none', color='k')
        plt.text(key_mag_large['text_x'], key_mag_large['text_y'],
                 f'M{key_mag_large["size"]}',  # ', MT\'s use GCMT\ndepth & magnitude',
                 fontsize=label_font_size, color='k', va='center', transform=ax1.transAxes)

        # The main event.

        size = dp_lib.event_marker_size(mainshock['magnitude'], aftershock_scale_factor)
        x = 0
        y = mainshock['depth']
        plt.plot(x, y, 'o',
                 markersize=size, markeredgecolor='k', markeredgewidth=marker_edge_width,
                 markerfacecolor=mainshock_day_color, )

        # Plot the aftershocks.
        for _index, _elon in enumerate(events['longitude']):
            size = dp_lib.event_marker_size(events['magnitude'][_index], aftershock_scale_factor)
            x = events['distance_along_strike'][strike_index][_index]
            y = events['depth'][_index]
            depth = float(events['depth'][_index])
            plt.plot(x, y, 'o',
                     markersize=size, markeredgecolor='k', markeredgewidth=marker_edge_width,
                     markerfacecolor=aftershock_time_color[_index])

        plt.setp(ax1.get_yticklabels()[0], visible=False)

        # Set the plot title.
        plt.suptitle(title, fontsize=title_font_size)
        plot_text = f'Other events within {int(max_radius_km)} km radius; {date_range_label}; ' \
            f'within {aftershocks_depth_radius} km depth; ' \
            f'M ≥ {min_mag}; strike {int(azimuth[strike_index])} deg'

        # Set up the axes and for each one set [left, bottom, width, height]
        # The elevation plot.
        ax1t = fig.add_axes([0.08, 0.75, 0.82, 0.15], sharex=ax1)
        ax1t.spines['right'].set_visible(False)
        ax1t.spines['top'].set_visible(False)
        ax1t.spines['left'].set_visible(False)
        ax1t.spines['bottom'].set_visible(False)
        plt.setp(ax1t.get_xticklabels(), visible=False)
        ax1t.set_title(plot_text, fontsize=label_font_size)

        ax1t.plot(distance_along_strike_km, elevation[strike_index], '-', alpha=elevation_bar_alpha,
                  color=elevation_bar_axes_color, lw=0.5)
        ax1t.fill_between(distance_along_strike_km, elevation[strike_index], ax1t.get_ylim()[0],
                          alpha=elevation_bar_alpha,
                          color=elevation_bar_color)
        ax1t.set_ylabel(f'{param.topo_file_name}\nElevation (m)', fontsize=6, color=elevation_bar_axes_color)
        ax1t.tick_params(axis='y', colors=elevation_bar_axes_color, labelsize=6, labelcolor=elevation_bar_axes_color)
        ax1t.hlines(0, min(distance_along_strike_km), max(distance_along_strike_km),
                    colors='k', linestyles='dashed', lw=0.5, color=elevation_bar_axes_color)
        ax1t.text(min(distance_along_strike_km) + 0.1, 0.1, 'sea level', color=elevation_bar_axes_color,
                  fontsize=label_font_size)

        # Plot the N/S/E/W on the plot.
        coord1, coord2 = dp_lib.get_directions(azimuth[strike_index])
        ax1t.text(0.97, 0.75, coord1,
                  fontsize=label_font_size,
                  color=elevation_bar_axes_color, transform=ax1t.transAxes, clip_on=False)
        ax1t.text(0.01, 0.75, coord2,
                  fontsize=label_font_size,
                  color=elevation_bar_axes_color, transform=ax1t.transAxes, clip_on=False)

        # Plot the colorbar.
        # Set up the axes [left, bottom, width, height] (set the same as ax1)
        ax = fig.add_axes([0.915, 0.15, 0.01, .6])
        cmap, norm = dp_lib.get_rgba(plt, - days_before, days_after, None)
        cb = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='vertical')
        # cb.ax.invert_yaxis()
        tick_locator = ticker.MaxNLocator(nbins=5)
        cb.locator = tick_locator
        cb.update_ticks()
        cb.set_label(f"days from {event_type.lower()}")

        # Plot the logo.
        if os.path.isfile(param.logo_image):
            im = Image.open(param.logo_image)
            im.thumbnail((param.logo_width, param.logo_height), Image.ANTIALIAS)  # resizes image in-place
            fig.figimage(im, param.logo_x, param.logo_y)

        # Print the figure to a .png file
        matplotlib.pyplot.savefig(plot_file_name)

        plt.close('all')

        print(f'[INFO] Done with figure {map_tag}: {plot_file_name} ', flush=True, file=log_file)


def location_map(magnitude, map_tag='aftershocks', map_type='location_map', q_start=None, q_end=None):
    """

     Group 4: Location map of aftershocks.

    """
    # Clean up first!
    plot_file_name = f'{mainshock_id}_{map_tag}_{map_type}.png'
    plot_file_name = os.path.join(image_dir, plot_file_name)
    if os.path.exists(plot_file_name):
        os.remove(plot_file_name)

    # Some map parameters.
    map_size_lon = map_size_lat / math.cos(mainshock['latitude'] * math.pi / 180)
    lat_bounds = mainshock['latitude'] - map_size_lat, mainshock['latitude'] + map_size_lat
    lon_bounds = mainshock['longitude'] - map_size_lon, mainshock['longitude'] + map_size_lon

    # Generate and save  earthquake browser link for the event.
    if None not in (q_start, q_end):
        browser_start = q_start.split('T')[0]
        browser_end = q_end.split('T')[0]
        eq_browser_seimicity_url, eq_browser_aftershocks_url = dp_lib.get_eq_browser_url(browser_start,
                                                                                         browser_end,
                                                                                         lat_bounds[0], lat_bounds[1],
                                                                                         lon_bounds[0],
                                                                                         lon_bounds[1], min_mag)
        fp = open(os.path.join(data_dir,
                               f'{mainshock_id}_events_and_iris_earthquake_browser_links.html'), 'w')

        _url = fdsn_url
        _url = _url.replace('format=xml', 'format=text')
        _url = _url.replace('&includeallorigins=true', '')
        _url = _url.replace('&includeallmagnitudes=true', '')

        fp.write(f'<h3>{title.replace("°", "&deg;")}</h3>'
                 f'<br /><br />'
                 f'<a href="{eq_browser_aftershocks_url}" target="_NEW">'
                 f'IRIS Earthquake Browser URL (M &ge; {magnitude} from {browser_start} to {browser_end})</a>'
                 f'<br /><br />'
                 f'<a href="{eq_browser_seimicity_url}" target="_NEW">'
                 f'IRIS Earthquake Browser URL (M &ge; {magnitude})</a>'
                 f'<br /><br />'
                 f'<a href="{_url}" target="_NEW">'
                 f'FDSN Event Service URL (M &ge; {magnitude} from {browser_start} to {browser_end})</a>'
                 f'<br /><br />'
                 f'<a href="{gcmt_url}" target="_NEW">'
                 f'GCMT catalog URL (M &ge; {magnitude} from {browser_start} to {browser_end})</a>'
                 )
        fp.close()

    # Set up the figure.
    fig = matplotlib.pyplot.figure(figsize=(figure_size[0], figure_size[1]), dpi=param.dpi)

    gs1 = gridspec.GridSpec(3, 3)
    gs1.update(left=0.05, right=0.9, bottom=0.05, top=0.9, wspace=0.01)
    ax1 = plt.subplot(gs1[:, :])

    # Production date time stamp.
    production_date = dp_lib.version_timestamp(script_version)
    ax1.text(0.01, 0.02, production_date, horizontalalignment='left', fontsize=5,
             verticalalignment='top', transform=ax1.transAxes)

    # The basemap.
    m = Basemap(projection='merc', llcrnrlat=lat_bounds[0], urcrnrlat=lat_bounds[1], llcrnrlon=lon_bounds[0],
                urcrnrlon=lon_bounds[1], lat_ts=mainshock['latitude'],
                resolution=basemap_resolution)  # lat_ts = 20, lat of true scale
    # m.shadedrelief(alpha=location_map_relief_alpha, zorder=10)
    try:
        m.drawcoastlines(color='gray')
    except Exception as ex:
        print(f'[ERR] draw error {ex}')
        pass

    m.drawcountries(color='gray')
    m.drawstates(color='gray')
    m.etopo(scale=etopo_scale, alpha=location_map_relief_alpha, zorder=10)

    m.drawparallels(np.arange(-90., 91., lat_line_spacing), labels=[1, 0, 0, 0], fontsize=6, linewidth=0.3, alpha=0.5)
    m.drawmeridians(np.arange(-180., 181., lat_line_spacing), labels=[0, 0, 0, 1], fontsize=6, linewidth=0.3, alpha=0.5)
    # Set the plot title.
    plt.suptitle(title, fontsize=title_font_size)

    # Ratio of y to x axes with 1.6 empirical corrections to keep beach balls as circles.
    ratio = 1.6 * (ax1.get_ylim()[1] - ax1.get_ylim()[0]) / (ax1.get_xlim()[1] - ax1.get_xlim()[0])

    # Plot label on the upper right.
    plot_text = f'Other events within {int(max_radius_km)} km radius; {date_range_label}; ' \
        f'within {aftershocks_depth_radius} km depth; M ≥ {magnitude}'
    ax1.set_title(plot_text, fontsize=label_font_size)

    # Plot M7 & M5 dots for reference.
    this_size = (key_mag_small['size'] / aftershock_scale_factor) ** 2.5
    plt.plot(key_mag_small['x'], key_mag_small['y'], 'o', markersize=this_size,
             transform=ax1.transAxes,
             fillstyle='none', color='k')
    plt.text(key_mag_small['text_x'], key_mag_small['text_y'],
             f'M{key_mag_small["size"]}',
             fontsize=label_font_size, color='k', va='center', transform=ax1.transAxes)

    this_size = (key_mag_large['size'] / aftershock_scale_factor) ** 2.5
    plt.plot(key_mag_large['x'], key_mag_large['y'], 'o', markersize=this_size,
             transform=ax1.transAxes,
             fillstyle='none', color='k')
    plt.text(key_mag_large['text_x'], key_mag_large['text_y'],
             f'M{key_mag_large["size"]}, MT\'s use GCMT\nlocation, depth & magnitude',
             fontsize=label_font_size, color='k', va='center', transform=ax1.transAxes)

    # Plot the aftershocks.
    warn_issued = False
    for lon_index, elon in enumerate(events['longitude']):
        xmap, ymap = m(events['longitude'][lon_index], events['latitude'][lon_index])
        size = (events['magnitude'][lon_index] / aftershock_scale_factor) ** 2.5

        # GCMT event.
        not_gcmt = True
        try:
            if events['mt'][lon_index] is not None:
                not_gcmt = False
                b = beach(events['mt'][lon_index], xy=(xmap, ymap),
                          width=((events['magnitude'][lon_index] * mag_scale['location_map']
                                  / (4.6 * aftershock_scale_factor)) ** 2.5),
                          linewidth=0.5, alpha=1,
                          facecolor=aftershock_color[lon_index], edgecolor='k', axes=ax1)
                b.set_zorder(5000)
                ax1.add_collection(b)
        except:
            if not warn_issued:
                print(f'[WARN] GCMP beach ball failed for event {lon_index}')
                warn_issued = True
            not_gcmt = True

        if not_gcmt:
            m.plot(xmap, ymap, markerfacecolor=aftershock_color[lon_index], markeredgecolor='k',
                   markeredgewidth=marker_edge_width,
                   marker='o', zorder=500,
                   markersize=size)

    # -------------- plot MainShock GCMT (slightly transparent)
    if mainshock['mt'] is not None:
        b = beach(mainshock['mt'], xy=m(mainshock['longitude'], mainshock['latitude']),
                  width=map_size_lat * ((aftershock_scale_factor * 3.3 * mainshock['magnitude']) ** 2.5), linewidth=0.5,
                  alpha=0.3,
                  facecolor=mainshock_color, edgecolor='k')
        b.set_zorder(10000)
        ax1.add_collection(b)
    else:
        xmap, ymap = m(mainshock['longitude'], mainshock['latitude'])
        size = (mainshock['magnitude'] / aftershock_scale_factor) ** 2.5
        m.plot(xmap, ymap, marker='*', markerfacecolor=mainshock_color, markeredgecolor='k', markersize=size,
               markeredgewidth=marker_edge_width, zorder=10000, alpha=0.6)

    # Draw the search area circle.
    x0, y0 = m(mainshock['longitude'], mainshock['latitude'])
    x1, y1 = m(mainshock['longitude'], mainshock['latitude'] + max_radius)
    radius_map = y1 - y0
    circle = plt.Circle((x0, y0), radius_map, color=elevation_bar_color, alpha=elevation_bar_alpha,
                        lw=2, fill=False)
    ax1.add_patch(circle)

    xmap, ymap = m(mainshock['longitude'], mainshock['latitude'])
    size = (mainshock['magnitude'] / aftershock_scale_factor) ** 2.5
    if mainshock['mt'] is None:
        m.plot(xmap, ymap, marker='*', markerfacecolor='y', markeredgecolor='k', markersize=size, alpha=0.6)

    # ---- if using projection = 'merc'
    scale_bar_lon = lon_bounds[0] + 1.6 * map_size_lon
    scale_bar_lat = lat_bounds[0] + map_size_lat * 0.15
    m.drawmapscale(scale_bar_lon, scale_bar_lat, 0, 0, scalebar_length_km, units='km',
                   labelstyle='simple',
                   fillcolor1='w', fillcolor2=elevation_bar_axes_color,
                   fontcolor='#555555',
                   zorder=5, fontsize=label_font_size)
    # Plot the colorbar.
    # on the figure total in percent [left, bottom, width, height]
    ax = fig.add_axes([0.82, 0.2, 0.01, 0.55])
    cmap, norm = dp_lib.get_rgba(plt, min_event_depth, max_event_depth, None)
    cb = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='vertical')
    cb.ax.invert_yaxis()
    tick_locator = ticker.MaxNLocator(nbins=5)
    cb.locator = tick_locator
    cb.update_ticks()
    cb.set_label("depth (km)")

    # Plot the logo.
    if os.path.isfile(param.logo_image):
        im = Image.open(param.logo_image)
        im.thumbnail((param.logo_width, param.logo_height), Image.ANTIALIAS)  # resizes image in-place
        fig.figimage(im, 0.79 * param.logo_x, param.logo_y, zorder=1000)

    # fig.tight_layout()

    # Print the figure to a .png file
    matplotlib.pyplot.savefig(plot_file_name, bbox_inches='tight')

    plt.close('all')

    print(f'[INFO] Done with figure {map_tag}: {plot_file_name} ', flush=True, file=log_file)


def heatmap(magnitude, map_tag='aftershocks'):
    """

     Group 5: Heatmap of aftershocks.

    """
    # Clean up first!
    plot_file_name = f'{mainshock_id}_{map_tag}_{map_index[5]}.png'
    plot_file_name = os.path.join(image_dir, plot_file_name)
    if os.path.exists(plot_file_name):
        os.remove(plot_file_name)

    # Some map parameters.
    map_size_lon = map_size_lat / math.cos(mainshock['latitude'] * math.pi / 180)
    lat_bounds = mainshock['latitude'] - map_size_lat, mainshock['latitude'] + map_size_lat
    lon_bounds = mainshock['longitude'] - map_size_lon, mainshock['longitude'] + map_size_lon

    # Set up the figure.
    fig = matplotlib.pyplot.figure(figsize=(figure_size[0], figure_size[1]), dpi=param.dpi)

    gs1 = gridspec.GridSpec(3, 3)
    gs1.update(left=0.05, right=0.9, bottom=0.05, top=0.9, wspace=0.01)
    ax1 = plt.subplot(gs1[:, :])

    # Production date time stamp.
    production_date = dp_lib.version_timestamp(script_version)
    ax1.text(0.01, 0.02, production_date, horizontalalignment='left', fontsize=5,
             verticalalignment='top', transform=ax1.transAxes)

    # The basemap.
    m = Basemap(projection='merc', llcrnrlat=lat_bounds[0], urcrnrlat=lat_bounds[1], llcrnrlon=lon_bounds[0],
                urcrnrlon=lon_bounds[1], lat_ts=mainshock['latitude'],
                resolution=basemap_resolution)  # lat_ts = 20, lat of true scale

    # m.shadedrelief(alpha=location_map_relief_alpha, zorder=10)
    try:
        m.drawcoastlines(color='gray')
    except Exception as ex:
        print(f'[ERR] draw error {ex}')
        pass

    m.drawcountries(color='gray')
    m.drawstates(color='gray')
    m.etopo(scale=etopo_scale, alpha=heatmap_relief_alpha, zorder=10)

    m.drawparallels(np.arange(-90., 91., lat_line_spacing), labels=[1, 0, 0, 0], fontsize=6, linewidth=0.3, alpha=0.5)
    m.drawmeridians(np.arange(-180., 181., lat_line_spacing), labels=[0, 0, 0, 1], fontsize=6, linewidth=0.3, alpha=0.5)

    # Set the plot title.
    plt.suptitle(title, fontsize=title_font_size)

    # Ratio of y to x axes with 1.6 empirical corrections to keep beach balls as circles.
    ratio = 1.6 * (ax1.get_ylim()[1] - ax1.get_ylim()[0]) / (ax1.get_xlim()[1] - ax1.get_xlim()[0])

    # Plot label on the upper right.
    plot_text = f'Other events within {int(max_radius_km)} km radius; {date_range_label}; ' \
        f'within {aftershocks_depth_radius} km depth; ' \
        f'M ≥ {magnitude}'
    ax1.set_title(plot_text, fontsize=label_font_size)

    # Plot M7 & M5 dots for reference.
    this_size = (key_mag_small['size'] / aftershock_scale_factor) ** 2.5
    plt.plot(key_mag_small['x'], key_mag_small['y'], 'o', markersize=this_size,
             transform=ax1.transAxes,
             fillstyle='none', color='k')
    plt.text(key_mag_small['text_x'], key_mag_small['text_y'],
             f'M{key_mag_small["size"]}',
             fontsize=label_font_size, color='k', va='center', transform=ax1.transAxes)

    this_size = (key_mag_large['size'] / aftershock_scale_factor) ** 2.5
    plt.plot(key_mag_large['x'], key_mag_large['y'], 'o', markersize=this_size,
             transform=ax1.transAxes,
             fillstyle='none', color='k')
    plt.text(key_mag_large['text_x'], key_mag_large['text_y'],
             f'M{key_mag_large["size"]}, MT\'s use GCMT\nlocation, depth & magnitude',
             fontsize=label_font_size, color='k', va='center', transform=ax1.transAxes)

    # Plot heatmap on Basemap.

    lon_list = list()
    lon_list.append(mainshock['longitude'])
    lat_list = list()
    lat_list.append(mainshock['latitude'])
    for _index, _elon in enumerate(events['longitude']):
        lon_list.append(_elon)
        lat_list.append(events['latitude'][_index])

    start_lon = math.floor(min(lon_list))
    end_lon = math.ceil(max(lon_list))
    start_lat = math.floor(min(lat_list))
    end_lat = math.ceil(max(lat_list))

    # compute appropriate bins to aggregate data
    # nx is number of bins in x-axis, i.e. longitude
    # ny is number of bins in y-axis, i.e. latitude
    # form the bins
    ll_x, ll_y = m(start_lon, start_lat)
    ur_x, ur_y = m(end_lon, end_lat)
    x_km = (ur_x - ll_x) / 1000.0
    y_km = (ur_y - ll_y) / 1000.0

    nx = math.ceil(x_km / d_km)
    ny = math.ceil(y_km / d_km)
    lon_bins = np.linspace(start_lon, end_lon, nx)
    lat_bins = np.linspace(start_lat, end_lat, ny)

    # aggregate the number of earthquakes in each bin, we will only use the density
    density, lat_edges, lon_edges = np.histogram2d(lat_list, lon_list, [lat_bins, lon_bins])

    # get the mesh for the lat and lon
    lon_bins_2d, lat_bins_2d = np.meshgrid(lon_bins, lat_bins)
    # convert the bin mesh to map coordinates:
    xs, ys = m(lon_bins_2d, lat_bins_2d)  # will be plotted using pcolormesh

    cmap = plt.get_cmap(param.heatmap_color_map)

    # Here adding one row and column at the end of the matrix, so that
    # density has same dimension as xs, ys, otherwise, using shading='gouraud'
    # will raise error
    density = np.hstack((density, np.zeros((density.shape[0], 1))))
    density = np.vstack((density, np.zeros((density.shape[1]))))

    # Plot heatmap with the color map
    m.pcolormesh(xs, ys, density, cmap=cmap, shading='gouraud')

    # -------------- plot yellow star MainShock GCMT (slightly transparent)

    if mainshock['mt'] is not None:
        if mainshock['mt'][0] != 0 and mainshock['mt'][1] != 0:
            b = beach(mainshock['mt'], xy=m(mainshock['longitude'], mainshock['latitude']),
                      width=map_size_lat * ((aftershock_scale_factor * mag_scale['heatmap']
                                             * mainshock['magnitude']) ** 2.5), linewidth=0.5,
                      alpha=0.4,
                      facecolor=mainshock_color, edgecolor='k')
            b.set_zorder(10000)
            ax1.add_collection(b)
    else:
        xmap, ymap = m(mainshock['longitude'], mainshock['latitude'])
        size = (mainshock['magnitude'] / aftershock_scale_factor) ** 2.5
        m.plot(xmap, ymap, marker='*', markerfacecolor=mainshock_color, markeredgecolor='k', markersize=size,
               markeredgewidth=marker_edge_width, zorder=10000, alpha=0.6)

    # Draw the search area circle.
    x0, y0 = m(mainshock['longitude'], mainshock['latitude'])
    x1, y1 = m(mainshock['longitude'], mainshock['latitude'] + max_radius)
    radius_map = y1 - y0
    circle = plt.Circle((x0, y0), radius_map, color=elevation_bar_color, alpha=elevation_bar_alpha,
                        lw=2, fill=False)
    ax1.add_patch(circle)

    # ---- if using projection = 'merc'
    scale_bar_lon = lon_bounds[0] + 1.6 * map_size_lon
    scale_bar_lat = lat_bounds[0] + map_size_lat * 0.15
    m.drawmapscale(scale_bar_lon, scale_bar_lat, 0, 0, scalebar_length_km, units='km',
                   labelstyle='simple',
                   fillcolor1='w', fillcolor2=elevation_bar_axes_color,
                   fontcolor='#555555',
                   zorder=5, fontsize=label_font_size)

    # Plot the colorbar.
    # on the figure total in percent [left, bottom, width, height]
    ax = fig.add_axes([0.82, 0.2, 0.01, 0.55])
    norm_min = min(map(min, density))
    norm_max = max(map(max, density))
    norm = matplotlib.colors.Normalize(vmin=norm_min, vmax=norm_max)
    cb = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='vertical')
    tick_locator = ticker.MaxNLocator(nbins=5)
    cb.locator = tick_locator
    cb.update_ticks()
    cb.set_label('Number of earthquakes')

    # Plot the logo.
    if os.path.isfile(param.logo_image):
        im = Image.open(param.logo_image)
        im.thumbnail((param.logo_width, param.logo_height), Image.ANTIALIAS)  # resizes image in-place
        fig.figimage(im, 0.80 * param.logo_x, param.logo_y, zorder=1000)

    # Post the bin size
    # on the figure total in percent [left, bottom, width, height]
    ax2 = fig.add_axes([0.8, 0.85, 0.1, 0.1])
    ax2.axis('off')
    plt.text(0, 0, f'{d_km} km bins', horizontalalignment='left', fontsize=8)

    # fig.tight_layout()

    # Print the figure to a .png file
    matplotlib.pyplot.savefig(plot_file_name, bbox_inches='tight')

    plt.close('all')

    print(f'[INFO] Done with figure {map_tag}: {plot_file_name} ', flush=True, file=log_file)

    return norm_min, norm_max


def location_animation_frame(magnitude, map_tag='aftershocks', map_type='location_animation'):
    """

     Group 6: Hourly animation of aftershocks.

    """
    # Some map parameters.
    map_size_lon = map_size_lat / math.cos(mainshock['latitude'] * math.pi / 180)
    lat_bounds = mainshock['latitude'] - map_size_lat, mainshock['latitude'] + map_size_lat
    lon_bounds = mainshock['longitude'] - map_size_lon, mainshock['longitude'] + map_size_lon

    frame_end_time = dp_time_end
    if frame_end_time > datetime.utcnow():
        frame_end_time = datetime.utcnow() + timedelta(seconds=0, minutes=0, hours=0, days=1)
    frame_interval = timedelta(seconds=0, minutes=0, hours=1)
    frame_start_time = dp_time_animation_start
    frame_time = dp_time_animation_start
    files_to_load = f'{mainshock_id}_{map_tag}_{map_type}'

    mainshock_frame_count = -1
    frame_count = -1
    # Remove image files when done.
    try:
        files_to_remove = glob.glob(f'{os.path.join(scratch_dir, files_to_load)}*.png')
        for this_file in files_to_remove:
            os.remove(this_file)
    except Exception as _er:
        print('[ERR] Failed to remove\n {_er}', flush=True, file=log_file)

    # This Basemap takes a while to draw, so save it as a pickle file and reload it in the loop.
    # File for use by Python's object serialization module pickle. Added time tag to make it unique to this run.
    pickle_file = os.path.join(scratch_dir, f'{mainshock_id}_{map_type}_{int(datetime.now().timestamp())}.pickle')

    fig0 = matplotlib.pyplot.figure(figsize=(figure_size[0], figure_size[1]), dpi=param.dpi)

    gs1 = gridspec.GridSpec(3, 3)
    gs1.update(left=0.05, right=0.9, bottom=0.05, top=0.9, wspace=0.01)
    ax1 = plt.subplot(gs1[:, :])
    m2 = Basemap(projection='merc', llcrnrlat=lat_bounds[0], urcrnrlat=lat_bounds[1], llcrnrlon=lon_bounds[0],
                 urcrnrlon=lon_bounds[1], lat_ts=mainshock['latitude'],
                 resolution=basemap_resolution)  # lat_ts = 20, lat of true scale

    pickle.dump(m2, open(pickle_file, 'wb'), -1)

    # clear the figure
    plt.close('all')

    warn_issued = False

    while frame_start_time <= frame_time <= frame_end_time:
        frame_count += 1

        # Set up the figure.
        fig = matplotlib.pyplot.figure(figsize=(figure_size[0], figure_size[1]), dpi=param.dpi)

        gs1 = gridspec.GridSpec(3, 3)
        gs1.update(left=0.05, right=0.9, bottom=0.05, top=0.9, wspace=0.01)
        ax1 = plt.subplot(gs1[:, :])

        # Production date time stamp.
        production_date = dp_lib.version_timestamp(script_version)
        ax1.text(0.01, 0.02, production_date, horizontalalignment='left', fontsize=5,
                 verticalalignment='top', transform=ax1.transAxes)

        # Read pickle back in and plot it again (should be much faster).
        m = pickle.load(open(pickle_file, 'rb'))
        m.ax1 = ax1

        # m.shadedrelief(alpha=location_map_relief_alpha, zorder=10)
        try:
            m.drawcoastlines(color='gray')
        except Exception as ex:
            print(f'[ERR] draw error {ex}')
            pass

        m.drawcountries(color='gray')
        m.drawstates(color='gray')
        m.etopo(scale=etopo_scale, alpha=location_map_relief_alpha, zorder=10)

        m.drawparallels(np.arange(-90., 91., lat_line_spacing), labels=[1, 0, 0, 0], fontsize=6, linewidth=0.3,
                        alpha=0.5)
        m.drawmeridians(np.arange(-180., 181., lat_line_spacing), labels=[0, 0, 0, 1], fontsize=6, linewidth=0.3,
                        alpha=0.5)
        # Set the plot title.
        plt.suptitle(title, fontsize=title_font_size)

        # Plot label on the upper right.
        date_range_label = f'from {dp_time_animation_start.strftime("%Y-%m-%d")} to ' \
            f'{query_time_end.split("T")[0].strip()}'
        plot_text = f'Other events within {int(max_radius_km)} km radius; {date_range_label}; ' \
            f'within {aftershocks_depth_radius} km depth; ' \
            f'M ≥ {magnitude}'
        ax1.set_title(plot_text, fontsize=label_font_size)

        # Plot M7 & M5 dots for reference.
        this_size = (key_mag_small['size'] / aftershock_scale_factor) ** 2.5
        plt.plot(key_mag_small['x'], key_mag_small['y'], 'o', markersize=this_size,
                 transform=ax1.transAxes,
                 fillstyle='none', color='k')
        plt.text(key_mag_small['text_x'], key_mag_small['text_y'],
                 f'M{key_mag_small["size"]}',
                 fontsize=label_font_size, color='k', va='center', transform=ax1.transAxes)

        this_size = (key_mag_large['size'] / aftershock_scale_factor) ** 2.5
        plt.plot(key_mag_large['x'], key_mag_large['y'], 'o', markersize=this_size,
                 transform=ax1.transAxes,
                 fillstyle='none', color='k')
        plt.text(key_mag_large['text_x'], key_mag_large['text_y'],
                 f'M{key_mag_large["size"]}, MT\'s use GCMT\nlocation, depth & magnitude',
                 fontsize=label_font_size, color='k', va='center', transform=ax1.transAxes)

        # Plot the aftershocks.
        for lon_index, elon in enumerate(events['longitude']):
            this_event_time = events['datetime'][lon_index]
            this_hour = (frame_time - this_event_time) / 3600.0

            if this_hour == 0:
                alpha = 1
            elif this_hour < 0:
                continue
            else:
                alpha = 1.0 - (this_hour / fade_duration)

            if alpha <= 0.001:
                continue
            elif alpha > 1:
                alpha = 1

            xmap, ymap = m(events['longitude'][lon_index], events['latitude'][lon_index])
            size = (events['magnitude'][lon_index] / aftershock_scale_factor) ** 2.5

            # GCMT events.
            not_gcmt = True
            if events['mt'][lon_index] is not None:
                try:
                    b = beach(events['mt'][lon_index], xy=(xmap, ymap),
                              width=((events['magnitude'][lon_index] * mag_scale['location_animation']
                                      / (4.6 * aftershock_scale_factor)) ** 2.5),
                              linewidth=0.5, alpha=alpha,
                              facecolor=aftershock_color[lon_index], edgecolor='k', axes=ax1)
                    b.set_zorder(5000)
                    ax1.add_collection(b)
                    not_gcmt = False
                except:
                    if not warn_issued:
                        print(f'[WARN] GCMP beach ball failed for event {lon_index}')
                        warn_issued = True
                    not_gcmt = True
            if not_gcmt:
                m2.plot(xmap, ymap, markerfacecolor=aftershock_color[lon_index], markeredgecolor='k',
                        markeredgewidth=marker_edge_width, alpha=alpha,
                        marker='o', zorder=500,
                        markersize=size)

        # -------------- plot yellow star MainShock GCMT (slightly transparent)

        if frame_time >= mainshock['datetime'] != 0:
            mainshock_frame_count += 1
            if mainshock_frame_count <= fade_duration:
                alpha = 1
            else:
                alpha = 0.3
            if mainshock['mt'] is not None:
                if mainshock['mt'][1] != 0:
                    b = beach(mainshock['mt'], xy=m(mainshock['longitude'], mainshock['latitude']),
                              width=map_size_lat * ((aftershock_scale_factor * 3.3 * mainshock['magnitude']) ** 2.5),
                              linewidth=0.5,
                              alpha=alpha,
                              facecolor=mainshock_color, edgecolor='k')
                    b.set_zorder(10000)
                    ax1.add_collection(b)
            else:
                xmap, ymap = m(mainshock['longitude'], mainshock['latitude'])
                size = (mainshock['magnitude']/ aftershock_scale_factor) ** 2.5
                m2.plot(xmap, ymap, markerfacecolor=mainshock_color, markeredgecolor='k',
                        markeredgewidth=marker_edge_width, alpha=alpha,
                        marker='o', zorder=10000,
                        markersize=size)

        # Draw the search area circle.
        x0, y0 = m(mainshock['longitude'], mainshock['latitude'])
        x1, y1 = m(mainshock['longitude'], mainshock['latitude'] + max_radius)
        radius_map = y1 - y0
        circle = plt.Circle((x0, y0), radius_map, color=elevation_bar_color, alpha=elevation_bar_alpha,
                            lw=2, fill=False)
        circle.set_zorder(10000)
        ax1.add_patch(circle)

        xmap, ymap = m(mainshock['longitude'], mainshock['latitude'])
        size = (mainshock['magnitude'] / aftershock_scale_factor) ** 2.5
        if mainshock['latitude'] is None:
            m.plot(xmap, ymap, marker='*', markerfacecolor='y', markeredgecolor='k', markersize=size, alpha=0.6)

        # ---- if using projection = 'merc'
        scale_bar_lon = lon_bounds[0] + 1.6 * map_size_lon
        scale_bar_lat = lat_bounds[0] + map_size_lat * 0.15
        m.drawmapscale(scale_bar_lon, scale_bar_lat, 0, 0, scalebar_length_km, units='km',
                       labelstyle='simple',
                       fillcolor1='w', fillcolor2=elevation_bar_axes_color,
                       fontcolor='#555555',
                       zorder=5, fontsize=label_font_size)

        # Plot the colorbar.
        # on the figure total in percent [left, bottom, width, height]
        ax = fig.add_axes([0.82, 0.2, 0.01, 0.55])
        cmap, norm = dp_lib.get_rgba(plt, min_event_depth, max_event_depth, None)
        cb = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='vertical')
        tick_locator = ticker.MaxNLocator(nbins=5)
        cb.ax.invert_yaxis()
        cb.locator = tick_locator
        cb.update_ticks()
        cb.set_label("depth (km)")

        # Plot the logo.
        if os.path.isfile(param.logo_image):
            im = Image.open(param.logo_image)
            im.thumbnail((param.logo_width, param.logo_height), Image.ANTIALIAS)  # resizes image in-place
            fig.figimage(im, 0.80 * param.logo_x, param.logo_y, zorder=1000)

        # Post the frame time
        # on the figure total in percent [left, bottom, width, height]
        ax2 = fig.add_axes([0.8, 0.85, 0.1, 0.1])
        ax2.axis('off')
        plt.text(0, 0, f'{frame_time.strftime("%Y-%m-%d")}{line_break}'
        f'{frame_time.strftime("%H")}:00 UTC', horizontalalignment='left', fontsize=8)

        fig.tight_layout()

        # Print the figure to a .png file
        file_name = f'{files_to_load}_{frame_count:05}.png'
        file_name = os.path.join(scratch_dir, file_name)
        matplotlib.pyplot.savefig(file_name, bbox_inches='tight')

        # Keep the first image for screen display
        if frame_count == 0:
            file_name = f'{os.path.join(video_dir, files_to_load)}.png'
            matplotlib.pyplot.savefig(file_name, bbox_inches='tight')

        plt.close('all')

        print(f'[INFO] Done with frame {frame_count:>5}', end='\r', flush=True, file=log_file)

        frame_time += frame_interval

    # Create video
    # Remove the previous video file if exists.
    try:
        files_to_remove = glob.glob(f'{os.path.join(video_dir, files_to_load)}*.mp4')
        for this_file in files_to_remove:
            os.remove(this_file)
    except Exception as _er:
        print(f'[ERR] Failed to remove file {this_file}\n {_er}', flush=True, file=log_file)

    # create the video.
    try:
        print('\n[INFO] Creating the video:', flush=True, file=log_file)
        # Apply -vf pad=ceil(iw/2)*2:ceil(ih/2)*2 filter to avoid eight not divisible by 2 (1644x1491)
        # error without rescaling.
        command = f'ffmpeg -i {os.path.join(scratch_dir, files_to_load)}_%05d.png -c:v libx264 -pix_fmt ' \
            f'yuv420p -crf 23 ' \
            f'-r {param.frames_per_second} ' \
            f'-vf pad=ceil(iw/2)*2:ceil(ih/2)*2 ' \
            f'-y {os.path.join(video_dir, files_to_load)}.mp4'.split()
        subprocess.call(command)
    except Exception as _er:
        print(f'[ERR] Command {command} failed\n {_er}', flush=True, file=log_file)

    # Remove image files when done.
    try:
        files_to_remove = glob.glob(f'{os.path.join(scratch_dir, files_to_load)}*.png')
        for this_file in files_to_remove:
            os.remove(this_file)
    except Exception as _er:
        print(f'[ERR] Failed to remove file {this_file}\n {_er}', flush=True, file=log_file)

    # Remove pickle files when done.
    try:
        os.remove(pickle_file)
    except Exception as _er:
        print(f'[ERR] Failed to remove pickle file {pickle_file}\n {_er}', flush=True, file=log_file)


def heatmap_animation_frame(magnitude, map_tag='aftershocks', map_type='heatmap_animation',
                            norm_min=None, norm_max=None):
    """

     Group 7: Heatmap animation of aftershocks.

    """

    draw_error = False
    # Some map parameters.
    map_size_lon = map_size_lat / math.cos(mainshock['latitude'] * math.pi / 180)
    lat_bounds = mainshock['latitude'] - map_size_lat, mainshock['latitude'] + map_size_lat
    lon_bounds = mainshock['longitude'] - map_size_lon, mainshock['longitude'] + map_size_lon

    frame_end_time = dp_time_end
    if frame_end_time > datetime.utcnow():
        frame_end_time = datetime.utcnow() + timedelta(seconds=0, minutes=0, hours=0, days=1)
    frame_interval = timedelta(seconds=0, minutes=0, hours=1)
    frame_start_time = dp_time_animation_start
    frame_time = dp_time_animation_start
    files_to_load = f'{mainshock_id}_{map_tag}_{map_type}'
    pickle_file = os.path.join(scratch_dir, f'{mainshock_id}_{map_type}_{int(datetime.now().timestamp())}.pickle')

    mainshock_frame_count = -1
    frame_count = -1

    # Remove image files when done.
    try:
        files_to_remove = glob.glob(f'{os.path.join(scratch_dir, files_to_load)}*.png')
        for this_file in files_to_remove:
            os.remove(this_file)
    except Exception as _er:
        print(f'[ERR] Failed to remove\n {_er}', flush=True, file=log_file)

    # This Basemap takes a while to draw, so save it as a pickle file and reload it in the loop.
    fig0 = matplotlib.pyplot.figure(figsize=(figure_size[0], figure_size[1]), dpi=param.dpi)

    gs1 = gridspec.GridSpec(3, 3)
    gs1.update(left=0.05, right=0.9, bottom=0.05, top=0.9, wspace=0.01)
    ax1 = plt.subplot(gs1[:, :])
    m2 = Basemap(projection='merc', llcrnrlat=lat_bounds[0], urcrnrlat=lat_bounds[1], llcrnrlon=lon_bounds[0],
                 urcrnrlon=lon_bounds[1], lat_ts=mainshock['latitude'],
                 resolution=basemap_resolution)  # lat_ts = 20, lat of true scale

    pickle.dump(m2, open(pickle_file, 'wb'), -1)

    # clear the figure
    plt.close('all')

    if norm_min is None or norm_max is None:
        density_val, norm_min, norm_max = event_density_range(m2, events['latitude'], events['longitude'], d_km)

    while frame_start_time <= frame_time <= frame_end_time:
        frame_count += 1
        # Set up the figure.
        fig = matplotlib.pyplot.figure(figsize=(figure_size[0], figure_size[1]), dpi=param.dpi)

        gs1 = gridspec.GridSpec(3, 3)
        gs1.update(left=0.05, right=0.9, bottom=0.05, top=0.9, wspace=0.01)
        ax1 = plt.subplot(gs1[:, :])

        # Production date time stamp.
        production_date = dp_lib.version_timestamp(script_version)
        ax1.text(0.01, 0.02, production_date, horizontalalignment='left', fontsize=5,
                 verticalalignment='top', transform=ax1.transAxes)

        # Read pickle back in and plot it again (should be much faster).
        m = pickle.load(open(pickle_file, 'rb'))
        m.ax1 = ax1

        # m.shadedrelief(alpha=location_map_relief_alpha, zorder=10)
        try:
            m.drawcoastlines(color='gray')
        except Exception as ex:
            print(f'[ERR] draw error {ex}')
            pass

        m.drawcountries(color='gray')
        m.drawstates(color='gray')
        m.etopo(scale=etopo_scale, alpha=heatmap_relief_alpha, zorder=10)

        m.drawparallels(np.arange(-90., 91., lat_line_spacing), labels=[1, 0, 0, 0], fontsize=6, linewidth=0.3,
                        alpha=0.5)
        m.drawmeridians(np.arange(-180., 181., lat_line_spacing), labels=[0, 0, 0, 1], fontsize=6, linewidth=0.3,
                        alpha=0.5)

        # Set the plot title.
        plt.suptitle(title, fontsize=title_font_size)

        # Ratio of y to x axes with 1.6 empirical corrections to keep beach balls as circles.
        ratio = 1.6 * (ax1.get_ylim()[1] - ax1.get_ylim()[0]) / (ax1.get_xlim()[1] - ax1.get_xlim()[0])

        # Plot label on the upper right.
        date_range_label = f'from {dp_time_animation_start.strftime("%Y-%m-%d")} to ' \
            f'{query_time_end.split("T")[0].strip()}'
        plot_text = f'Other events within {int(max_radius_km)} km radius; {date_range_label}; ' \
            f'within {aftershocks_depth_radius} km depth; ' \
            f'M ≥ {magnitude}'
        ax1.set_title(plot_text, fontsize=label_font_size)

        # Plot M7 & M5 dots for reference.
        this_size = (key_mag_small['size'] / aftershock_scale_factor) ** 2.5
        plt.plot(key_mag_small['x'], key_mag_small['y'], 'o', markersize=this_size,
                 transform=ax1.transAxes,
                 fillstyle='none', color='k')
        plt.text(key_mag_small['text_x'], key_mag_small['text_y'],
                 f'M{key_mag_small["size"]}',
                 fontsize=label_font_size, color='k', va='center', transform=ax1.transAxes)

        this_size = (key_mag_large['size'] / aftershock_scale_factor) ** 2.5
        plt.plot(key_mag_large['x'], key_mag_large['y'], 'o', markersize=this_size,
                 transform=ax1.transAxes,
                 fillstyle='none', color='k')
        plt.text(key_mag_large['text_x'], key_mag_large['text_y'],
                 f'M{key_mag_large["size"]}, MT\'s use GCMT\nlocation, depth & magnitude',
                 fontsize=label_font_size, color='k', va='center', transform=ax1.transAxes)

        # Plot heatmap on Basemap.

        lon_list = list()
        lat_list = list()
        density = None
        for lon_index, elon in enumerate(events['longitude']):
            if events['datetime'][lon_index] <= frame_time:
                lon_list.append(elon)
                lat_list.append(events['latitude'][lon_index])
        if lon_list:
            start_lon = math.floor(min(events['longitude']))
            end_lon = math.ceil(max(events['longitude']))
            start_lat = math.floor(min(events['latitude']))
            end_lat = math.ceil(max(events['latitude']))

            # form the bins
            ll_x, ll_y = m(start_lon, start_lat)
            ur_x, ur_y = m(end_lon, end_lat)
            x_km = (ur_x - ll_x) / 1000.0
            y_km = (ur_y - ll_y) / 1000.0

            # compute appropriate bins to aggregate data
            # nx is number of bins in x-axis, i.e. longitude
            # ny is number of bins in y-axis, i.e. latitude
            nx = math.ceil(x_km / d_km)
            ny = math.ceil(y_km / d_km)
            lon_bins = np.linspace(start_lon, end_lon, nx)
            lat_bins = np.linspace(start_lat, end_lat, ny)

            # aggregate the number of earthquakes in each bin, we will only use the density
            density, lat_edges, lon_edges = np.histogram2d(lat_list, lon_list, [lat_bins, lon_bins])

            # get the mesh for the lat and lon
            lon_bins_2d, lat_bins_2d = np.meshgrid(lon_bins, lat_bins)
            # convert the bin mesh to map coordinates:
            xs, ys = m(lon_bins_2d, lat_bins_2d)  # will be plotted using pcolormesh

            cmap = plt.get_cmap(param.heatmap_color_map)
            # plt.register_cmap(cmap=custom_map)
            # Here adding one row and column at the end of the matrix, so that
            # density has same dimension as xs, ys, otherwise, using shading='gouraud'
            # will raise error
            density = np.hstack((density, np.zeros((density.shape[0], 1))))
            density = np.vstack((density, np.zeros((density.shape[1]))))

            # Plot heatmap with the color map
            m.pcolormesh(xs, ys, density, cmap=cmap, vmin=norm_min, vmax=norm_max, shading='gouraud')

        # -------------- plot yellow star MainShock GCMT (slightly transparent)
        if mainshock_frame_count <= fade_duration:
            mainshock_alpha = 0.3
        else:
            mainshock_alpha = 0.3

        if mainshock['mt'] is None:
            xmap, ymap = m(mainshock['longitude'], mainshock['latitude'])
            size = (mainshock['magnitude']/ aftershock_scale_factor) ** 2.5
            m2.plot(xmap, ymap, markerfacecolor=mainshock_color, markeredgecolor='k',
                    markeredgewidth=marker_edge_width, alpha=mainshock_alpha,
                    marker='*', zorder=10000,
                    markersize=size)
        elif mainshock['datetime'] <= frame_time and \
                mainshock['mt'][0] != 0 and mainshock['mt'][1] != 0:

            mainshock_frame_count += 1
            b = beach(mainshock['mt'], xy=m(mainshock['longitude'], mainshock['latitude']),
                      width=map_size_lat * ((aftershock_scale_factor * mag_scale['heatmap_animation']
                                             * mainshock['magnitude']) ** 2.5), linewidth=0.5,
                      alpha=mainshock_alpha,
                      facecolor=mainshock_color, edgecolor='k')
            b.set_zorder(10000)
            ax1.add_collection(b)

        # Draw the search area circle.
        x0, y0 = m(mainshock['longitude'], mainshock['latitude'])
        x1, y1 = m(mainshock['longitude'], mainshock['latitude'] + max_radius)
        radius_map = y1 - y0
        circle = plt.Circle((x0, y0), radius_map, color=elevation_bar_color, alpha=elevation_bar_alpha,
                            lw=2, fill=False)
        ax1.add_patch(circle)

        xmap, ymap = m(mainshock['longitude'], mainshock['latitude'])
        size = (mainshock['magnitude'] / aftershock_scale_factor) ** 2.5
        if events['mt'][lon_index] is None:
            m.plot(xmap, ymap, marker='*', markerfacecolor='y', markeredgecolor='k', markersize=size, alpha=0.6)

        # ---- if using projection = 'merc'
        scale_bar_lon = lon_bounds[0] + 1.6 * map_size_lon
        scale_bar_lat = lat_bounds[0] + map_size_lat * 0.15
        m.drawmapscale(scale_bar_lon, scale_bar_lat, 0, 0, scalebar_length_km, units='km',
                       labelstyle='simple',
                       fillcolor1='w', fillcolor2=elevation_bar_axes_color,
                       fontcolor='#555555',
                       zorder=5, fontsize=label_font_size)

        # Plot the colorbar.
        # on the figure total in percent [left, bottom, width, height]
        if density is not None:
            ax = fig.add_axes([0.82, 0.2, 0.01, 0.55])
            norm = matplotlib.colors.Normalize(vmin=norm_min, vmax=norm_max)
            cb = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='vertical')
            tick_locator = ticker.MaxNLocator(nbins=5)
            cb.locator = tick_locator
            cb.update_ticks()
            cb.set_label('Number of earthquakes')

        # Plot the logo.
        if os.path.isfile(param.logo_image):
            im = Image.open(param.logo_image)
            im.thumbnail((param.logo_width, param.logo_height), Image.ANTIALIAS)  # resizes image in-place
            fig.figimage(im, 0.80 * param.logo_x, param.logo_y, zorder=1000)

        # Post the frame time and number of events
        # on the figure total in percent [left, bottom, width, height]
        ax2 = fig.add_axes([0.8, 0.80, 0.1, 0.1])
        ax2.axis('off')
        plt.text(0, 0, f'{frame_time.strftime("%Y-%m-%d")}{line_break}'
        f'{frame_time.strftime("%H")}:00 UTC{line_break}{line_break}{len(lon_list)} events{line_break}{d_km} km bins',
                 horizontalalignment='left', fontsize=8)

        fig.tight_layout()

        # Print the figure to a .png file
        file_name = f'{files_to_load}_{frame_count:05}.png'
        file_name = os.path.join(scratch_dir, file_name)
        matplotlib.pyplot.savefig(file_name, bbox_inches='tight')

        # Keep the first image for screen display
        if frame_count == 0:
            file_name = f'{os.path.join(video_dir, files_to_load)}.png'
            matplotlib.pyplot.savefig(file_name, bbox_inches='tight')

        plt.close('all')

        print(f'[INFO] Done with frame {frame_count:>5}', end='\r', flush=True, file=log_file)

        frame_time += frame_interval

    # Create video
    # Remove the previous video file if exists.
    try:
        files_to_remove = glob.glob(f'{os.path.join(video_dir, files_to_load)}*.mp4')
        for this_file in files_to_remove:
            os.remove(this_file)
    except Exception as _er:
        print(f'[ERR] Failed to remove file {this_file}\n {_er}', flush=True, file=log_file)

    # create the video.
    try:
        print('\nCreating the video:', flush=True, file=log_file)
        file_name = files_to_load
        # Apply -vf pad=ceil(iw/2)*2:ceil(ih/2)*2 filter to avoid eight not divisible by 2 (1644x1491)
        # error without rescaling.
        command = f'ffmpeg -i {os.path.join(scratch_dir, file_name)}_%05d.png -c:v libx264 -pix_fmt ' \
            f'yuv420p -crf 23 ' \
            f'-r {param.frames_per_second} ' \
            f'-vf pad=ceil(iw/2)*2:ceil(ih/2)*2 ' \
            f'-y {os.path.join(video_dir, files_to_load)}.mp4'.split()
        subprocess.call(command)
    except Exception as _er:
        print(f'[ERR] Command {command} failed\n {_er}', flush=True, file=log_file)

    # Remove image files when done.
    try:
        files_to_remove = glob.glob(f'{os.path.join(scratch_dir, files_to_load)}*.png')
        for this_file in files_to_remove:
            os.remove(this_file)
    except Exception as _er:
        print(f'[ERR] Failed to remove file {this_file}\n {_er}', flush=True, file=log_file)

    # Remove pickle files when done.
    try:
        os.remove(pickle_file)
    except Exception as _er:
        print(f'[ERR] Failed to remove pickle file {pickle_file}\n {_er}', flush=True, file=log_file)


"""Main code"""

try:
    options, remainder = getopt.getopt(sys.argv[1:], 'hvxr:e:m:b:a:l:p:s:T:M:',
                                       ['help', 'verbose', 'refid=', 'xdate', 'eid=', 'minmag=', 'before=', 'after=',
                                        'label=', 'plots=', 'scale=', 'title=', 'rmag='])
    for opt, arg in options:
        if opt in ('-h', '--help'):
            usage()
            sys.exit(2)
        elif opt in ('-v', '--verbose'):
            verbose = True
        elif opt in ('-r', '--refid'):
            reference_event_id = arg.strip()
            event_type = 'Reference event'
        elif opt in ('-x', '--xdate'):
            date_xaxis = True
        elif opt in ('-e', '--eid'):
            mainshock_id = arg.strip()
            event_type = 'Event'
        elif opt in ('-m', '--minmag'):
            min_mag = float(arg.strip())
        elif opt in ('-T', '--title'):
            title_text = arg.replace('"', '')
        elif opt in ('-b', '--before'):
            days_before = arg.strip()
            if days_before.lower() == 'none':
                days_before = None
            else:
                days_before = int(days_before)
        elif opt in ('-a', '--after'):
            days_after = int(arg.strip())
        elif opt in ('-l', '--label'):
            map_tag = arg.strip()
        elif opt in ('-s', '--scale'):
            scale = float(arg.strip())
        elif opt in ('-p', '--plots'):
            # Seismicity map.
            if arg.strip() == '0':
                plots = [0]
            # Regular map(s).
            else:
                plot = arg.strip().split(',')
                plots = list()
                for p in plot:
                    plots.append(int(p))
        elif opt in ('-M', '--rmag'):
            reference_mag = float(arg.strip())

except getopt.GetoptError as er:
    usage()
    print(f'\n\n{60 * "="}\n{er}\n{60 * "="}\n\n', sep='\n', flush=True)
    sys.exit(3)

# User has option of defining a mainshock or a reference event. The only difference is how plots are labeled.
# W decide on what type of event is based on the user inputs.
is_reference_event = False
if mainshock_id is not None and reference_event_id is not None:
    er = f'[ERR] Cannot define both mainshock ID and reference event ID.'
    usage()
    print(f'\n\n{60 * "="}\n{er}\n{60 * "="}\n\n', sep='\n', flush=True)
    sys.exit(4)
elif reference_event_id is not None:
    mainshock_id = reference_event_id
    is_reference_event = True
elif mainshock_id is None and reference_event_id is None:
    er = f'[ERR] Must provide either mainshock ID or reference event ID.'
    usage()
    print(f'\n\n{60 * "="}\n{er}\n{60 * "="}\n\n', sep='\n', flush=True)
    sys.exit(5)
elif plots != [0] and days_before is None:
    er = f'[ERR] a days before option "-b" is only valid with plot option "-p" of 0.\n' \
         f'     (-p 0 -b none)'
    usage()
    print(f'\n\n{60 * "="}\n{er}\n{60 * "="}\n\n', sep='\n', flush=True)
    sys.exit(6)

mainshock_fdsn, mainshock_gcmt = event_lib.get_event_by_id(mainshock_id, catalog_dc, scratch_dir, log_file,
                                                           gcmt_dc=gcmt_dc)

# Time windows before and after the mainshock.
if days_before is None:
    time_before = None
else:
    time_before = timedelta(seconds=0, minutes=0, hours=0, days=days_before)
    time_before_animation = timedelta(seconds=0, minutes=0, hours=0, days=days_before_animation)

time_after = timedelta(seconds=0, minutes=0, hours=0, days=days_after)

# If there is a GCMT solution for the mainshock, it overrides the FDSN solution. Save the needed info.
mainshock = dict()
if mainshock_gcmt is None:
    selected_mainshock = mainshock_fdsn
    mainshock['mt'] = None
    mainshock['strikes'] = None
else:
    selected_mainshock = mainshock_gcmt
    mainshock['mt'] = mainshock_gcmt['mt']
    mainshock['strikes'] = mainshock_gcmt['strikes']

# User has the option of not using the GCMT coordinates for the mainshock.
if not use_gcmt_for_mainshock:
    mainshock['latitude'] = mainshock_fdsn['origin'].latitude
    mainshock['longitude'] = mainshock_fdsn['origin'].longitude
    mainshock['depth'] = mainshock_fdsn['origin'].depth / 1000.0
else:
    mainshock['latitude'] = selected_mainshock['origin'].latitude
    mainshock['longitude'] = selected_mainshock['origin'].longitude
    mainshock['depth'] = selected_mainshock['origin'].depth / 1000.0

mainshock['magnitude'] = selected_mainshock['magnitude'].mag
# QuakeML depths are in meters.
mainshock['datetime'] = mainshock_fdsn['origin'].time
mainshock['date_time'] = mainshock_fdsn['origin'].time.strftime('%Y-%m-%dT%H:%M:%S')
mainshock['name'] = mainshock_fdsn['description'].text
mainshock['id'] = mainshock_fdsn['id']

# Map parameters.
map_size_lat, lat_line_spacing, scalebar_length_km, max_radius, min_mag, \
              max_mag, aftershock_scale_factor = dp_lib.set_map_parameters(mainshock['magnitude'])

max_radius_km = degrees2kilometers(max_radius)

# Adjust the heatmap bin width accordingly
d_km = max_radius_km / 20.0

# Etopo uses scale to reduce map resolution for maps that are zoomed in, we need to
# increase the scale. But this will use to much memory.
if scale is None:
    if max_radius_km < 200.0:
        etopo_scale = 2.6
    else:
        etopo_scale = 1.2
else:
    etopo_scale = scale

# Calculate points along the 2 nodal plane strikes.
if mainshock['strikes'] is not None:
    azimuth = mainshock['strikes']
    distance_along_strike = np.linspace(- map_size_lat, map_size_lat, param.num_points_along_strike)
    distance_along_strike_km = list()
    for deg in distance_along_strike:
        distance_along_strike_km.append(degrees2kilometers(deg))

    strike_line_1, strike_line_2, elevation = dp_lib.get_strike_elevation(mainshock, distance_along_strike, azimuth)
else:
    azimuth = None

# Time windows before and after the mainshock.
day_index = 0
query_time_start = time_before

if time_before is not None:
    dp_time_start = mainshock['datetime'] - time_before
    dp_time_animation_start = mainshock['datetime'] - time_before_animation
    query_time_start = dp_time_start.strftime('%Y-%m-%dT%H:%M:%S')

dp_time_end = mainshock['datetime'] + time_after
if dp_time_end > datetime.utcnow():
    dp_time_end = datetime.utcnow() + timedelta(seconds=0, minutes=0, hours=0, days=1)

query_time_end = dp_time_end.strftime('%Y-%m-%dT%H:%M:%S')

query_depth_start = mainshock['depth'] - aftershocks_depth_radius
query_depth_end = mainshock['depth'] + aftershocks_depth_radius
if query_depth_start < 0:
    query_depth_start = 0

# Set the plot title.
if title_text is None:
    title = f'{event_type}: M{mainshock["magnitude"]} {mainshock["name"]}\n{mainshock["date_time"]} UTC Lat: ' \
        f'{mainshock["latitude"]:.2f}° ' \
        f'Lon: {mainshock["longitude"]:.2f}° ' \
        f'depth: {mainshock["depth"]:.1f} km'

else:
    title = f'{title_text}{em_dash}{event_type.lower()}: M{mainshock["magnitude"]} {mainshock["name"]}' \
        f'\n{mainshock["date_time"]} UTC Lat: {mainshock["latitude"]:.2f}° ' \
        f'Lon: {mainshock["longitude"]:.2f}° ' \
        f'depth: {mainshock["depth"]:.1f} km'

if mainshock['mt'] is None:
    title = f'{title} (no focal mechanism)'

# Do maps for different map time intervals.
do_zero_mag = True

if days_before is None:
    print(f'[INFO] Working on historic seismicity', flush=True, file=log_file)
else:
    print(f'[INFO] Working {days_before} days before to {days_after} days after', flush=True, file=log_file)

dp_days_before = days_before
dp_days_after = days_after
if dp_days_before is not None:
    if dp_days_before + dp_days_after <= 14:
        plot_days_tick_interval = 1
    else:
        plot_days_tick_interval = math.ceil((dp_days_before + dp_days_after) / 10) + 1
date = None
if query_time_start is None:
    date_range_label = f'all times'
else:
    date_range_label = f'from {query_time_start.split("T")[0].strip()} to ' \
        f'{query_time_end.split("T")[0].strip()}'
print(f'[INFO] Working on {date_range_label}', flush=True, file=log_file)

fdsn_events, fdsn_url = event_lib.get_fdsn_events(catalog_dc, query_time_start, query_time_end, [], [], min_mag,
                                                  log_file,
                                                  min_depth=query_depth_start, max_depth=query_depth_end,
                                                  center_lat=mainshock['latitude'], center_lon=mainshock['longitude'],
                                                  radius=max_radius, return_url=True)

# Since GCMT search area is a rectangle and FDSN does it based on radius, we look for GCMT events in a larger area.
lat_limits = [mainshock['latitude'] - 2.0 * max_radius, mainshock['latitude'] + 2.0 * max_radius]

# For GCMT service the start time is 1976-1-1
this_start = query_time_start
if query_time_start is None:
    this_start = catalog_begin

lon_limits = [mainshock['longitude'] - 2.0 * max_radius, mainshock['longitude'] + 2.0 * max_radius]
gcmt_events, gcmt_url = event_lib.get_gcmt_events(this_start, query_time_end, lat_limits, lon_limits, min_mag, log_file,
                                                  max_depth=query_depth_end, min_depth=query_depth_start,
                                                  scratch_dir=scratch_dir,
                                                  return_url=True)
if gcmt_events is not None:
    # For now we are using the values given in the association variable to associate FDSN and GCMT events. If
    # there are too many missed event, increase the association_factor to reduce the sensitivity.
    association_index = event_lib.get_association_index(fdsn_events['origin'], 'neic', gcmt_events['origin'],
                                                        'gcmt', log_file, association_factor=1.0)
else:
    association_index = dict()


events = dict()
# Get the foreshocks, aftershocks and the mainshock.  Preference given to GCMT events.
events['date_time'], events['datetime'], events['latitude'], events['longitude'], events['depth'], \
    events['magnitude'], events['mt'] = \
    dp_lib.parse_event_data(fdsn_events, gcmt_events, association_index, log_file, event_id=mainshock_id,
                            unique_id=True, save_event_list=False)

if plot_unassociated_gcmts and mainshock_gcmt is not None and gcmt_events is not None:
    aux_events = dict()
    # These are GCMT events that did not get associated to an FDSN events but are still within the search radius.
    # We will just plot them.
    aux_events['date_time'], aux_events['datetime'], aux_events['latitude'], aux_events['longitude'], aux_events[
        'depth'], aux_events['magnitude'], aux_events['mt'] = \
        dp_lib.get_unassociated_gcmts(mainshock_gcmt, gcmt_events, association_index, max_radius, log_file)

    for aux_index, aux in enumerate(aux_events['latitude']):
        events['date_time'].append(aux_events['date_time'][aux_index])
        events['datetime'].append(aux_events['datetime'][aux_index])
        events['latitude'].append(aux_events['latitude'][aux_index])
        events['longitude'].append(aux_events['longitude'][aux_index])
        events['depth'].append(aux_events['depth'][aux_index])
        events['magnitude'].append(aux_events['magnitude'][aux_index])
        events['mt'].append(aux_events['mt'][aux_index])

print(f'[INFO] Found total of {len(events["date_time"])} events', flush=True, file=log_file)
events['time_diff_sec'] = list()
events['time'] = list()
events['distance'] = list()
events['distance_along_strike'] = [[], []]

# Event IDs with mag < min_mag.
# if len(events_id) == 1:
#    print('[WARN] EXITING- there are no aftershocks', flush=True, file=log_file)
#    sys.exit()

for event_index, lat in enumerate(events['latitude']):

    this_eventtime = events['datetime'][event_index]

    events['time_diff'] = this_eventtime - mainshock['datetime']
    seconds = events['time_diff']
    events['time_diff_sec'].append(seconds)

    events['time'].append(event_index)

    distance_degrees_1 = list()
    distance_degrees_2 = list()
    if mainshock['strikes']:
        for d_index, dist in enumerate(distance_along_strike):
            strike_1_lat_lon = strike_line_1[d_index]
            strike_2_lat_lon = strike_line_2[d_index]
            distance_degrees_1.append(dp_lib.great_circle_distance(strike_1_lat_lon[0], strike_1_lat_lon[1],
                                                                   events['latitude'][event_index],
                                                                   events['longitude'][event_index]))
            distance_degrees_2.append(dp_lib.great_circle_distance(strike_2_lat_lon[0], strike_2_lat_lon[1],
                                                                   events['latitude'][event_index],
                                                                   events['longitude'][event_index]))

        events['distance_along_strike'][0].append(distance_along_strike_km[np.argmin(distance_degrees_1)])
        events['distance_along_strike'][1].append(distance_along_strike_km[np.argmin(distance_degrees_2)])

# Set the range for depth based on the mainshock's depth and the search radius.
if len(events['depth']) > 0:
    max_event_depth = math.ceil(max(mainshock['depth'], max(events['depth'])) * 1.1)
    min_event_depth = math.floor(min(mainshock['depth'], min(events['depth'])) * 0.9)
else:
    max_event_depth = mainshock['depth'] + aftershocks_depth_radius
    min_event_depth = mainshock['depth'] - aftershocks_depth_radius

if min_event_depth < 0:
    min_event_depth = 0

print(f'[INFO] Depth range is from {min_event_depth} km to {max_event_depth} km', flush=True, file=log_file)

# Set up color maps.
mainshock_color = dp_lib.get_rgba(plt, min_event_depth, max_event_depth, mainshock['depth'])
if days_before is not None:
    mainshock_day_color = dp_lib.get_rgba(plt, - days_before, days_after, 0)

aftershock_color = list()
aftershock_time_color = list()

for lon_index, elon in enumerate(events['longitude']):
    aftershock_color.append(dp_lib.get_rgba(plt, min_event_depth, max_event_depth, events['depth'][lon_index]))
    if days_before is not None:
        aftershock_time_color.append(dp_lib.get_rgba(plt, - days_before, days_after,
                                                     events['time_diff_sec'][lon_index] / 86400))

# Create the requested aftershock maps and videos.
for i, plot in enumerate(plots):
    if int(plot) == 0:
        location_map(min_mag, map_tag=tag, map_type=map_index[0])
    if int(plot) == 1:
        magnitude_vs_days(min_mag, map_tag=tag)
    elif int(plot) == 2:
        if mainshock['mt'] is None:
            print(f'[WARN] No GCMT solution found, plot of {param.map_description[plot]} not generated', flush=True,
                  file=log_file)
            print(f'[WARN] No GCMT solution found, plot of {param.map_description[plot]} not generated', flush=True)
            continue
        strike_vs_day(events['distance_along_strike'], map_tag=tag)
    elif int(plot) == 3:
        if mainshock['mt'] is None:
            print(f'[WARN] No GCMT solution found, plot of {param.map_description[plot]} not generated', flush=True,
                  file=log_file)
            print(f'[WARN] No GCMT solution found, plot of {param.map_description[plot]} not generated', flush=True)
            continue
        strike_vs_depth(events['distance_along_strike'], map_tag=tag)
    elif int(plot) == 4:
        location_map(min_mag, q_start=query_time_start, q_end=query_time_end, map_tag=tag)
    elif int(plot) == 5:
        if len(events['latitude']) <= 1:
            print(f'[WARN] No GCMT solution found, plot of {param.map_description[plot]} not generated', flush=True,
                  file=log_file)
            print(f'[WARN] No GCMT solution found, plot of {param.map_description[plot]} not generated', flush=True)
            continue
        norm_min, norm_max = heatmap(min_mag, map_tag=tag)
    elif int(plot) == 6:
        location_animation_frame(min_mag, map_tag=tag, map_type=map_index[6])
    elif int(plot) == 7:
        if len(events['latitude']) <= 1:
            print(f'[WARN] No GCMT solution found, plot of {param.map_description[plot]} not generated', flush=True,
                  file=log_file)
            print(f'[WARN] No GCMT solution found, plot of {param.map_description[plot]} not generated', flush=True)
            continue
        heatmap_animation_frame(min_mag, map_tag=tag, map_type=map_index[7], norm_min=norm_min,
                                norm_max=norm_max)
