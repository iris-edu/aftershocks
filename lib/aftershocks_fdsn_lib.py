#!/usr/bin/env python
import os
import sys
import math

import numpy as np
from datetime import datetime

import netCDF4

import requests
import urllib
from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError

import matplotlib
import matplotlib.cm as cmx

import collections

from obspy.geodetics import degrees2kilometers, kilometer2degrees

# Import the aftershocks parameters and libraries.
_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
param_dir = os.path.join(_dir, 'param')

sys.path.append(param_dir)

import aftershock_fdsn_maps_param as param

"""
    Description:

    A Python utility library used by the main script.

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
        2020-09-16 Manoch: V.2020.260 public release
        2020-08-22 Manoch: V.2020.236 FDSN support
        2020-08-01 Manoch: V.2020.214 release.

"""

spud_service_gcmt_url = param.spud_service_gcmt_url

radius_deg = param.radius_deg

catalog_dc = param.catalog_dc
gcmt_dc = param.gcmt_dc


def mkdir(target_directory):
    """ Make a directory if it does not exist."""
    directory = None
    try:
        directory = target_directory
        if not os.path.exists(directory):
            os.makedirs(directory)
        return directory
    except Exception as _er:
        print(f'[ERR] failed to create directory {directory}\n{_er}', flush=True)
        return None


def get_rgba(plt, min_val, max_val, val, cmap=None):
    """Return colormap or color value."""

    if cmap is None:
        cmap = cm = plt.get_cmap(param.map_color_map)

    norm = matplotlib.colors.Normalize(vmin=min_val, vmax=max_val)
    if val is None:
        return cmap, norm
    scalar_color_map = cmx.ScalarMappable(norm=norm, cmap=cmap)
    return scalar_color_map.to_rgba(val)


def now():
    """Return the current time"""

    return  datetime.now().strftime('%Y-%m-%d %H:%M:%S')


def number_of_days(from_date, to_date):
    """Return number of days between two dates."""
    return (to_date - from_date).days


def event_marker_size(event_mag, scale_factor):
    """Calculates an event marker size based on the magnitude of the event and the scale factor"""
    marker_size = (event_mag / scale_factor) ** 2.5

    return marker_size


def get_index(value, value_array):
    """
    Search for the nearest (value) in an array of values and return
    the position index of the closest value.
     """
    value_index = (np.abs(value_array - value)).argmin()
    return value_index


def get_directions(az):
    """set the the N/S/E/W directions based on the given azimuth."""
    direction_1 = ''
    direction_2 = ''

    if az < 23 or az > 338:
        direction_1 = 'N'
        direction_2 = 'S'
    elif 68 > az >= 23:
        direction_1 = 'NE'
        direction_2 = 'SW'
    elif 113 > az >= 68:
        direction_1 = 'E'
        direction_2 = 'W'
    elif 158 > az >= 113:
        direction_1 = 'SE'
        direction_2 = 'NW'
    elif 203> az >= 158:
        direction_1 = 'S'
        direction_2 = 'N'
    elif 248 > az >= 203:
        direction_1 = 'SW'
        direction_2 = 'NE'
    elif 294 > az >= 248:
        direction_1 = 'W'
        direction_2 = 'E'
    elif 338 >= az >= 294:
        direction_1 = 'NW'
        direction_2 = 'SE'

    return direction_1, direction_2
        

def days_ticks(use_ratio, ticks_before, ticks_after, ticks_interval):
    """Create tick marks for the days axis."""
    ticks_1 = np.dot(use_ratio, np.arange(ticks_interval, ticks_before + ticks_interval, step=ticks_interval))
    ticks_1 *= -1
    ticks_2 = np.dot(use_ratio, np.arange(0, ticks_after + ticks_interval, step=ticks_interval))
    ticks = np.concatenate((ticks_1, ticks_2))

    tick_locations = ticks

    labels_1 = np.arange(ticks_interval, ticks_before + ticks_interval, step=ticks_interval)
    labels_1 *= -1
    labels_2 = np.arange(0, ticks_after + ticks_interval, step=ticks_interval)
    labels = np.concatenate((labels_1, labels_2))

    tick_labels = labels
    return tick_locations, tick_labels


def get_color(this_value, color_range):
    """Find a color for a given value."""

    color_index = int(254 * (this_value - color_range[0]) / (color_range[1] - color_range[0]))
    color_index = min(color_index, 253)
    red_color = param.map_color_map_r[color_index]
    green_color = param.map_color_map_g[color_index]
    blue_color = param.map_color_map_b[color_index]

    return red_color, green_color, blue_color


class Event(object):
    def __init__(self, e_id=None, e_time=None, e_hour=None, e_datetime=None, e_latitude=None, e_longitude=None,
                 e_depth=None, e_country=None, e_region=None, e_cat=None, e_auth=None, e_mag=None,
                 g_id=None, g_time=None, g_lat=None, g_lon=None, g_depth=None, g_mag=None, g_tensor=None,
                 g_np1=None, g_np2=None):
        self.id = e_id
        self.date_time = e_time
        self.date_hour = e_hour
        self.datetime = e_datetime
        self.latitude = e_latitude
        self.longitude = e_longitude
        self.depth = e_depth
        self.country = e_country
        self.region = e_region
        self.catalog = e_cat
        self.author = e_auth
        self.magnitude = e_mag
        self.g_id = g_id
        self.g_time = g_time
        if g_lat is None:
            self.g_lat = e_latitude
            self.g_lon = e_longitude
            self.g_depth = e_depth
            self.g_mag = e_mag
        else:
            self.g_lat = g_lat
            self.g_lon = g_lon
            self.g_depth = g_depth
            self.g_mag = g_mag
        self.g_tensor = g_tensor
        self.g_np1 = g_np1
        self.g_np2 = g_np2

    def __str__(self):
        return str(self.__class__) + ": " + str(self.__dict__)


def version_timestamp(version):
    current_time = datetime.utcnow()
    timestamp = f'{param.production_label} {version}: {current_time.strftime("%Y-%m-%d %H:%M")} UTC'
    return timestamp


def parse_event_data(fdsn_events, gcmt_events, association, log_file, event_id=None,
                     unique_id=False, save_event_list=False):
    """Parse the lines returned from event service"""

    e_date_time = list()
    e_datetime = list()
    lat = list()
    lon = list()
    depth = list()
    mag = list()
    mt = list()

    fp = None
    if save_event_list and event_id:
        fp = open(os.path.join(param.data_dir, f'{event_id}_events_list.txt'), 'w')
        fp.write('# EventID|Time|Latitude|Longitude|Depth/km|Author|Catalog|Contributor|ContributorID|'
                 'MagType|Magnitude|MagAuthor|EventLocationName\n')
    for this_index, line in enumerate(fdsn_events['origin']):
        if event_id in fdsn_events['id'][this_index]:
            print(f'[INFO] Skipped event {event_id} as requested!', flush=True, file=log_file)
            continue
        if gcmt_events is not None and unique_id and this_index in association['neic']:
            _index = association['neic'].index(this_index)
            this_origin = gcmt_events['origin'][association['gcmt'][_index]]
            print(f'[INFO] Selecting GCMT origin\n\n{this_origin}\n\n'
                  f'Over FDSN origin\n\n{fdsn_events["origin"][this_index]}', flush=True, file=log_file)
            this_magnitude = gcmt_events['magnitude'][association['gcmt'][_index]]
            mt.append(gcmt_events['mt'][association['gcmt'][_index]])
        else:
            this_origin = fdsn_events['origin'][this_index]
            this_magnitude = fdsn_events['magnitude'][this_index]
            mt.append(None)

        # Save the event list
        if fp is not None:
            fp.write(f'{line.strip()}\n')

        e_date_time.append(this_origin.time.strftime('%Y-%m-%dT%H:%M:%S'))
        e_datetime.append(this_origin.time)
        lat.append(this_origin.latitude)
        lon.append(this_origin.longitude)
        depth.append(this_origin.depth / 1000.0)
        mag.append(this_magnitude.mag)

    if fp is not None:
        fp.close()

    return e_date_time, e_datetime, lat, lon, depth, mag, mt


def get_unassociated_gcmts(mainshock_gcmt, gcmt_events, association_index, max_radius, log_file):
    """Parse the lines returned from event service and select the unassociated GCMT events"""

    e_date_time = list()
    e_datetime = list()
    lat = list()
    lon = list()
    depth = list()
    mag = list()
    mt = list()

    fp = None

    for this_index, line in enumerate(gcmt_events['origin']):
        if this_index not in association_index['gcmt']:
            gcmt_origin = gcmt_events['origin'][this_index]
            dist_deg = great_circle_distance(mainshock_gcmt['origin'].latitude, mainshock_gcmt['origin'].longitude,
                                             gcmt_origin.latitude, gcmt_origin.longitude, mode='degrees')
            # Exclude the mainshock
            if dist_deg <= max_radius and gcmt_origin != mainshock_gcmt['origin']:
                print(f'[INFO] Selecting unassociated GCMT origin\n\n{gcmt_origin}\n\n', flush=True, file=log_file)
                this_magnitude = gcmt_events['magnitude'][this_index]
                mt.append(gcmt_events['mt'][this_index])

                e_date_time.append(gcmt_origin.time.strftime('%Y-%m-%dT%H:%M:%S'))
                e_datetime.append(gcmt_origin.time)
                lat.append(gcmt_origin.latitude)
                lon.append(gcmt_origin.longitude)
                depth.append(gcmt_origin.depth / 1000.0)
                mag.append(this_magnitude.mag)

    return e_date_time, e_datetime, lat, lon, depth, mag, mt


def define_color_map():
    """Define a custom color map."""

    # --- Rainbow_r
    color_map_r = param.color_map_r
    color_map_g = param.color_map_g
    color_map_b = param.color_map_b

    a = np.array(color_map_r)
    b = np.array(color_map_g)
    c = np.array(color_map_b)

    color_map = np.vstack((a, b, c))

    return color_map, color_map_r, color_map_g, color_map_b


def set_map_parameters(mag):
    """Set the map parameters based on the event magnitude."""

    # Min mag for aftershocks.
    min_mag = 3.1
    max_mag = 9.8
    aftershock_scale_factor = 2
    if mag > 9:
        map_size_lat = 10
        lat_line_spacing = 5
        scalebar_length_km = 500
        max_radius = 1000
    elif mag >= 8.5:
        map_size_lat = 5
        lat_line_spacing = 2
        scalebar_length_km = 300
        max_radius = 500
    elif mag >= 8:
        map_size_lat = 4
        lat_line_spacing = 2
        scalebar_length_km = 200
        max_radius = 400
    elif mag >= 7:
        map_size_lat = 3
        lat_line_spacing = 2
        scalebar_length_km = 100
        max_radius = 300
    elif mag >= 6:
        map_size_lat = 2
        lat_line_spacing = 1
        scalebar_length_km = 100
        max_radius = 200
        min_mag = 1.8
        max_mag = 8.5
        aftershock_scale_factor = 1.9
    else:
        map_size_lat = 1
        lat_line_spacing = 0.5
        scalebar_length_km = 50
        max_radius = 100
        max_mag = 8.0
        min_mag = 1.3
        aftershock_scale_factor = 1.8

    max_radius = kilometer2degrees(max_radius)
    return map_size_lat, lat_line_spacing, scalebar_length_km, max_radius, min_mag, max_mag, aftershock_scale_factor


def get_page(url, code=False):
    """Extract the content of HTML page or return HTML code if code=True"""
    req = Request(url)
    try:
        this_response = urlopen(req)
        if code:
            return this_response.code, requests.get(url)
    except HTTPError as error:
        print('URL: {}\n'.format(url))
        print('The server couldn\'t fulfill the request.\n')
        print('Error code: {}\n'.format(error.code))
        if code:
            return None, None
        else:
            return None
    except URLError as error:
        print('URL: {}\n'.format(url))
        print('Failed to reach server.\n')
        print('Reason: {}\n'.format(error.reason))
        if code:
            return None, None
        else:
            return None
    except Exception as er:
        print('URL: {}\n'.format(url))
        print('Error: {}\n'.format(er))
        if code:
            return None, None
        else:
            return None
    else:
        print(f'[INFO] connected to {this_response.geturl()}\n')

    mybytes = this_response.read()
    this_response.close()

    this_string = mybytes.decode('ascii').split('\n')
    return this_string


def spud_response_to_text(this_response):
    """Get the formal title and Event Name from SPUD, however, SPUDiSERVICE does not
       return a true XML! So, we get it as a text and get the info from it. If SPUD changes,
       we need to update this. Manoch 2020-06-17
    """
    response_list = list()
    try:
        for this_line in this_response.content.splitlines():
            # Replace the &lt; and &gt; that used to be < and > of the XML tags with | so
            # we can separate tags and values
            this_data = this_line.decode('ascii').replace('&lt;', '|').replace('&gt;', '|').strip()
            response_list.append(this_data)
    except Exception as ex:
        print(f'[EXP] spud_response_to_text: {ex}')

    return response_list


def get_spud_moment_tensor_id(this_data):
    """Search SPUD's data and look for the Moment tensor ID.
    """
    for item in this_data:
        if 'MomentTensor' in item:
            moment_tensor_id = item.split('"')[1]

            return moment_tensor_id
    return None


def parse_spud_response(this_data):
    """Get the formal title and Event Name from SPUD, however, SPUDiSERVICE does not
       return a true XML! So, we get it as a text and get the infor from it. If SPUD changes,
       we need to update this. Manoch 2018-07-16
    """

    tag_list = ['|Time|', '|EventName|', '|NP1|', '|NP2|', '|Dip|', '|Strike|', '|Rake|', '|Mpp', '|Mrp',
                '|Mrr', '|Mrt', '|Mtp', '|Mtt', 'MomentTensor Exponent', '|Description|', '|Magnitude|',
                '|Latitude|', '|Longitude|', '|Depth']
    moment_items = collections.OrderedDict()
    nodal_plane = {'np1': {'dip': None, 'rake': None, 'strike': None},
                   'np2': {'dip': None, 'rake': None, 'strike': None}}
    for this_line in this_data:
        for tag in tag_list:
            if tag in this_line:
                if tag == '|Description|':
                    this_description = this_line.strip().split('|')[2].replace('Moment Tensor for', '')
                    this_description = this_description.strip().replace('  ', ' ')
                    moment_items['Description'] = this_description
                elif tag == '|NP1|':
                    plane = 'np1'
                elif tag == '|NP2|':
                    plane = 'np2'
                elif tag == '|Dip|':
                    nodal_plane[plane]['dip'] = this_line.strip().split('|')[2]
                elif tag == '|Rake|':
                    nodal_plane[plane]['rake'] = this_line.strip().split('|')[2]
                elif tag == '|Strike|':
                    nodal_plane[plane]['strike'] = this_line.strip().split('|')[2]
                elif tag == 'MomentTensor Exponent':
                    this_key = this_line.replace('"', '').strip()
                    moment_items['exponent'] = this_key.split('=')[1]
                elif tag in this_line:
                    value = this_line.strip().split('|')[2]
                    key = tag.replace('|', '').strip()
                    if key not in moment_items.keys():
                        moment_items[key] = value
    return moment_items, nodal_plane


def great_circle_distance(lat_1, lon_1, lat_2, lon_2, mode='degrees'):
    """Calculate the Great Circle ditsnce"""
    
    theta_1 = (90.0 - lat_1) * radius_deg
    theta_2 = (90.0 - lat_2) * radius_deg
    phi_1 = lon_1 * radius_deg
    phi_2 = lon_2 * radius_deg
    angle = math.acos(min(1.0, math.sin(theta_1) * math.sin(theta_2) *
                          math.cos(phi_2 - phi_1) + math.cos(theta_1) * math.cos(theta_2)))
    distance_deg = angle / radius_deg
    distance_km = degrees2kilometers(distance_deg)
    if mode == 'degrees':
        return distance_deg
    if mode == 'km':
        return distance_km

    return distance_deg, distance_km


def great_circle_points(lat_1, lon_1, distance_deg, azi):
    """ Get the Lat/Lon of points along a great circle given azi + dist."""
    if distance_deg == 0:
        lat_2 = lat_1
        lon_2 = lon_1
        return lat_2, lon_2

    distance_deg_radius = distance_deg * radius_deg
    azimuth_radius = azi * radius_deg
    theta_1 = (90. - lat_1) * radius_deg
    phi_1 = lon_1 * radius_deg
    ctheta_2 = math.sin(distance_deg_radius) * math.sin(theta_1) * math.cos(azimuth_radius) + \
               math.cos(theta_1) * math.cos(distance_deg_radius)
    theta_2 = math.acos(ctheta_2)
    if theta_1 == 0:
        phi_2 = azimuth_radius
    elif theta_2 == 0:
        phi_2 = 0.
    else:
        sphi_2 = math.sin(distance_deg_radius) * math.sin(azimuth_radius) / math.sin(theta_2)
        cphi_2 = (math.cos(distance_deg_radius) - math.cos(theta_1) *
                  ctheta_2) / (math.sin(theta_1) * math.sin(theta_2))
        phi_2 = phi_1 + math.atan2(sphi_2, cphi_2)
    lat_2 = 90. - theta_2 / radius_deg
    lon_2 = phi_2 / radius_deg
    #if lon_2 > 360:
    #    lon_2 -= 360.0
    #if lon_2 < 0:
     #   lon_2 += 360.0

    return lat_2, lon_2


def get_event_info_from_file(event_file, verbose=False):
    """Read an event file and extract the event information"""

    if not os.path.exists(event_file):
        print(f'[ERR] Could not find the event file: {event_file}', flush=True)
        return Event()

    with open(event_file, 'r') as fp:
        lines = fp.read().split('\n')

    line = lines[0].strip()
    if verbose:
        print(f'[INFO] Event line: {lines}', flush=True)

    if not len(lines):
        print(f'[ERR] Event file {event_file} is empty!', flush=True)
        return Event()

    # Split the event file line to get the event information.
    try:
        items = line.split('|')
        event_id, year, month, day, hour, minute, second, event_latitude, \
        event_longitude, event_depth, event_country, not_used, not_used, \
        event_area, event_cat, event_auth, mag_info = items[0:17]

        mag = mag_info.split(',')
        e_date_time = f'{year}-{month}-{day}T{hour}:{minute}:{second}'
        e_date_hour = f'{year}-{month}-{day}T{hour}'
        mag_info = {'type': mag[0], 'value': float(mag[1]), 'author': mag[2]}

        event = Event(e_id=event_id, e_time=e_date_time, e_hour=e_date_hour, e_latitude=float(event_latitude),
                      e_longitude=float(event_longitude),
                      e_depth=float(event_depth), e_country=event_country.strip(), e_region=event_area.strip(),
                      e_cat=event_cat.strip(), e_auth=event_auth.strip(), e_mag=mag_info['value'])
    except Exception as ex:
        print(f'[ERR] Bad event record {line}\n{ex}', flush=True)
        event = Event()

    return event


def get_event_info_from_web(base_url, event_id):
    """Requests data from the event services."""

    event_url = f'{base_url}format=text&eventid={event_id}&nodata=404'

    print(f'[INFO] requesting events from {event_url}', flush=True)
    try:
        link = urllib.request.urlopen(event_url)
        event_data = link.read().decode().split('\n')
        link.close()
    except Exception as ex:
        print(f'[ERR] Request failed {ex}')
        return None

    event_id, event_date_time, event_datetime, event_lat, event_lon, event_depth, event_author, \
    event_catalog, event_contributor, \
    event_contributor_id, event_mag_type, event_mag, \
    vents_mag_author, event_location = parse_event_data(event_data, unique_id=True)

    event_date_hour = list()
    for t_index, t in enumerate(event_date_time):
        t = t.split(':')[0]
        event_date_hour.append(f'{t}:00:00')

    event = Event(e_id=event_id[0], e_time=event_date_time[0], e_datetime=event_datetime,
                  e_hour=event_date_hour, e_latitude=float(event_lat[0]),
                  e_longitude=float(event_lon[0]),
                  e_depth=float(event_depth[0]), e_country=None, e_region=event_location[0].strip(),
                  e_cat=event_catalog[0].strip(), e_auth=event_author[0].strip(), e_mag=event_mag[0])
    
    return event


def get_gcmt_info(event_id, verbose=False):

    # Get the GCMT information for the event.
    spud_event_url = f'{spud_service_gcmt_url}?eventid={event_id}'
    if verbose:
        print(f'[INFO] Requesting: {spud_event_url}')
    mystr, response = get_page(spud_event_url, code=True)
    spud_event_data = spud_response_to_text(response)
    spud_moment_tensor_id = get_spud_moment_tensor_id(spud_event_data)

    # Get the MT data.
    if spud_moment_tensor_id is None:
        print(f'[WARN] get_gcmt_info: SPUD does not have the GCMT information for event {event_id} yet!')
        return None, None
    spud_mt_url = f'{spud_service_gcmt_url}/{spud_moment_tensor_id}'
    if verbose:
        print(f'[INFO] Requesting: {spud_mt_url}')
    mystr, response = get_page(spud_mt_url, code=True)
    spud_mt_data = spud_response_to_text(response)
    moment_items, nodal_plane = parse_spud_response(spud_mt_data)

    return moment_items, nodal_plane


def get_strike_elevation(event, distance_deg, azim):
    """Extract the elevation data from the topo file."""

    topo_file = param.topo_file
    if not os.path.isfile(topo_file):
        print(f'[ERR] Topo file {topo_file} does not exist')
        sys.exit()

    nc_data = netCDF4.Dataset(topo_file)
    nc_lons = nc_data.variables[param.topo_file_x][:]
    nc_lats = nc_data.variables[param.topo_file_y][:]
    nc_elev = nc_data.variables[param.topo_file_z][:]
    nc_data.close()

    strike_line_1 = list()
    strike_line_2 = list()
    elevation_line_1 = list()
    elevation_line_2 = list()
    elevation = list()

    for deg in distance_deg:
        lat_1, lon_1 = great_circle_points(event['latitude'], event['longitude'], deg, azim[0])
        strike_line_1.append((lat_1, lon_1))
        lat_idx = get_index(lat_1, nc_lats)
        lon_idx = get_index(lon_1, nc_lons)
        elevation_line_1.append(nc_elev[lat_idx, lon_idx])

        lat_2, lon_2 = great_circle_points(event['latitude'], event['longitude'], deg, azim[1])
        strike_line_2.append((lat_2, lon_2))
        lat_idx = get_index(lat_2, nc_lats)
        lon_idx = get_index(lon_2, nc_lons)
        elevation_line_2.append(nc_elev[lat_idx, lon_idx])
    elevation.append(elevation_line_1)
    elevation.append(elevation_line_2)

    return strike_line_1, strike_line_2, elevation


def get_eq_browser_url(start_time, end_time, min_lat, max_lat, min_lon, max_lon, min_mag):
    """Mint an earthquake browser URL for the event."""

    eq_browser_aftershocks_url = f'https://ds.iris.edu/ieb/index.html?format=text&nodata=404' \
        f'&starttime={start_time}&endtime={end_time}&minmag={min_mag}' \
        f'&maxmag=10&mindepth=0&pbl=1' \
        f'&maxdepth=900&orderby=time-desc&&maxlat={max_lat}&minlat={min_lat}&maxlon={max_lon}' \
        f'&minlon={min_lon}&zm=5&mt=ter'

    aftershocks_url = f'<a href="{eq_browser_aftershocks_url} target="_NEW">{eq_browser_aftershocks_url}</a>'

    eq_browser_seismicity_url = f'https://ds.iris.edu/ieb/index.html?format=text&nodata=404' \
        f'&minmag={min_mag}' \
        f'&maxmag=10&mindepth=0&pbl=1' \
        f'&maxdepth=900&orderby=time-desc&src=usgs&limit=1000&maxlat={max_lat}&minlat={min_lat}&maxlon={max_lon}' \
        f'&minlon={min_lon}&zm=5&mt=ter'

    seismicity_url = f'<a href="{eq_browser_seismicity_url} target="_NEW">{eq_browser_seismicity_url}</a>'

    return eq_browser_seismicity_url, eq_browser_aftershocks_url

