import sys
import os

from datetime import datetime, timedelta
from obspy.core.event import read_events
import obspy
from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError
import uuid

"""
    Description:

    A Python FDSN/GCMT event library for use with the Aftershocks data product.

    Copyright and License:

    This software Copyright (c) 2021 IRIS (Incorporated Research Institutions for Seismology).

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
        2021-03-04 Manoch: v.2021.063 check for lat/lon limits when making GCMT requests.
                                      check for missing description in QuakeML.
                                      Improved GCMT event association checks.
        2021-02-09 Manoch: v.2021.040  r2.1 public release
        2020-08-22 Manoch: v.2020.236 FDSN
"""

# Association parameters
association_threshold = {'latitude': 0.5, 'longitude': 0.5, 'seconds': 30.0, 'depth': 30.0, 'magnitude': 0.2,
                         'gcmt_factor': 2.0}
event_service_url = {'ISC': 'http://www.isc.ac.uk/fdsnws/event/1/query?',
                     'NEIC': 'https://earthquake.usgs.gov/fdsnws/event/1/query?',
                     'IRIS': 'http://service.iris.edu/fdsnws/event/1/query?',
                     'GCMT': 'https://www.globalcmt.org/cgi-bin/globalcmt-cgi-bin/CMT5/form?'}


def get_cmt_id_from_resource_id(res_id):
    """Parse a QuakeML resource ID to get CMT ID"""
    tag = 'smi:local/cmtsolution/'
    cmt_id = None
    if res_id.strip().startswith(tag):
        res_id = res_id.replace(tag, '')
        cmt_id = (res_id.split('/'))[0]
    if cmt_id is None:
        print(f'[ERR] failed to get CMT ID from {res_id}')

    return cmt_id


def retrieve_gcmt_events(gcmt_url, log_file, eid=None):
    """Extract the GCMT event list from url or single line that contains the id provided"""
    req = Request(gcmt_url)
    try:
        this_response = urlopen(req)
    except HTTPError as ex:
        print(f'[ERR] URL: {gcmt_url}\nserver couldn\'t fulfill the request.\n{ex}', flush=True, file=log_file)
        return None
    except URLError as ex:
        print(f'[ERR] URL: {gcmt_url}\nfailed to reach server.\n{ex}', flush=True, file=log_file)
        return None
    except Exception as ex:
        print(f'[ERR] URL: {gcmt_url}\nfailed\n{ex}', flush=True, file=log_file)
        return None
    else:
        print(f'[INFO] Connected to: {this_response.geturl()}\n', flush=True, file=log_file)

    _data = this_response.read()
    this_response.close()
    try:
        lines = _data.decode('ascii').split('\n')
    except Exception as ex:
        print(f'[ERR] _data.decode\nfailed to convert {_data}.\n{ex}', flush=True, file=log_file) 
        return None

    events = ''

    # Looking for the id line.
    if eid is not None:
        print(f'[INFO] Looking for GCMT ID {eid}', flush=True, file=log_file)
        cmt_id = get_cmt_id_from_resource_id(eid)
        for line in lines:
            if cmt_id in line:
                line = line.strip()
                lon, lat, str1, dip1, rake1, str2, dip2, rake2, sc, iexp, name = line.split()
                return float(str1), float(str2)
        return None, None

    # Parsing the page
    _start = False
    _stop = False
    next_page = None
    for line in lines:
        if 'Output in CMTSOLUTION format' in line:
            _start = True
        elif 'End of events found with given criteria' in line:
            _start = False
            _stop = True
        elif 'More solutions' in line:
            next_page = line.split('start=')[1].split('"')[0]
            return events, next_page
        elif _start:
            # Skip lines with HTML tags.
            if '<' in line and '>' in line:
                continue

            # Skip the beginning blank lines.
            if not events.strip():
                events = line
            # Add this even information.
            else:
                events = f'{events}\n{line}'
        elif _stop:
            if not events:
                return None, None
            return events, next_page
    # Return None if no events found.
    if not events:
        return None, None
    else:
        return events, next_page


def get_catalog_event(in_catalog, log_file, output=False, verbose=False):
    """"Print quakeML catalog in readable format"""
    event = dict()
    # Event description.
    event['id'] = in_catalog.resource_id.id

    # Extract event info.
    # Event:   2021-03-04T13:27:36.464000Z | -37.563, +179.444 | 7.3 mww | manual
    _event = str(in_catalog)
    _event = _event.split('\n')[0]
    _event = _event.strip()
    items = _event.split('|')
    _date_time = items[0].split()[1].replace('Z', '')
    _datetime = datetime.strptime(_date_time, '%Y-%m-%dT%H:%M:%S.%f')
    event['datetime'] = _datetime

    # Watch for empty descriptions.
    if len(in_catalog.event_descriptions) > 0:
        event['description'] = in_catalog.event_descriptions[0]

        if 'Engdahl' not in event['description'].type:
            for desc in in_catalog.event_descriptions:
                if 'Engdahl' in desc.type:
                    event['description'] = desc
                    if verbose:
                        print(f'[INFO] Event description set to {in_catalog.desc.type} description.',
                              flush=True, file=log_file)
                    # break
    else:

        event['description'] = '-'
        if verbose:
            print(f'[WARN] Missing event description.',
                  flush=True, file=log_file)

    # Event origin.
    event['origin'] = in_catalog.origins[0]
    for orig in in_catalog.origins:
        if orig.resource_id.id == in_catalog.preferred_origin_id:
            event['origin'] = orig
            if verbose:
                print(f'[INFO] Origin set to: {in_catalog.preferred_origin_id}', flush=True, file=log_file)
            break

    # Event magnitude.
    event['magnitude'] = in_catalog.magnitudes[0]
    for mag in in_catalog.magnitudes:
        if mag.resource_id.id == in_catalog.preferred_magnitude_id:
            event['magnitude'] = mag
            if verbose:
                print(f'[INFO] Magnitude set to: {in_catalog.preferred_magnitude_id}', flush=True, file=log_file)
            break

    # Event faocal mechanism.
    event['focal_mechanism'] = None
    event['mt'] = None
    for fm in in_catalog.focal_mechanisms:
        if fm.resource_id.id == in_catalog.preferred_focal_mechanism_id:
            event['focal_mechanism'] = fm
            if verbose:
                print(f'[INFO] Focal Mechanism set to: {in_catalog.preferred_focal_mechanism_id}', flush=True,
                      file=log_file)
            break

    if output:
        print('\n')
        for key in event.keys():
            print(event[key])

    return event
    

def check_date_format(date_string, logfile, time_string=''):
    """Check date_string format and convert it to datetime depending on time."""
    if date_string is None:
        return None, None
    date = date_string.strip()
    while '  ' in date:
        date = date.replace('  ', ' ')
    date = date.replace(' ', 'T')

    time = time_string.strip()
    if 'T' not in date and time:
        date_time = f'{date}T{time}'
    elif 'T' in date and time:
        _date_only = (date.split('T'))[0]
        date = f'{_date_only}T{time}'
        date_time = datetime.strptime(date, '%Y-%m-%dT%H:%M:%S')
        return date, date_time
    else:
        date_time = date

    if 'T' in date_time:
        try:
            date_time = datetime.strptime(date_time, '%Y-%m-%dT%H:%M:%S')
        except Exception as ex:
            print(f'[ERR] Failed to convert "{date_string}" to datetime\n{ex}', flush=True, file=logfile)
            return date_time, None
    else:
        try:
            date_time = datetime.strptime(date_time, '%Y-%m-%d')
        except Exception as ex:
            print(f'[ERR] Failed to convert "{date_string}" to datetime\n{ex}', flush=True, file=logfile)
            return date_time, None
    return date, date_time


def associate(origin1, origin2, association_factor):
    """Check to event origins for possible association."""
    # Origin difference parameters.
    del_lat = abs(origin1['latitude'] - origin2['latitude'])
    del_lon = abs(origin1['longitude'] - origin2['longitude'])
    del_t = abs(origin1['time'] - origin2['time'])

    # Check for association.
    if del_lat <= association_threshold['latitude'] * association_factor and \
            del_lon <= association_threshold['longitude'] * association_factor\
            and del_t <= association_threshold['seconds'] * association_factor:
        return True

    return False


def get_association_index(origin1, orig1_auth, origin2, orig2_auth, log_file, association_factor=1.0):
    """ Check two origin lists for association and return association indexes"""

    assoc_index = {orig1_auth: [], orig2_auth: []}
    if origin1 is None or origin2 is None:
        return assoc_index
    for ind1, orig1 in enumerate(origin1):
        for ind2, orig2 in enumerate(origin2):
            assoc = associate(orig1, orig2, association_factor)

            if assoc:
                if ind1 in assoc_index[orig1_auth]:
                    print(f'[ERR] {orig1_auth} index {ind1} is already associated', flush=True, file=log_file)
                    continue
                if ind2 in assoc_index[orig2_auth]:
                    print(f'[ERR] {orig2_auth} index {ind1} is already associated', flush=True, file=log_file)
                    continue
                assoc_index[orig1_auth].append(ind1)
                assoc_index[orig2_auth].append(ind2)
                break
    return assoc_index


def get_fdsn_events(fdsn_dc, start_str, end_str, lat_limits, lon_limits, min_mag, log_file, max_mag=10, min_depth=None,
                    max_depth=None, center_lat=None, center_lon=None, radius=None, return_url=False):
    """get query-based FDSN events"""

    # Historic seismicity requests may have a start_str set to None
    if start_str is not None:
        if 'T' in start_str:
            start_str, s_date = check_date_format(start_str, log_file)
        else:
            start_str, s_date = check_date_format(start_str, log_file, time_string='00:00:00')

    if 'T' in end_str:
        end_str, e_date = check_date_format(end_str, log_file)
    else:
         end_str, e_date = check_date_format(end_str, log_file, time_string='00:00:00')

    if radius is None:
        _location = f'minlatitude={lat_limits[0]}&maxlatitude={lat_limits[1]}&minlongitude={lon_limits[0]}' \
            f'&maxlongitude={lon_limits[1]}'
    else:
        _location = f'lat={center_lat}&lon={center_lon}&minradius=0&maxradius={radius}'

    _depth = ''
    if min_depth is not None:
        _depth = f'&mindepth={min_depth}'
    if max_depth is not None:
        _depth = f'{_depth}&maxdepth={max_depth}'

    if start_str is None:
        _start = f'&starttime=1976-01-01'
    else:
        _start = f'&starttime={start_str}'

    fdsn_query = (f'{_location}{_start}{_depth}'
                  f'&endtime={end_str}'
                  f'&minmagnitude={min_mag}&maxmagnitude={max_mag}'
                  f'&format=xml&includeallorigins=true&includeallmagnitudes=true')

    url = f'{event_service_url[fdsn_dc]}{fdsn_query}'
    print(f'[INFO] Requesting: {url}', flush=True, file=log_file)
    try:
       catalog = read_events(url)
    except Exception as er:
        print(f'[ERR] Request failed\n{er}', flush=True, file=log_file)
        sys.exit(1)
    
    fdsn_description = list()
    fdsn_origin = list()
    fdsn_magnitude = list()
    fdsn_focal_mechanism = list()
    fdsn_id = list()
    fdsn_mt = list()

    for cat in catalog:
        this_event = get_catalog_event(cat, log_file)

        # Prevent crashes when event description is missing.
        if 'description' in this_event:
            fdsn_description.append(this_event['description'])
        else:
            fdsn_description.append(this_event['-'])

        fdsn_origin.append(this_event['origin'])
        fdsn_magnitude.append(this_event['magnitude'])
        fdsn_focal_mechanism.append(this_event['focal_mechanism'])
        fdsn_mt.append(this_event['mt'])
        fdsn_id.append(this_event['id'])
    fdsn_event = {'description': fdsn_description, 'origin': fdsn_origin, 'magnitude': fdsn_magnitude,
                  'focal_mechanism': fdsn_focal_mechanism, 'mt': fdsn_mt, 'id': fdsn_id}
    if return_url:
        return fdsn_event, url
    else:
        return fdsn_event


def get_event_by_id(event_id, fdsn_dc, scratch_dir, log_file, gcmt_dc=None):
    """From FDSN get an event based on ID and then get the corresponding GCMT event, if any"""

    fdsn_query = f'eventid={event_id}&nodata=404'
    _url = f'{event_service_url[fdsn_dc]}{fdsn_query}'
    print(f'[INFO] Requesting:\n{_url}', flush=True, file=log_file)
    try:
        catalog = read_events(_url)
    except Exception as er:
        print(f'[ERR] Request failed\n{er}', flush=True, file=log_file)
        sys.exit(1)
    fdsn_event = get_catalog_event(catalog[0], log_file)
    if gcmt_dc is None:
        return fdsn_event

    fdsn_origin = fdsn_event['origin']
    search_start = fdsn_origin.time - timedelta(seconds=association_threshold['seconds'],
                                                minutes=0, hours=0, days=0) * association_threshold['gcmt_factor']
    search_end = fdsn_origin.time + timedelta(seconds=association_threshold['seconds'],
                                              minutes=0, hours=0, days=0) * association_threshold['gcmt_factor']
    s_y, s_m, s_d = search_start.strftime('%Y-%m-%d').split('-')
    e_y, e_m, e_d = search_end.strftime('%Y-%m-%d').split('-')

    # GCMT form is not very sensitive to lat/lon, need to increase the threshold for search.
    # Watch for validity.
    lat_min = fdsn_origin.latitude - association_threshold['latitude'] * association_threshold['gcmt_factor']
    if lat_min < -90.0:
        lat_min = -90.0
    lat_max = fdsn_origin.latitude + association_threshold['latitude'] * association_threshold['gcmt_factor']
    if lat_max > 90.0:
        lat_max = 90.0

    lon_min = fdsn_origin.longitude - association_threshold['longitude'] * association_threshold['gcmt_factor']
    if lon_min < -180.0:
        lon_min = -180.0
    lon_max = fdsn_origin.longitude + association_threshold['longitude'] * association_threshold['gcmt_factor']
    if lon_max > 180.0:
        lon_max = 180.0

    depth_min = (fdsn_origin.depth / 1000.0) - association_threshold['depth'] * association_threshold['gcmt_factor']
    if depth_min < 0:
        depth_min = 0.0
    depth_max = (fdsn_origin.depth / 1000.0) + association_threshold['depth'] * association_threshold['gcmt_factor']

    mag_min = fdsn_event['magnitude'].mag - association_threshold['magnitude']
    mag_max = fdsn_event['magnitude'].mag + association_threshold['magnitude']

    gcmt_query = (f'itype=ymd&yr={s_y}&mo={int(s_m)}&day={int(s_d)}'
                  f'&oyr={e_y}&omo={int(e_m)}&oday={int(e_d)}&jyr=1976&jday=1&ojyr=1976'
                  f'&ojday=1&otype=ymd&nday=1&lmw={mag_min}&umw={mag_max}&lms=0&ums=10&lmb=0&umb=10&llat={lat_min}'
                  f'&ulat={lat_max}&llon={lon_min}&ulon={lon_max}'
                  f'&lhd={depth_min}&uhd={depth_max}&lts=-9999&uts=9999&lpe1=0'
                  f'&upe1=90&lpe2=0&upe2=90&list=4')
    _url = f'{event_service_url[gcmt_dc]}{gcmt_query}'
    start_flag = 0
    gcmt_events = ''
    while start_flag is not None:
        if start_flag == 0:
            these_events, start_flag = retrieve_gcmt_events(_url, log_file)
        else:
            these_events, start_flag = retrieve_gcmt_events(f'{_url}&start={start_flag}', log_file)

        if not gcmt_events:
            gcmt_events = these_events
        else:
            gcmt_events = f'{gcmt_events}\n{these_events}'

    if gcmt_events is None:
        print(f'[WARN] No GCMT events found!', flush=True, file=log_file)
        return fdsn_event, None
    else:
        tf_name = os.path.join(scratch_dir, f'{str(uuid.uuid4())}.txt')
        with open(tf_name, 'w') as fp:
            fp.write(gcmt_events)
        fp.close()
        catalog = obspy.read_events(tf_name)
        os.remove(tf_name)

        if catalog:
            for cat in catalog:
                _gcmt_event = get_catalog_event(cat, log_file)

                # Make sure event time is acceptable.
                if search_start <= _gcmt_event['datetime'] <= search_end:
                    _tensor = _gcmt_event['focal_mechanism'].moment_tensor.tensor
                    gcmt_event = _gcmt_event
                    gcmt_event['mt'] = [_tensor.m_rr, _tensor.m_tt, _tensor.m_pp,
                                        _tensor.m_rt, _tensor.m_rp, _tensor.m_tp]
                    break
                else:
                    gcmt_event = None
            if gcmt_event is None:
                print(f'[WARN] No GCMT events returned!', flush=True, file=log_file)
                return fdsn_event, None
        else:
            print(f'[WARN] No GCMT events returned!', flush=True, file=log_file)
            return fdsn_event, None

    gcmt_query = gcmt_query.replace('list=4', 'list=2')
    gcmt_id = gcmt_event['id']

    _url = f'{event_service_url[gcmt_dc]}{gcmt_query}'
    print(f'[INFO] Requesting: {_url}', flush=True, file=log_file)
    strikes = retrieve_gcmt_events(_url, log_file, eid=gcmt_id)
    print(f'[INFO] Strikes are {strikes}', flush=True, file=log_file)
    gcmt_event['strikes'] = strikes

    return fdsn_event, gcmt_event


def get_gcmt_events(start_str, end_str, lat_limits, lon_limits, min_mag, log_file, max_depth=1000, min_depth=0,
                    gcmt_dc='GCMT', scratch_dir=None, return_url=False):
    """get GCMT events"""

    start_str = start_str.strip()
    start_str.replace(' ', 'T')
    if 'T' in start_str:
        start_str = start_str.split('T')[0]

    end_str = end_str.strip()
    end_str.replace(' ', 'T')
    if 'T' in end_str:
        end_str = end_str.split('T')[0]
    s_y, s_m, s_d = start_str.split('-')
    e_y, e_m, e_d = end_str.split('-')

    gcmt_description = list()
    gcmt_origin = list()
    gcmt_magnitude = list()
    gcmt_focal_mechanism = list()
    gcmt_mt = list()
    gcmt_id = list()

    lat_min = lat_limits[0]
    if lat_min < -90.0:
        lat_min = -90.0

    lat_max = lat_limits[1]
    if lat_max > 90.0:
        lat_max = 90.0

    lon_min = lon_limits[0]
    if lon_min < -180.0:
        lon_min = -180.0

    lon_max = lon_limits[1]
    if lon_max > 180.0:
        lon_max = 180.0

    gcmt_query = (f'itype=ymd&yr={s_y}&mo={int(s_m)}&day={int(s_d)}'
                  f'&oyr={e_y}&omo={int(e_m)}&oday={int(e_d)}&jyr=1976&jday=1&ojyr=1976'
                  f'&ojday=1&otype=ymd&nday=1&lmw={min_mag}&umw=10&lms=0&ums=10&lmb=0&umb=10&llat={lat_min}'
                  f'&ulat={lat_max}&llon={lon_min}&ulon={lon_max}'
                  f'&lhd={min_depth}&uhd={max_depth}&lts=-9999&uts=9999&lpe1=0'
                  f'&upe1=90&lpe2=0&upe2=90&list=4')
    
    url = f'{event_service_url[gcmt_dc]}{gcmt_query}'
    start_flag = 0
    gcmt_events = None
    while start_flag is not None:
        if start_flag == 0:
            these_events, start_flag = retrieve_gcmt_events(url, log_file)
        else:
            these_events, start_flag = retrieve_gcmt_events(f'{url}&start={start_flag}', log_file)
    
        if gcmt_events is None:
            gcmt_events = these_events
        else:
            gcmt_events = f'{gcmt_events}\n{these_events}'
    
    if gcmt_events is None:
        print(f'[WARN] No GCMT events found!', flush=True, file=log_file)
        gcmt_event = None
    else:
        tf_name = os.path.join(scratch_dir, f'{str(uuid.uuid4())}.txt')
        # print(f'[INFO] TF:{tf_name}')
        with open(tf_name, 'w') as fp:
            fp.write(gcmt_events)
        fp.close()
        catalog = obspy.read_events(tf_name)
        os.remove(tf_name)

        for cat in catalog:
            this_event = get_catalog_event(cat, log_file)
            gcmt_id.append(cat.resource_id.id)
            gcmt_description.append(this_event['description'])
            gcmt_origin.append(this_event['origin'])
            gcmt_magnitude.append(this_event['magnitude'])
            gcmt_focal_mechanism.append(this_event['focal_mechanism'])
            tensor = this_event['focal_mechanism'].moment_tensor.tensor
            gcmt_mt.append([tensor.m_rr, tensor.m_tt, tensor.m_pp, tensor.m_rt, tensor.m_rp, tensor.m_tp])
        gcmt_event = {'description': gcmt_description, 'origin': gcmt_origin, 'magnitude': gcmt_magnitude,
                      'focal_mechanism': gcmt_focal_mechanism, 'mt': gcmt_mt, 'id': gcmt_id}

    if return_url:
        return gcmt_event, url
    else:
        return gcmt_event


def get_iris_id(fdsn_id, dc):
    """Using an existing FDSN event, get the ccrresponding IRIS event ID """

    search_event = get_event_by_id(fdsn_id, dc, None, None)


    # Now we go for the complete solution based on the event information.
    this_event = dict()
    this_event['latitude'] = search_event['origin'].latitude
    this_event['longitude'] = search_event['origin'].longitude
    this_event['magnitude'] = search_event['magnitude'].mag
    # QuakeML depths are in meters.
    this_event['depth'] = search_event['origin'].depth / 1000.0
    this_event['datetime'] = search_event['origin'].time
    this_event['date_time'] = search_event['origin'].time.strftime('%Y-%m-%dT%H:%M:%S')
    this_event['name'] = search_event['description'].text
    this_event['mt'] = search_event['mt']

    start_str = (this_event['datetime'] - association_threshold['seconds']).strftime('%Y-%m-%dT%H:%M:%S')
    end_str = (this_event['datetime'] + association_threshold['seconds']).strftime('%Y-%m-%dT%H:%M:%S')
    lat_min = this_event['latitude'] - association_threshold['latitude']
    lat_max = this_event['latitude'] + association_threshold['latitude']
    lon_min = this_event['longitude'] - association_threshold['longitude']
    lon_max = this_event['longitude'] + association_threshold['longitude']
    min_mag = this_event['magnitude'] - association_threshold['magnitude']
    max_mag = this_event['magnitude'] + association_threshold['magnitude']
    min_depth = this_event['depth'] - association_threshold['depth']
    max_depth = this_event['depth'] + association_threshold['depth']
    iris_event = get_fdsn_events('IRIS', start_str, end_str, (lat_min, lat_max), (lon_min, lon_max),
                                 min_mag, max_mag=max_mag, min_depth=min_depth, max_depth=max_depth,
                                 center_lat=None, center_lon=None, radius=None)
    if iris_event['id'][0] is None:
        iris_id = None
    elif 'eventid=' not in iris_event['id'][0] :
        iris_id = None
    else:
        iris_id = (iris_event['id'][0].split('eventid='))[1].strip()

    if iris_id is None:
        print(f'[ERR] Was not able to find IRIS ID for event ID {fdsn_id} from {dc} data center!')

    return iris_id


