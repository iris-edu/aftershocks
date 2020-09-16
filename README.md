 
 Incorporated Research Institutions for Seismology (IRIS)
 Data Management Center (DMC)
 Data Products Team
 Aftershocks Data Product

 2020-09-16

------------------------------------------------------------------------------------------------------------------------

 DESCRIPTION:

This is an IRIS DMC Python code bundle to produce animations and figures of seismicity similar to those of the
Aftershocks data product (http://ds.iris.edu/ds/products/aftershocks/). The main script in this bundle
(src/aftershock_fdsn_maps.py) can be configured via its parameter file (param/aftershock_fdsn_maps_param.py) or  command
line arguments.

Currently bundle is configured to use the FDSN event service from USGS to get the base event information. Then it
queries the Global CMT (GCMT) catalog and will attempt to associate GCMT events with the FDSN event. If an association
is established, the FDSN event is replaced with the GCMT event. User may also select to  plot the remaining unassociated
GCMT events (if any). Changing the FDSN (https://www.fdsn.org/webservices/) event service provider should be as simple
as changing the parameters. However, some data centers may have slightly different configuration and as a result,
switching to a different FDSN data event service provider may require further tuning.

 PYTHON REQUIREMENTS:

       - this bundle has been tested under Python 3.8.3, Anaconda3-2020.07 on macOS 10.14.6, Windows 10.1903,
          and Linux CentOS7.

       - additional required Python module(s) with the minimum version:
             . obspy                1.2.2
             . scipy                1.5.2
             . numpy                1.19.1
             . matplotlib           3.3.1
             . basemap              1.2.2
             . ffmpeg               4.3.1
             . basemap-data-hires   1.2.2

 ETOPO REQUIREMENTS:

       This bundle uses ETOPO1 (a 1 arc-minute global relief model of Earth's surface) to draw topography when needed.
       For this, it uses the ETOPO1 built using GMT 4.3.1 (ETOPO1_Ice_g_gmt4.grd, available from the NOAA website
       https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/grid_registered/netcdf/). The
       ETOPO1 file (~1 GB) should be download from NOAA at above URL and saved under the "assets" directory.
       If you have the file saved in another location, make sure to point the parameter file to that location
       (see the "topo_file" parameter in the parameter file).

 BUNDLE INSTALLATION:

       unpack the bundle
            . "src" directory contains the main code.
            . The parameter file is under the "param" directory.
            . "lib" directory has the library Python files.
            . "assets" empty reserved for the ETOPO1 file (~1 GB), should be downloaded from NOAA (see ETOPO REQUIREMENTS).

 CONFIGURE THE PACKAGE:

        With Python configured, you should be able to run the package examples without further modifications. However:
        - if necessary, update the Python path on the first line of the src/aftershock_fdsn_maps.py
        - if desired, configure the package by updating param/aftershock_fdsn_maps_param.py file
        - For the conda version of Basemap, it needs the "PROJ_LIB" variable to be set so it can find the epsg data.
          You can either set the environment variable "PROJ_LIB" to ""{your anaconda path}/share/proj" or place the
          "os.environ['PROJ_LIB'] = '{your anaconda path}/share/proj'" statement just before importing the Basemap in the
          main code.

 PACKAGE TEST:

    Run src/python aftershock_fdsn_maps.py with "-h" option (assuming you have properly configured your Python
    installation) and it will print a USAGE messages.

    Run  "src/gmv_generalized.py -e us7000asvb -p 1 -x"  to create a Seismicity map within 10 days of the M 7.8
    Alaska event of July 22. The plot will appear under the image directory.






