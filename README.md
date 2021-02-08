 Incorporated Research Institutions for Seismology (IRIS)\
 Data Management Center (DMC)\
 Data Products Team\
 Aftershocks Data Product

 2020-09-16

------------------------------------------------------------------------------------------------------------------------

 DESCRIPTION:

Following all global earthquakes with magnitudes greater than 7, IRIS DMC Aftershocks data product 
(http://ds.iris.edu/ds/products/aftershocks/) automatically generates a sequence animations and figures and makes them 
 available to public via DMC's Searchable Product Depository, SPUD (http://ds.iris.edu/spud/aftershock). Starting on August 1, 
2020, events from the FDSN (USGS) web service (https://earthquake.usgs.gov/fdsnws/event/1/) and the GCMT catalog 
(https://www.ldeo.columbia.edu/~gcmt/GCMT_Website_Additions/www.globalcmt.org.ghost.html) are used for this purpose. 
For each qualified event, aftershocks figures and animations are updated hourly for 10 days post event. 


This Python package contains the code behind the creation of Aftershocks data product's plots and animations. The main 
Python code in this package (_aftershock_fdsn_maps.py_) can be configured via its parameter file 
(_aftershock\_fdsn\_maps\_param.py_) or through the command line arguments. Currently, the package is configured to use 
the FDSN event service from USGS to get the base event information. Then it queries the Global CMT (GCMT) catalog and 
will attempt to associate each GCMT event with the corresponding FDSN event. If association is successful, the FDSN 
event is replaced with the corresponding GCMT event. User may select to also plot the unassociated GCMTs (if any). 
Changing the FDSN event service provider (https://www.fdsn.org/webservices/) should be as simple as changing the 
parameters. However, some data centers may have slightly different configuration and therefore change the service 
provider may require further parameter tuning.

 CONTENT:

This package contains the following files:

     src/
       aftershock_fdsn_maps.py
           This is the main Python code behind the production of plots and animations. Calling the code with -h  
           displays a list of other options available to tune plot and animation production. It also provides test 
           examples to run.
     
     param/
       aftershock_fdsn_maps_param.py
           A Python file that contains all Aftershocks data product parameters. You may modify this file to 
           customize plot and animation production. All parameter definitions in this file must follow Python 
           rules. Each parameter group in this file is commented for clarification.
     
     lib/
       - event_fdsn_lib.py
           A Python FDSN/GCMT event library.

       - aftershocks_fdsn_lib.py
           A Python utility library used by the main script.

    CHANGES.txt
       A text file containing history of changes to this package.

    INSTALL.txt
       The installation notes

    README.md
       This file


 INSTALLATION:

    see the INSTALL.txt file


USAGE:
   
    aftershock_fdsn_maps.py (V.2020.266):

    This is the Python code behind the IRIS DMC's Aftershocks Data product. It is capable of producing plots and 
    animations that are part of the IRIS DMC's Aftershocks Data product (http://ds.iris.edu/spud/aftershock).

    The code can be configured via its parameter file "aftershock_maps_param.par" or via the command 
    line arguments. Currently parameters are optimized for use with the Mercator map projection, NEIC/USGS FDSN 
    event services and GCMT event catalog. Changes in projection, resolution or data provider may require 
    parameter tuning and/or code update.

    command line options:
	    -h --help		   this message
	    -v --verbose	   run in verbose mode
	    -e --eid [event ID]	   process this USGS event ID as the mainshock
	    -b --before [integer]  number of days before the event to display 
	    -a --after [integer]   number of days after the event to display 
	    -m --minmag [float]	   minimum event magnitude to consider
	    -l --label [string]	   include this label in the file name
        -p --plot [1,...]	   plots to create. Use one or more (comma separated) indices from the list of 
                               possible plots (see the "map_description" variable in the param file). 

	    	index		description
		=====		===========
	    	0		Background seismicity within 100 km depth
		1		Seismicity  within 10 days of mainshock
	    	2		[requires GCMT solution] Aftershock distance along strike (nodal plane) vs day
		3		[requires GCMT solution] Aftershock distance along strike (nodal plane) vs depth
	    	4		Location of events within 10 days of mainshock
		5		heatmap of the seismicity  within 10 days of mainshock
	    	6		Animation of the seismicity  within 10 days of mainshock
		7		Animation of the seismicity heatmap within 10 days of mainshock

        -r --refid [event ID]       process this FDSN event ID as a reference event
	    -s --scale [float]          ETOPO uses scale to reduce map resolution for maps that are zoomed in. 
	                                The default scale is automatically set between 1.2 and 2.6. The larger scale values 
	                                require more memory
	    -x --xdate	                if -x is present, the x-axis will have date labels rather than day labels
        -T --title [double-quoted]	use this plot title instead of the default title

    NOTE: either eid or refid should be provided

    Example:   
    
        From: https://earthquake.usgs.gov/fdsnws/event/1/query?format=text&starttime=2020-07-22&endtime=2020-07-23&minmag=6.8&nodata=404

        we obtain the event ID and then run:
             aftershock_fdsn_maps.py -v -e us7000asvb
        to create aftershock plots for event us7000asvb


For additional configuration, please see the parameter file (aftershock_maps_param.py)


EXAMPLES:

    To make plots for an event of interest, you need to go to the USGS earthquake catalog search page 
    (https://earthquake.usgs.gov/fdsnws/event/1/)
    and query for the region of interest and find the ID for that event.

    Here is a sample URL to search for Alaska events with magnitude 6 and larger starting July 19, 2020 with 
    the output format as text.

    https://earthquake.usgs.gov/fdsnws/event/1/query?lat=54.84&lon=-159.27&minradius=0&maxradius=2.698
    &starttime=2020-07-19T06:12:44&mindepth=0&maxdepth=136.81&endtime=2020-08-01T06:12:44
    &format=text&minmagnitude=6&maxmagnitude=10&nodata=404

    We select the main event of M 7.8 on July 22 with event ID us7000asvb in the first column:

    - To create a Seismicity map within 10 days of mainshock (index 1 in the usage above) run:
       src/gmv_generalized.py -e us7000asvb -p 1 -x 
          . -e us7000asvb to specify the event ID, 
          . -p to plot the seismicity map,
          . -x to select date label for x-axis

       The code should log its actions and at the end provide the path to the created image under the image 
       directory. If you want the code to log to a file, then set the "log_to_screen" parameter to False in 
       the parameter file.

    - To create more than one plot, run:
       src/gmv_generalized.py -e us7000asvb -p 1,4 -x
      This command also plots a location map (index 4) in addition the the seismicity may (index 1). 
      
    - Omit the -p option and have the script create all the indices listed by "map_list = [1, 2, 3, 4, 5, 6, 7]"
      in the parameter file
    

		
SELECTING PARAMETERS:

 This package produces event-based plots and animations and therefore you need to start with an event. The main 
  event could be identified as either a "mainshock" (by providing its ID via _-e_ or _--eid_ option) or as a 
  "reference event" (by providing its ID via _-r_ or _--refid_ option). Plots for both type of events are 
  identical except for the  title and some labeling. Use the reference event when you want to use the code to 
  create seismic activity plots for an area (for example, swarms in an area). 
  
  You should obtain the Event ID or the reference ID from the FDSN data center that the code queries. Currently 
  the package uses USGS event service and therefore, you need to go to the USGS earthquake catalog search page 
  (https://earthquake.usgs.gov/fdsnws/event/1/) and query for the event of interest and find the ID for that 
  event (see the EXAMPLES section above).
 
All other parameters are optional. However, you may tune some of the parameters via the call arguments:
 - change the number of days before (_-b_ option) the event for which the seismicity should be included. 
 Use "none" if you want to get the historic seismicity
 
 - change the number of days after (_-a_ option) the event for which the seismicity should be included.
 
 - use _-m_ to set the minimum magnitude for the events that are included. Currently the minimum magnitude 
 depends on the main or reference event's magnitude. The minimum magnitude and other plot parameters are set 
 by the _set\_map\_parameters_ function in the aftershocks library file.
 
 - to label plot files with a customized tag, use the _-l_ option. This tag is included in the file names.
 
 - to plot specific map(s) and/or animation(s), use the _-p_ option and provide index of one or more plots 
 to creates (see the _map\_description_ in the parameter file for a complete list of plot options 
 (also see below). 
 
 		    0		Background seismicity within 100 km depth
		    1		Seismicity  within 10 days of mainshock
		    2		Aftershock distance along strike (nodal plane) vs day
		    3		Aftershock distance along strike (nodal plane) vs depth
		    4		Location of events within 10 days of mainshock
		    5		heatmap of the seismicity  within 10 days of mainshock
		    6		Animation of the seismicity  within 10 days of mainshock
		    7		Animation of the seismicity heatmap within 10 days of mainshock

 - omit the _-p_ option and have the script create maps/animations for all indices listed by  
 _map\_list = [1, 2, 3, 4, 5, 6, 7]_ in the parameter file.
 		    
 - use  _-s_ (_--scale_) option to control the resolution of the map's backgrounds. ETOPO uses the "scale" to 
 reduce map resolution for maps that are zoomed in. The default is selected automatically between 1.2 and 2.6 
 with larger values require more memory. 
 - by default, the day axis of the plots are labeled by the number of days to or from the event origin. 
 Use -x to change the labeling to date.
 - to customize the plot labels, use the _-T_ option and place the plot title inside double quotes


CITATION:

To cite the use of this software reference:

Hutko, A. R., M. Bahavar, C. Trabant, R. T. Weekly, M. Van Fossen, T. Ahern (2017), Data Products at the IRIS‐DMC: \
Growth and Usage, Seismological Research Letters, 88, no. 3, https://doi.org/10.1785/0220160190.\
\
Or cite the following DOI:

    doi:10.17611/dp/as.code.1


 HISTORY
- 2020-09-16 V.2020.260 FDSN R.2.1 public release
- 2020-08-22 V.2020.236 FDSN support
- 2020-08-01 V.2020.214 R2 release.
- 2014-12-17 R1, development and initial release.

 
 COMMENTS/QUESTIONS:

    Please contact manoch@iris.washington.edu


