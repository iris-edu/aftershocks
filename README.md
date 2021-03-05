 Incorporated Research Institutions for Seismology (IRIS)\
 Data Management Center (DMC)\
 Data Products Team\
 Aftershocks Data Product

 2021-02-09

------------------------------------------------------------------------------------------------------------------------

 DESCRIPTION:

Following all global earthquakes with magnitudes greater than 7, IRIS DMC Aftershocks data product 
(http://ds.iris.edu/ds/products/aftershocks/) automatically generates a sequence of animations and figures and makes them 
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
parameters. However, some data centers may have slightly different configuration and therefore changing the service 
provider may require further parameter tuning.

 CONTENT:

This package contains the following files:

     src/
       aftershock_fdsn_maps.py
           This is the main Python code behind the production of plots and animations. Calling the code with -h  
           option displays a list of other options available to tune plot and animation production. It also 
           provides test examples to run.
     
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
       The package README file 


Visit the Wiki pages (https://github.com/iris-edu/aftershocks/wiki) for more information.


CITATION:

To cite the use of this software reference:

Hutko, A. R., M. Bahavar, C. Trabant, R. T. Weekly, M. Van Fossen, T. Ahern (2017), Data Products at the IRIS‐DMC: \
Growth and Usage, Seismological Research Letters, 88, no. 3, https://doi.org/10.1785/0220160190.\
\
Or cite the following DOI:

    doi:10.17611/dp/as.code.1


 HISTORY
- 2021-03-04 v.2021.063 r2.4 check for lat/lon limits when making GCMT requests.
                                   check for missing description in QuakeML.
- 2021-02-17 v.2021.048 r2.3 updated the usage message
- 2021-02-13 v.2021.044 r2.2 updated heatmap scale label
- 2021-02-09 v.2021.040 FDSN r2.1 public release
- 2020-08-22 v.2020.236 FDSN support
- 2020-08-01 v.2020.214 R2 in production.
- 2014-12-17 r1, development and initial release.

 
 COMMENTS/QUESTIONS:

    Please contact manoch@iris.washington.edu


