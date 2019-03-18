
This codebase includes the code for extracting and processing the AIMS manta tow, modelled cyclone
and early bleaching data as well as all analyses and summarizing routines used in the paper.

This codebase includes code for extracting data from AIMS' databases.  Currently, this operation can
only be performed by AIMS staff with adequately credentials. The java applet used to perform the
data extraction dumps the extraction as a flat csv file.  These csv files (along with other files)
are available within a publically accessible data holding hosted by AIMS (link provided *HERE*).
                                                                                                   
The csv data sets (and the locations of where they should be located in the working directory) are:
 - manta.csv (data/primary)
 - bleaching-manta-aesthetics.csv (data/primary)
 - cots.csv (data/primary)
 - 170909 Cyclone wave data from Marji.csv (data/primary)
                                                                                                   
In addition, there are a number of datasets that were provided as Excel sheets (and the locations of
where they should be located in the working directory) are:
  - 2016 Bleaching Response AIMS Site Selection_1998 2002 Repeated Aerial Su....xlsx (data/primary)
                                                                                                   
Some of the mapping routines require specific shapefiles.  The code for the production of these
figures is provided here for transparency and illustrative purposes.  However, the shapefiles are
not publically available.

The file CoralTrends.R acts as an index to all analyses.  This scripts loads a script of helper
functions, loads all required libraries (or indicates which are missing), and runs all other R
scripts. 
