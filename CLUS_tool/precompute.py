#*********************************
#            precompute          *
#*********************************

# Standard packages
import sys
import os
import numpy as np
import warnings
warnings.filterwarnings("ignore")

print('*********************************************************************************')
print('Running {0}'.format(sys.argv[0]))
print('*********************************************************************************')
dir_WRtool    = sys.argv[1]  # WRTOOL directory
dir_OUTPUT    = sys.argv[2]  # OUTPUT directory
name_outputs  = sys.argv[3]  # name of the outputs
varname       = sys.argv[4]  # variable name
level         = int(sys.argv[5])  # level to select
wnd           = int(sys.argv[6])  # running mean filter time window
numens        = int(sys.argv[7])  # total number of ensemble members
season        = sys.argv[8]  # season
area          = sys.argv[9]  # area
filenames     = sys.argv[10:len(sys.argv)+1]  # file names (absolute paths)

# User-defined libraries
sys.path.insert(0,dir_WRtool+'WRtool/')
from readsavencfield import read4Dncfield, save3Dncfield
from manipulate import sel_season, anomalies, sel_area

# OUTPUT DIRECTORY
OUTPUTdir=dir_OUTPUT+'OUTPUT/'
OUTfig=OUTPUTdir+'OUTfig/'
OUTnc=OUTPUTdir+'OUTnc/'
OUTtxt=OUTPUTdir+'OUTtxt/'
OUTpy=OUTPUTdir+'OUTpy/'
if not os.path.exists(OUTPUTdir):
    os.mkdir(OUTPUTdir)
    print('The output directory {0} is created'.format(OUTPUTdir))
    os.mkdir(OUTfig)
    os.mkdir(OUTnc)
    os.mkdir(OUTtxt)
    os.mkdir(OUTpy)
else:
    print('The output directory {0} already exists'.format(OUTPUTdir))

#____________Reading the netCDF file of 4Dfield, for one or more ensemble members
for ens in range(numens):
    ifile=filenames[ens]
    if numens>1:
        print('\n/////////////////////////////////////\nENSEMBLE MEMBER %s' %ens)
    elif numens==1:
        print('\n/////////////////////////////////////')
    var, level, lat, lon, dates, time_units, var_units, time_cal = read4Dncfield(ifile,lvel=level)
    var_season,dates_season,var_filt = sel_season(var,dates,season,wnd) # Select data for that season and average data for that season after a running mean of window=wnd
    var_anom = anomalies(var_filt,var_season, dates_season,season)#,freq)                        # Compute anomalies by removing the time-mean for each time (day, month,...)
    namefile='glob_anomaly_{0}_{1}.nc'.format(name_outputs,ens) #like anomaly_zg500_day_ECEARTH31_base_T255_DJF_EAT_0
    ofile=OUTnc+namefile
    save3Dncfield(lat,lon,var_anom,varname,var_units,dates_season,time_units,time_cal,ofile)
    var_area, lat_area, lon_area = sel_area(lat,lon,var_anom,area)      # Selecting only [latS-latN, lonW-lonE] box region defineds by area
    namefile='area_anomaly_{0}_{1}.nc'.format(name_outputs,ens) #like anomaly_zg500_day_ECEARTH31_base_T255_DJF_EAT_0
    ofile=OUTnc+namefile
    save3Dncfield(lat_area,lon_area,var_area,varname,var_units,dates_season,time_units,time_cal,ofile)

print('\nINPUT DATA:')
print('Variable is {0}{1} ({2})'.format(varname,level,var_units))
print('Number of ensemble members is {0}'.format(numens))
print('Input data are read and anomaly fields for each ensemble member are computed and saved in {0}'.format(OUTnc))

print('\n*********************************************************************************')
print('END {0}'.format(sys.argv[0]))
print('*********************************************************************************')
