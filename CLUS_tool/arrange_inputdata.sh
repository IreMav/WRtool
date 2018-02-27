#!/bin/bash

# Information required by the WRtool:

#-------------------------------about data-------------------------------------------
# Write only letters or numbers, no punctuation marks!
# If you want to leave the field empty write 'no' 
varname=zg           #variable name in the file
level=500            #level to select
freq=day             #data frequency ('day','mon',year',...)
filterwnd=5          #running mean filter time window
varunits=m           #variable units
model=ECEARTH31      #model name ECEARTH31 NCEPNCAR ERAInterim
institute=CNR        #institute name
kind=base            #base: baseline, stoc: stochastic physics, obs: observations
res=T255             #T255 144x73 N216L85
numens=3            #total number of ensemble members (put numens=1 and enstoselect=no)
                     #if you have only one realization)
#If you have more than one ensemble member set enstoselect:
enstoselect=no         #single ensemble member to analyse (starting from 0!)
                     #leave it empty (enstoselect=no) if you do not want to select
                     #any particular ensemble member
season=DJF           #seasonal average
area=EAT             #regional average (examples:'EAT':Euro-Atlantic
                     #                           'PNA': Pacific North American
                     #                           'NH': Northern Hemisphere)
syr=1990
eyr=1992

# Previous of the weaher regime analysis,
# data should be remap r144x73 and the period of interest should be selected
# cdo examples below:

FINAL_PATHtmp=/group_workspaces/jasmin2/primavera1/tools/CNR_WeatherRegimes/CLUS_tool/DATAremap/tmp
FINAL_PATH=/group_workspaces/jasmin2/primavera1/tools/CNR_WeatherRegimes/CLUS_tool/DATAremap

if [ $model = ECEARTH31 ]; then
    echo $model
    ####EC-EARTH3.1
    lastens=2
    for ens in $(seq 0 ${lastens}); do
        ORIGINAL_PATH="/group_workspaces/jasmin2/primavera1/WP2/ATM/EC-EARTH3.1/HIST_${res}/${freq}/en${ens}/${varname}"
        for year in $(seq ${syr} ${eyr}); do
            ORIGINAL_FILENAME="${varname}_A${freq}_EC-EARTH31_${res}_${institute}_${ens}_${year}"
            exps_data="${ORIGINAL_PATH}/${ORIGINAL_FILENAME}.nc"
            echo "-------------------------------------------------"
            echo $exps_data
            echo "-------------------------------------------------"
            cdo -L remapcon,r144x73 -sellevel,50000 $exps_data ${FINAL_PATHtmp}/remap_${ORIGINAL_FILENAME}.nc
        done
        ncrcat $(find $FINAL_PATHtmp -type f -name "remap_${varname}_A${freq}_EC-EARTH31_${res}_${institute}_${ens}*.nc") "${FINAL_PATH}/${varname}_A${freq}_EC-EARTH31_${res}_${institute}_${ens}_${syr}-${eyr}.nc"
    done
    ls ${FINAL_PATH}
    rm -rf ${FINAL_PATHtmp}/*.nc
    
elif [ $model = METOFFICE ]; then
    echo $model
    ####METOFFICE
    lastens=3
    for ens in $(seq 1 ${lastens}); do
        ORIGINAL_PATH="/group_workspaces/jasmin2/primavera1/WP2/ATM/METOFFICE/HadGEM3-GA3/${res}/present_day/${freq}/geopotential_height/member${ens}"
        for year in $(seq ${syr} ${eyr}); do
            ORIGINAL_FILENAME="${varname}_${freq}_MetUM_GA3_${res}_r${ens}i1p1_${year}"
            exps_data="${ORIGINAL_PATH}/${ORIGINAL_FILENAME}.nc"
            echo "-------------------------------------------------"
            echo $exps_data
            echo "-------------------------------------------------"
            cdo -L remapcon,r144x73 -sellevel,500 $exps_data ${FINAL_PATHtmp}/remap_${ORIGINAL_FILENAME}.nc
        done
        ncrcat $(find $FINAL_PATHtmp -type f -name "remap_${varname}_${freq}_MetUM_GA3_${res}_r${ens}i1p1*.nc") "${FINAL_PATH}/${varname}_${freq}_MetUM_GA3_${res}_r${ens}i1p1_${syr}-${eyr}.nc"
    done
    ls ${FINAL_PATH}
    rm -rf ${FINAL_PATHtmp}/*.nc
fi

