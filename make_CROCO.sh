#!/bin/bash

# Select your year range from CROCO and the script will create the monthly means for each year
# and then concatenate the outputs into one year as well as create a climatology from your
# entire model period.


# Refer to your crocotools_param.m settings for refdate

REFDATE=1990
START=2002
END=2017

# Start the loop to avergae the data from the avg outputs into monthly form for each year
for year in `seq $START $END`
do
	for month in `seq 1 12`
	do	
	    cdo monavg -selname,temp  croco_avg_Y${year}M${month}.nc tmp$month.nc  
    	done	
    cdo cat tmp{1..12}.nc tmp${year}.nc
    cdo -setcalendar,standard -setreftime,${REFDATE}-01-01,00:00:00,seconds -settaxis,${year}-01-01,00:00:00,1mon tmp${year}.nc croco_avg_Y${year}.nc
    rm tmp*	
done 

# Concatenate the monthly file together
cdo cat croco_avg_Y20??.nc croco_avg_SST_${START}_${END}.nc

# Create the climatology for the model period
#cdo -f nc ymonmean croco_avg_${START}_${END}.nc croco_clim_${START}_${END}.nc 
