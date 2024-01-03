#!/bin/bash

## --------------------------------------------------------------------------------
## Script to download the 20CR V2/V3 ensemble mean data of U850 and MSLP 
## Zonal mean values and monthly/seasonl means over the N.Atlantic/Europe are saved
##
## Written by Richard Cornes (ricorne@noc.ac.uk), National Oceanography Centre.
## On 2023-07-06
##
## This script uses the CDO software (https://code.mpimet.mpg.de/projects/cdo)
## --------------------------------------------------------------------------------

module load cdo/1.9.5 nco

for y in {2015..2015}; do
    yearfile=UGRD850_${y}_daily.tar
    #wget https://portal.nersc.gov/archive/home/projects/incite11/www/20C_Reanalysis_version_3/everymember_anal_netcdf/daily/UGRD850/${yearfile} .
    
    for i in {001..080}; do

	selfile=${y}/UGRD850.${y}.daily_mem${i}.nc
	
	tar -xf $yearfile $selfile && mv $selfile tempfile.nc
	ncks -O -d lon,300,0 -d lat,105,165 tempfile.nc dummy.nc
	cdo -f nc2 zonmean dummy.nc 20CRv3_${y}_zonmean.nc
	rm -r tempfile.nc dummy.nc $y

    done
done

#cdo mergetime 20CRv3_*_zonmean.nc 20CRv3_zonmean.nc
