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

cd ../data

for y in {1851..2015}; do
    if(($y < 1981)); then
	prefix=SI
    else
	prefix=MO
    fi
    
    for v in uwnd.${y}.nc; do
	echo $y
	ncks -O -d level,5 -d lon,300,0 -d lat,105,165 https://psl.noaa.gov/thredds/dodsC/Datasets/20thC_ReanV3/Dailies/prs${prefix}/${v} dummy.nc
	cdo -f nc2 zonmean -selvar,uwnd dummy.nc 20CRv3_${y}_zonmean.nc
	rm dummy.nc
	ncks -O -d level,3 -d lon,150,0 -d lat,8,37  https://psl.noaa.gov/thredds/dodsC/Datasets/20thC_ReanV2c/Dailies/pressure/${v} dummy.nc
	cdo -f nc2 zonmean -selvar,uwnd dummy.nc 20CRv2c_${y}_zonmean.nc
	rm dummy.nc
    done
done

cdo mergetime 20CRv3_*_zonmean.nc 20CRv3_zonmean.nc
cdo mergetime 20CRv2c_*_zonmean.nc 20CRv2c_zonmean.nc

# Monthly means of u850 and SLP
ncks -O -d level,5 -d lon,300,60 -d lat,90,180 https://psl.noaa.gov/thredds/dodsC/Datasets/20thC_ReanV3/Monthlies/prsSI-MO/uwnd.mon.mean.nc 20CR_v3_u850_NA.nc
ncks -O -d lon,300,60 -d lat,90,180 https://psl.noaa.gov/thredds/dodsC/Datasets/20thC_ReanV3/Monthlies/miscSI-MO/prmsl.mon.mean.nc 20CR_v3_SLP_NA.nc

# seasonal means
cdo seasmean 20CR_v3_u850_NA.nc 20CR_v3_u850_NA_seas.nc
cdo seasmean 20CR_v3_SLP_NA.nc 20CR_v3_SLP_NA_seas.nc

# Monthly anomalies - seasonal means
cdo ymonmean -selvar,prmsl -selyear,1961/1990 20CR_v3_SLP_NA.nc 20CR_v3_SLP_NA_clim.nc
cdo seasmean -ymonsub -selvar,prmsl 20CR_v3_SLP_NA.nc  20CR_v3_SLP_NA_clim.nc 20CR_v3_SLP_NA_seas_anom.nc 
