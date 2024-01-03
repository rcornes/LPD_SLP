## --------------------------------------------------------------------------------
## Functions to generate geostrophic wind speed/direction from sea-level pressure triangle
##
## This follows the method described by Wang et al. (2009) and Krueger et al. (2009)
##
## Written by R. Cornes, National Oceanography Centre, Southampton (ricorne@noc.ac.uk)
##
## References
## Krueger, O., F. Feser, and R. Weisse, 2019: Northeast Atlantic
## Storm Activity and Its Uncertainty from the Late Nineteenth to the
## Twenty-First Century. J. Climate, 32, 1919–1931,
## https://doi.org/10.1175/JCLI-D-18-0505.1.
##
## Wang, X.L., Zwiers, F.W., Swail, V.R. et al. Trends and variability
## of storminess in the Northeast Atlantic region, 1874–2007. Clim Dyn
## 33, 1179–1195 (2009). https://doi.org/10.1007/s00382-008-0504-5
## --------------------------------------------------------------------------------

geowind.UV <- function(slp,lat.deg,lon.deg){
    
    stopifnot( length(slp)==3, length(lat.deg)==3, length(lon.deg)==3 )

    ## Convert hPa to Pa
    slp <- slp * 100

    ## Convert degrees to radians
    lon <- lon.deg*pi/180
    lat <- lat.deg*pi/180

    ## Define Constants
    R <- 6378100.                       #Radius of the earth (metres)
    rho <- 1.25                         #Density of air
    omega <- 2*pi/(24*3600)             #Rotation rate of earth (radians per second)
    f <- mean( 2*omega*sin(lat) )       #Coriolis parmeter (averaged)
    
    ## Solve the three simultaneous equations
    xcoef <- R * lon * cos(lat)
    ycoef <- R * lat

    A <- rbind(c(xcoef[1], ycoef[1], 1),
               c(xcoef[2], ycoef[2], 1),
               c(xcoef[3], ycoef[3], 1) )
    b <- matrix(data=slp, nrow=3, ncol=1, byrow=FALSE)

    coefs <- solve(A, b)

    ## Calculate u and v vectors
    Ug <- -coefs[2]/(rho*f)
    Vg <- coefs[1]/(rho*f)

    ## Return U and V parameters
    c("U"=Ug, "V"=Vg,"a"=coefs[1],"b"=coefs[2],"c"=coefs[3])
}

geowind.speed <- function(U,V) (U^2 + V^2)^0.5

geowind.dir <- function(U,V) (atan2(U,V) * 360/2/pi) + 180
