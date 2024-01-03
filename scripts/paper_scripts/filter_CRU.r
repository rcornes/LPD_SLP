filter.cru <- function(thalf,tsin,tslow.return=TRUE,tshigh.return=FALSE,nan=TRUE,weight.return=FALSE){
###################################################  
# Translated from Tim Osborn's filter_cru  
# Richard Cornes
# 30th August 2009
###################################################
#  
# Uses the CRU Gaussian weighted filter to separate the low and high frequency
# components of a timeseries.
#
# thalf          : period of oscillation that is reduced in amplitude by 50%
# tsin           : input timeseries (1D only).
# nan            : If this flag is set, missing values in tsin are
#              replaced by the timeseries local mean before filtering, and then
#                  re-replaced by the missing code afterwards (and in the high
#                  and low frequency series too)
# tslow,tshigh   : low and high frequency components
# weight         : optionally return weights used
#
#-----------------------------------------------------------------------------
#
# Define arrays
#
#tssize <- size(tsin)
if (is.vector(tsin)!=1) stop('tsin should be vector')
nt <- length(tsin)
#print(nt)
tslow <- numeric(nt)
tshigh <- numeric(nt)
#
# Compute number of weights required
#
nw <- as.integer(thalf/2.5+0.5)
if ((nw %% 2)==0) nw <- nw+2 else nw <- nw+1
if (nw<=7) nw <- 7
#print(nw)
weight <- numeric(nw)
#
# Compute weights
#
wfactor <- -18/(thalf*thalf)
wroot <- 1/sqrt(2.*pi)
weight[1] <- wroot
wsum <- weight[1]
for (i in 2:nw){
  weight[i] <- wroot*exp(wfactor*(i-1)*(i-1))
  wsum <- wsum+2.*weight[i]
}
weight <- weight/wsum
#
# If required, pad the timeseries with its local mean where values are missing
#
tspad <- tsin
if (nan){
  misslist <- which(is.na(tsin))
  nmiss <- length(misslist)
  if (nmiss>0 & nmiss<nt){
    for (i in 1:nmiss){
      ele1 <- (misslist[i]-nw+1)
      ifelse(ele1<1, ele1 <- 1,ele1 <- ele1)
      ele2 <- (misslist[i]+nw)
      ifelse(ele2>nt,ele2 <- nt,ele2 <- ele2)
      locvals <- tsin[ele1:ele2]
      locmean <- mean(locvals,na.rm=TRUE)
      tspad[misslist[i]] <- locmean
    }
  }
}
#
# Extend ends of timeseries by mean from each end
#
nend <- nw-1
meanst <- sum(tspad[1:nend])/nend
meanen <- sum(tspad[(nt+1-nend):nt])/nend
tspad <- c(rep(meanst,nend),tspad,rep(meanen,nend))
#
# Apply the filter
#
for (i in 1:nt){
  wsum <- weight[1]*tspad[i+nend]
  for (j in 2:nw){
    wsum <- wsum+weight[j]*(tspad[i+nend-j+1]+tspad[i+nend+j-1])
  }
  tslow[i] <- wsum
}
#
# Compute the residual (high-frequency) component
#
tshigh <- tsin-tslow
#
# Insert the missing value if required
#
if (nan){
  if (nmiss>0 & nmiss<nt){
    tslow[misslist] <- NA
    tshigh[misslist] <- NA
}
}

if(tslow.return) ret.list <- list(tslow=tslow)
if(tshigh.return) ret.list <- list(tshigh=tshigh)
if(weight.return) ret.list <- list(weight=weight)
return(ret.list)
}
