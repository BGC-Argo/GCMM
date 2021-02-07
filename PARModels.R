
dirProgram <- "D:/R-Program/" # your function folder

source(paste(dirProgram,"EstPAR.R",sep=""))


PARModelGCMMsurf <- function(PAR0, Chlasurf, Lon, Lat, Year, Month, Day, Hour, Minute, Alpha=0.98) {

DATE <- paste(Year, "-", Month, "-", Day, sep="")

JDay <- as.double(format(as.Date(DATE), "%j"))

Time <- Hour + Minute/60

z1 <- seq(1, 400, 1)

Chla1 <- z1*0 + Chlasurf

GCres <- EstPAR(JDay, Lon, Lat, Time, Chla1, z1, Alpha=Alpha)

GCPAR0 <- GCres$PAR0

GCPARProfile <- GCres$PARProfile

PARGCMM <- PAR0/GCPAR0*GCPARProfile

return(PARGCMM)

}



PARModelGCMMprof <- function(PAR0, Chlaprof, z, Lon, Lat, Year, Month, Day, Hour, Minute, Alpha=0.98) {

DATE <- paste(Year, "-", Month, "-", Day, sep="")

JDay <- as.double(format(as.Date(DATE), "%j"))

Time <- Hour + Minute/60

z1 <- seq(1, 400, 1)

Chla1 <-  approx(z[!is.na(Chlaprof)], Chlaprof[!is.na(Chlaprof)], z1, rule=2)$y

GCres <- EstPAR(JDay, Lon, Lat, Time, Chla1, z1, Alpha=Alpha)

GCPAR0 <- GCres$PAR0

GCPARProfile <- GCres$PARProfile

PARGCMM <- PAR0/GCPAR0*GCPARProfile

return(PARGCMM)

}




PARModelLee05 <- function(PAR0, a490, bb490, zenith) {

z1 <- seq(1, 400, 1)

Kd_Lee05 <- Lee05(a490, bb490, zenith, z1)

PAR_Lee05 <- PAR0*exp((-1)*Kd_Lee05*z1)

return(PAR_Lee05)

}




Lee05 <- function(a490, bb490, zenith, z) {

KdPARz <- z*NA

X0 <- (-0.057)
X1 <- 0.482
X2 <- 4.221
V0 <- 0.183
V1 <- 0.702
V2 <- (-2.567)
a0 <- 0.090
a1 <- 1.465
a2 <- (-0.667)

K1 <- (X0 + X1*a490^0.5 +X2*bb490)*(1+a0*sin(zenith/180*pi))

K2 <- (V0 + V1*a490 + V2*bb490)*(a1 + a2*cos(zenith/180*pi))

KdPARz <- K1 + K2/(1+z)^0.5

return(KdPARz)

}




