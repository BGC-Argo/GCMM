
dirProgram <- "D:/R-Program/" # your function folder

source(paste(dirProgram,"ModelKdMM01.R",sep=""))

require(atmos)

EstPAR <- function(JDay, Lon, Lat, Time, Chla, Depth, Alpha=0.98) {

ResGC <- GreggCarder.f(JDay, Lon, Lat, hr=Time)

planck <- 6.6260755*10^(-34)

lightspeed <- 299792458

Avogadro <- 6.02214086*10^23

Lam <- ResGC$lam

Edabove <- ResGC$Ed

Ed0 <- Edabove*Alpha

WL <- c(400:700)

Sum0 <- 0

for(iwav in 2:301) {

Sum0 <- (Ed0[iwav]*WL[iwav] + Ed0[iwav-1]*WL[iwav-1])/2 + Sum0

}

PAR0 <-  Sum0/planck/lightspeed/10^3/Avogadro 

EdProfile <- matrix(data=NA,nrow=length(Depth), ncol=301)

Kd490 <- 0.0166 + 0.07242*Chla^0.68955

KdSpectra <- matrix(data=NA,nrow=length(Depth), ncol=301)

for(idepth in 1:length(Depth)) {

KdSpectra[idepth,] <- ModelKdMM01(Kd490[idepth], wl=490, range=seq(400,700,1))

}

PARProfile <- 0*Depth

KdPARProfile <- 0*Depth

for (iwav in 1:301) {

IntegratedKd <- NA*Depth

IntegratedKd[1] <- Depth[1]*KdSpectra[1,iwav]

for (i in 2: length(Depth)) {

IntegratedKd[i] <- IntegratedKd[i-1] + (Depth[i]-Depth[i-1])*KdSpectra[i,iwav]

}

EdProfile[,iwav] <- Ed0[Lam==WL[iwav]]*exp((-1)*IntegratedKd)

}

for(idepth in 1:length(Depth)) {

EdSpectra <- EdProfile[idepth,]

Sum <- 0

for(iwav in 2:301) {

Sum <- (EdSpectra[iwav]*WL[iwav] + EdSpectra[iwav-1]*WL[iwav-1])/2 + Sum

}

PARProfile[idepth] <-  Sum/planck/lightspeed/10^3/Avogadro 

}

KdPARProfile[1] <- (-1)*log(PARProfile[1]/PAR0)/Depth[1]

for(idepth in 2:length(Depth)) {

KdPARProfile[idepth] <-  (-1)*log(PARProfile[idepth]/PARProfile[idepth-1])/(Depth[idepth] - Depth[idepth-1])

}


return(list(KdPARProfile=KdPARProfile, PARProfile=PARProfile, PAR0=PAR0))

}