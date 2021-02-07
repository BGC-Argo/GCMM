
ModelKdMM01 <- function(KdInput, wl=490, range=seq(400,700,5)) {

Table <- read.table(file="D:/R-Program/MM01.txt", header = TRUE, sep = "\t", as.is = TRUE) # Change the File Folder

lambda <- Table$Lambda

Kw <- Table$Kw

e <- Table$e

chi<- Table$chi

KdRes <- range*NA

if(!is.na(KdInput)) {

Kbio <- (KdInput - Kw[lambda==wl])

if (Kbio<0) Kbio <- 0

Chl <- (Kbio/chi[lambda==wl])^(1/e[lambda==wl])

Kd <- Kw + chi*Chl^e

KdRes <- approx(lambda, Kd, range)$y

}

return(KdRes)

}