#==============================================================================#
#                                                                              #
#                 Fitting red fox density estimates (La Nava)                  #
#                using distance sampling. R + unmarked script.                 #
#                                                                              #
#   GENERALIZED SPATIAL MARK-RESIGHT MODELS WITH INCOMPLETE IDENTIFICATION:    #
#               AN APPLICATION TO RED FOX DENSITY ESTIMATES                    #
#                                                                              #
#                           Ecology and Evolution                              #
#                                                                              #
#             José Jiménez(1), Richard Chandler(2), Jordi Tobajas(1),          #
#             Esther Descalzo(1), Rafael Mateo(1), Pablo Ferreras(1)           #
#                                                                              #
#                             13/02/2019 8:10:15                               #
#                                                                              #
#==============================================================================#

# (1) Instituto de Investigación en Recursos Cinegéticos (IREC, CSIC-UCLM-JCCM),
#   Ronda de Toledo 12, 13071 Ciudad Real, Spain.
# (2) University of Georgia, Warnell School of Forestry and Natural Resources.

setwd('C:/Users/Jose/OneDrive/00 Master Seguimiento de la Diversidad Biológica 2023/Curso/AulaVirtual/03 Muestreo de distancias/Ejercicio')
library(unmarked)

dat2<-read.table("ZorroData.txt", header=TRUE)
dat2$visit<-as.factor(dat2$visit)
dat2$site<-as.factor(dat2$site)
transect.length <- read.table("LTran.txt", header=FALSE)

range(dat2$d)
cutpt <- seq(0, 400, by=50)
y2 <- formatDistData(dat2, "distance", "site", cutpt, "visit")

y2<-y2[order(as.numeric(rownames(y2))),,drop=FALSE]

##
##
##
##  STACKED analysis of Red Fox
##
##

ystacked<- rbind(y2[,1:8], y2[,9:16], y2[,17:24], y2[,25:32], y2[,33:40])

# Package into unmarked GDS data frame and inspect the data
umf.stacked <- unmarkedFrameDS(y = ystacked, survey="line", unitsIn="m",
   dist.breaks=cutpt, tlength=c(transect.length[,1]),
   siteCovs = NULL)


# USING POISSON OBS.MODEL WITH distsamp
#=======================================

# half-normal detection function
fm.hn.p <- distsamp(~1 ~1, 
   keyfun = "halfnorm", output = "density", unitsOut = "kmsq", data = umf.stacked)
# half-normal detection function
fm.exp <- distsamp(~1 ~1,
   keyfun = "exp", output = "density", unitsOut = "kmsq", data = umf.stacked)
# half-normal detection function
fm.hz <- distsamp(~1 ~1,  start=c(0.2, 5.6, 2),
   keyfun = "hazard", output = "density", unitsOut = "kmsq", data = umf.stacked)

f1 <- fitList('Halfnormal'    =fm.hn.p,
              'Exponential'   =fm.exp,
              'Hazard rate'   =fm.hz)

modSel(f1)

##             nPars    AIC delta AICwt cumltvWt
## Hazard rate     3 166.24  0.00 0.754     0.75
## Halfnormal      2 168.67  2.43 0.224     0.98
## Exponential     2 173.31  7.07 0.022     1.00

backTransform(fm.hz, type="state")

hist(fm.hz, lwd=3, col="lightgrey", ylim=c(0,.0075))

# Goodness of fitting model (GoF)
fitstats <- function(fm) {
  observed <- getY(fm@data)
  expected <- fitted(fm, na.rm = TRUE)
  resids <- residuals(fm, na.rm = TRUE)
  sse <- sum(resids^2, na.rm = TRUE)
  chisq <- sum((observed - expected)^2 / expected, na.rm = TRUE)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2, na.rm = TRUE)
  out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
  return(out)
}

(pb <- parboot(fm.hn.p, fitstats, nsim=1000, report=1))
plot(pb, col="grey90")
# chisquare test statistic
(c.hat <- pb@t0[2] / mean(pb@t.star[,2])) # c-hat as ratio of observed/expected


# USING NB OBS. MODEL WITH gdistsamp
#====================================

# Package into unmarked GDS data frame and inspect the data
umf.stacked.v2 <- unmarkedFrameGDS(y = ystacked, survey="line", unitsIn="m",
   dist.breaks=cutpt, numPrimary=1, tlength=c(transect.length[,1]),
   siteCovs =NULL)  

# half-normal-Poisson detection function
fm0.hn.p <- gdistsamp(~1, ~1, ~1,
   keyfun = "halfnorm", output = "density", unitsOut = "kmsq", mixture="P", data = umf.stacked.v2)
fm0.hn.p
# half-normal-Negative Binomial detection function
fm0.hn.NB <- gdistsamp(~1, ~1, ~1,
   keyfun = "halfnorm", output = "density", unitsOut = "kmsq", mixture="NB", data = umf.stacked.v2)
fm0.hn.NB 
# Exponential-Poisson detection function
fm0.exp.p <- gdistsamp(~1, ~1, ~1,
   keyfun = "exp", output = "density", unitsOut = "kmsq", mixture="P", data = umf.stacked.v2)
fm0.exp.p
# Exponential-Negative Binomial detection function
fm0.exp.NB <- gdistsamp(~1, ~1, ~1,
   keyfun = "exp", output = "density", unitsOut = "kmsq", mixture="NB", data = umf.stacked.v2)
fm0.exp.NB
# Hazard-rate-Poisson detection function
fm0.hr.p <- gdistsamp(~1, ~1, ~1,  start=c(0.3, 5, 1),
   keyfun = "hazard", output = "density", unitsOut = "kmsq", mixture="P", data = umf.stacked.v2)
fm0.hr.p
# Hazard-rate-Negative Binomial detection function
fm0.hr.NB <- gdistsamp(~1, ~1, ~1, start=c(-0.43, 5, 1, 0.9),
   keyfun = "hazard", output = "density", unitsOut = "kmsq", mixture="NB", data = umf.stacked.v2)
fm0.hr.NB

fl<-fitList(fm0.hn.p, fm0.hn.NB, fm0.exp.p, fm0.exp.NB, fm0.hr.p, fm0.hr.NB)
modSel(fl)

##            nPars    AIC delta AICwt cumltvWt
## fm0.hr.NB      4 152.17  0.00 0.433     0.43
## fm0.hr.p       3 152.96  0.78 0.293     0.73
## fm0.hn.NB      3 154.60  2.43 0.129     0.85
## fm0.hn.p       2 155.38  3.21 0.087     0.94
## fm0.exp.NB     3 157.20  5.03 0.035     0.98
## fm0.exp.p      2 157.98  5.81 0.024     1.00


backTransform(fm0.hr.NB, type="lambda")

# Goodness of fitting model (GoF)
fitstats <- function(fm) {
  observed <- getY(fm@data)
  expected <- fitted(fm, na.rm = TRUE)
  resids <- residuals(fm, na.rm = TRUE)
  sse <- sum(resids^2, na.rm = TRUE)
  chisq <- sum((observed - expected)^2 / expected, na.rm = TRUE)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2, na.rm = TRUE)
  out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
  return(out)
}

fmP <- try(gdistsamp(~1, ~1, ~1, umf.stacked.v2, output="density", 
  keyfun ="hazard", mixture="NB", K=1000, method="Nelder-Mead", 
  starts=c(-4.3, 5, 1, 0.9),
  control=list(trace=FALSE, REPORT=1, maxit = 50000)))
  if (class(fmP) == "try-error") {v<-1} else {
  tmpP <- summary(fmP)
  (pbP <- parboot(fmP, fitstats, nsim=1000, report=1))
}

par(mfrow=c(2,2))
plot(pbP, col="grey90")

# chisquare test statistic
(c.hat <- pbP@t0[2] / mean(pbP@t.star[,2])) # c-hat as ratio of observed/expected
