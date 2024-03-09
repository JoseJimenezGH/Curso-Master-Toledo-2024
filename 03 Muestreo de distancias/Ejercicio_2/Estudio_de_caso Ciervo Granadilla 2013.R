#==============================================================================#
#                                                                              #
#                        MUESTREO DE DISTANCIAS SIMPLE                         #
#                     CASO DE ESTUDIO: CIERVO EN GRANADILLA                    #
#                           Jose Jimenez (CSIC-IREC)                           #
#                             10:44 09/03/2024                                 #
#                                                                              #
#==============================================================================#

library(unmarked)

setwd('C:/Users/Usuario/OneDrive - Universidad de Castilla-La Mancha/Escritorio')

granad13<-read.table("Granad2013.txt", header=TRUE)
x<-sin((granad13[,4]/360)*(2*pi))*granad13[,5]
L<-sum(unique(granad13$Long_TRANS))/100

dists<-granad13

dists$transect<-dists$ID_TRANS
levels(as.factor(unique(dists$ID_TRANS)))
dists$distance<-x

dataDist<-dists[,6:7]
dataDist<-within(dataDist, transect <- factor(transect,
  labels = c('a','b','c','d','e','f','g','h','i','j')))

dist.breaks<-seq(0,150,25)

yDat <- formatDistData(dataDist, distCol="distance",
  transectNameCol="transect", dist.breaks=dist.breaks)

yDat

#fix(yDat)   # Uso esto para quitar los datos de transectos a cero

umf <- unmarkedFrameDS(y=as.matrix(yDat), survey="line",
  dist.breaks=dist.breaks, tlength=10*c(unique(granad13$Long_TRANS)),
  unitsIn="m")
  
summary(umf)
hist(umf)

h_sn_Null <- distsamp(~1~1, umf, keyfun="halfnorm", output="density",
  unitsOut="ha")
h_exp_Null<-distsamp(~1~1, umf, keyfun="exp", output="density",
  unitsOut="ha")
h_hz_Null<-distsamp(~1~1, umf, keyfun="hazard", output="density",
  unitsOut="ha")
h_unif_Null<-distsamp(~1~1, umf, keyfun="uniform", output="density",
  unitsOut="ha")

fms <- fitList(m1=h_sn_Null, m2=h_exp_Null, m3=h_hz_Null, m4=h_unif_Null)
modSel(fms)

backTransform(h_hz_Null, type="state")

# Ploteado de la funcción de  detección
fmhz.shape <- exp(coef(h_hz_Null, type="det"))
fmhz.scale <- exp(coef(h_hz_Null, type="scale"))
plot(function(x) gxhaz(x, shape=fmhz.shape, scale=fmhz.scale), 0, 225,
	xlab="Distancia (m)", ylab="Probabilidad de detección")

# Para superponerlo sobre el histograma
hist(umf, col="lightgrey", border="white", xlab="Distancia", ylab="Frecuencias",
  main="Histograma de distancias")
a<-apply(yDat,2,sum)
lines(a[1]*gxhaz(seq(0:225), shape=fmhz.shape, scale=fmhz.scale), lwd=2)



# Vamos a ver el GOF con *parboot* (boostrapping)
# Esta función nos devuelve tres estadísticos.
fitstats <- function(fm) {
  observed <- getY(fm@data)
  expected <- fitted(fm)
  resids <- residuals(fm)
  sse <- sum(resids^2)
  chisq <- sum((observed - expected)^2 / expected, na.rm=TRUE)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2)
  out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
  return(out)
}

(pb <- parboot(h_hz_Null, fitstats, nsim=500, report=1))
plot(pb)

(c.hat <- pb@t0[2] / mean(pb@t.star[,2]))  # c-hat as ratio of observed 


# USANDO BINOMIAL NEGATIVA
#===========================
umf <- unmarkedFrameGDS(y=as.matrix(yDat), siteCovs=NULL,
  survey="line", dist.breaks=dist.breaks, numPrimary=1,
  tlength=10*c(unique(granad13$Long_TRANS)),
  unitsIn="m")


(ciervo.bn0 <- gdistsamp( ~1, ~1, ~1, data =umf,  
	keyfun="halfnorm", mixture="NB")) 
(ciervo.bn1 <- gdistsamp( ~1, ~1, ~1, data =umf,  
	keyfun="exp", mixture="NB"))
(ciervo.bn2 <- gdistsamp( ~1, ~1, ~1, data =umf, 
	start=c(3,4,0,1),
	keyfun="hazard", mixture="NB")) 


fms2 <- fitList(m1=ciervo.bn0, m2=ciervo.bn1, m3=ciervo.bn2)
modSel(fms2)


backTransform(ciervo.bn0, type="lambda")
backTransform(ciervo.bn0, type="det")


# Vamos a ver el GOF con *parboot* (boostrapping)
# Esta función nos devuelve tres estadísticos.
fitstats <- function(fm) {
  observed <- getY(fm@data)
  expected <- fitted(fm)
  resids <- residuals(fm)
  sse <- sum(resids^2)
  chisq <- sum((observed - expected)^2 / expected, na.rm=TRUE)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2)
  out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
  return(out)
}

(pb <- parboot(ciervo.bn1, fitstats, nsim=500, report=1))
plot(pb)


(c.hat <- pb@t0[2] / mean(pb@t.star[,2]))  # c-hat as ratio of observed 

# Estima empírica bayesiana de la distribución a posteriori:
re.ciervo <- ranef(ciervo.bn2, K=150)
mean(bup(re.ciervo, "mean"))
plot(re.ciervo, layout=c(5,2), xlim=c(0, 50))