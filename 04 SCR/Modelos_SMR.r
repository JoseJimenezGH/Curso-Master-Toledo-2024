#==============================================================================#
#                                                                              #
#                        MARCAJE-RECAPTURA ESPACIAL                            #
#                     Seguimiento de la Diversidad Biológica                   #
#                          José Jiménez (CSIC-IREC)                            #
#                      UNIVERSIDAD DE CASTILLA-LA MANCHA                       #
#                             27/02/2024 9:01:34                               #
#                                                                              #
#==============================================================================#

setwd('C:/...')
library(scrbook)
library(coda)
library(lattice)
library(mcmcOutput)
source("SCR_functions.R")
set.seed(1960)
N<-50
m<-15
sigma<-1.2
lam0<-0.2
K<-5

gx<-gy<-seq(0,6,1)
X<-as.matrix(expand.grid(gx,gy))
J<-dim(X)[1]

# Límites de S
xlims<-ylims<-c(-1.5, 7.5)

# Simulamos datos
dat<-sim.pID.Data(N=N, K=K, sigma=sigma, lam0=lam0, knownID=m, X=X,
  xlims=xlims, ylims=ylims, obsmod="pois", nmarked="known")

plot(X, pch=3, xlim=c(min(X[,1])-2,max(X[,1])+2), ylim=c(min(X[,2])-2,max(X[,2])+2))
points(dat$S, col="red", pch=1)
S<-dat$S
Y <- dat$Y; yt<-apply(Y, c(1,2), sum)
detections <- rowSums(yt)>0
points(S[detections,], pch=16, col="red")

# Marcados
yobs<- dat$Yobs
y<-apply(yobs, c(1,2), sum)
nMarked<-dim(y)[1]

# Preparamos datos
n<-dat$n-apply(dat$Yknown, 2:3, sum)
n1<-apply(n,1,sum)

# Desechamos las fotos que no sabemos si corresponden a animales identificados
# o no. SMR no usa los registros de animales de estatus de identificación
# desconocido

X<-matrix(X, ncol=2)

# Ploteado de capturas
plot(X, pch=3, col="blue", cex=0.5, xlim=c(min(X[,1])-2,max(X[,1])+2),
  ylim=c(min(X[,2])-2,max(X[,2])+2), main="")
tot<-apply(dat$n, 1,sum)
symbols(X, circles=tot/15, inches=F,bg="#00000022", fg=NULL, add=T)
points(X, pch=3, col="blue", cex=0.5)
points(S, col="red", pch=1)
points(S[detections,], pch=16, col="red")

M<-100        # Aumentado de datos


library(nimble)
## Definimos modelo
code <- nimbleCode({

  lam0 ~ dunif(0,5)
  sig ~ dunif(0,5)
  sig2 <- 2*sig^2
  psi ~ dbeta(1,1)

  for(i in 1:M) {
    s[i,1] ~ dunif(xlim[1], xlim[2])
    s[i,2] ~ dunif(ylim[1], ylim[2])
    z[i] ~ dbern(psi)
    d2[i,1:J] <- (s[i,1]-x[1:J,1])^2 + (s[i,2]-x[1:J,2])^2
    lam[i,1:J] <- lam0*exp(-d2[i,1:J]/sig2)*z[i]*K
  }
  
  for(i in 1:nMarked) {
    for(j in 1:J) {
      y[i,j] ~ dpois(lam[i,j])
    }
  }

  for(j in 1:J) {
    Lam[j] <- sum(lam[((nMarked+1):M),j])
    n[j] ~ dpois(Lam[j])
  }

  N <- sum(z[1:M])
  D <- N/A

})


## Preparamos constantes
constants<-list(J=J, M=M, K=K,
            nMarked=nMarked,
            xlim=xlims, ylim=ylims,
            A=diff(xlims)*diff(ylims))

## Preparamos datos
data <- list(y=y, n=n1, x=X)

## Preparamos inicios
sst<-cbind(runif(M, xlims[1], xlims[2]), runif(M, ylims[1], ylims[2]))
inits <- list(z=rep(1,M), s=sst, lam0=0.5, sig=0.5)
params <- c("psi", "lam0", "sig", "N", "D")


Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits, check=FALSE)
Cmodel <- compileNimble(Rmodel)
mcmcspec<-configureMCMC(Rmodel, monitors=params)

mcmcspec$removeSamplers('z')
for(node in Rmodel$expandNodeNames('z')) mcmcspec$addSampler(target = node, type = 'slice')

pumpMCMC <- buildMCMC(mcmcspec)

CpumpMCMC <- compileNimble(pumpMCMC, project = Rmodel)

## Output:
library(coda)
library(lattice)

# Resultados usando 1 cadena
CpumpMCMC$run(10000)
samples <- as.matrix(CpumpMCMC$mvSamples)
res<-mcmc.list(mcmc(samples))
summary(window(res, start=500), dig=3)

cat("Población que simulamos = ", N, "individuos", "\n")
cat("Fotografías de animales identificados)= ", sum(y), "\n")
cat("Fotografías de animales no identificados", sum(n), "\n")


xyplot(window(res, start=500))
dat<-data.matrix(window(res, start=1000))

# Coeficiente variaci—n
sd(dat[,1])/mean(dat[,1])

# Histograma
hist(dat[,1], breaks=20, col="grey90")


