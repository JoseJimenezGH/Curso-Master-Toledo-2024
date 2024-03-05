#==============================================================================#
#                                                                              #
#                    CAPTURA-RECAPTURA ESPACIAL (SCR)                          #
#                       José Jiménez (CSIC-IREC)                               #
#                         19/07/2022 22:38:14                                  #
#                                                                              #
#==============================================================================#
setwd('C:/Users/Jose/OneDrive/00 Master Seguimiento de la Diversidad Biológica 2023/Curso/AulaVirtual/04 SCR/Modelos avanzados')
source("SCR_functions.R")
library(scrbook)
library(lattice)
library(coda)
library(mcmcOutput)
library(nimble)
nimble:::setNimbleOption('useSafeDeparse', FALSE)
nimbleOptions('useSafeDeparse')
source("sSampler.R")


## Simulación de datos
N <- 50  # Tamaño de población
K <- 15  # Ocasiones de muestreo
J <- 100 # número de trampas

data <- SimSCR0(N=N, K=K, array3d = TRUE, discard0=TRUE, rnd=2013)

p0<-data$p0
sigma<-data$sigma
N<-data$N

S<-data$S
y3d<-data$Y
(detections <- sum(y3d))

(nind <- dim(y3d)[1])
X <- data$traplocs
K <- data$K
J <- nrow(X) 
M <- 150

# Espacio de estados
xlim <- data$xlim
ylim <- data$ylim
area <- diff(xlim)*diff(ylim)

y <- apply(y3d,c(1,2),sum); sum(y)
yaug<-array(0,c(M,J))
yaug[1:nind,]<-y

# Ploteamos los capturados vs no capturados (puntos sólidos vs círculos)
plot(X, xlim=xlim, ylim=ylim, pch="+", asp=TRUE)
points(S, pch=1, col="red", cex=1.5)
captured<-apply(y3d,1,sum)
captured[captured>1]<-1
captured<-(1:nind)*captured

points(S[captured,], pch=16, col="red", cex=1.5)
datn<-apply(y3d, c(2,3), sum)
tot<-apply(datn, 1,sum)
symbols(X, circles=tot/5, inches=F, bg="#00000022", fg=NULL, add=T)
points(X, pch="+", cex=0.5)


# Spiderplot
spiderplotJJ4(y3d, X, buffer=2, lwd=2)


## definimos el modelo
code <- nimbleCode({
  
  p0 ~ dbeta(1,1)
  sigma ~ dunif(0,10)
  psi ~ dunif(0,1)
  
  for(i in 1:M){
    z[i] ~ dbern(psi) 
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    p[i,1:J] <- GetDetectionRate(s = s[i,1:2], X = X[1:J,1:2], J=J, sigma=sigma, lam0=p0, z=z[i])
    y[i,1:J] ~ dBinomVector(P = p[i,1:J], KT=KT[1:J])  
  }
  N <- sum(z[1:M])
  D <- N/area
})

# Function to calculate detection rate, but skip when z=0
GetDetectionRate <- nimbleFunction(
  run = function(s = double(1), lam0=double(0), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
     d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
     ans <- lam0*exp(-d2/(2*sigma^2))
     return(ans)
    }
  }
)

# Vectorized binomial
dBinomVector <- nimbleFunction(
  run = function(x = double(1), P = double(1), KT = double(1), log = integer(0, default = 0)) {
    J <- length(x)
    ans <- 0.0
    for(j in 1:J)
      ans <- ans + dbinom(x[j], KT[j], P[j], 1)
    returnType(double())
    if(log) return(ans)
    else return(exp(ans))
  })

rBinomVector  <- nimbleFunction(
  run = function(n = integer(), P = double(1), KT = double(1)) {
    J <- length(P)
    ans<- numeric(J)
    for(j in 1:J)
      ans[j] <- rbinom(1, KT[j], P[j])
    returnType(double(1))
    return(ans)
  })

registerDistributions(list(
  dBinomVector = list(
    BUGSdist = "dBinomVector(P, KT)",
    Rdist = "dBinomVector(P, KT)",
    discrete = TRUE,
    range = c(0, Inf),
    types = c('value = double(1)', 'P = double(1)', 'KT = double(1)'))
))

str(constants <- list(M=M,        # aumentado de datos 
                      KT=rep(K,J),# ocasiones de muestreo
                      J=J,        # número de trampas
                      area=area)) # área del espacio de estados

str( data   <-   list(y=yaug,     # matriz de capturas reducida sobre K
                      X=X,        # matriz de coordenadas de las trampas
                      xlim=xlim,  # extremos x del espacio de estados
                      ylim=ylim)) # extremos y del espacio de estados

# Inicios para las ubicaciones latentes
sst <- cbind(runif(M,xlim[1],xlim[2]),runif(M,ylim[1],ylim[2]))
for(i in 1:nind){
  sst[i,1] <- mean( X[yaug[i,]>0,1] )
  sst[i,2] <- mean( X[yaug[i,]>0,2] )
}
zd<-c(rep(1, nind), rbinom((M-nind),1,0.2))

str( inits   <-  list(p0=runif(1,0,1),    # probabilidad basal de detección
                      sigma=runif(1,0,2), # parametrización de sigma
                      psi=runif(1,0,1),   # parámetro de aumentado de datos
                      s=sst,              # ubicaciones de inicio
                      z=zd))              # está o no en la población

params <- c('N', 'D', 'sigma','psi', 'p0','s','z')

Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits)
Rmodel$initializeInfo()  # Conveniente para ver los inicios
Rmodel$calculate()
Cmodel <- compileNimble(Rmodel)
conf<-configureMCMC(Rmodel, monitors=params)

conf$removeSampler(paste("s[1:",M,", 1:2]", sep="")) 
for(i in 1:M){
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                    type = 'sSampler',
                    control=list(i=i, 
                                 xlim=Rmodel$xlim, 
                                 ylim=Rmodel$ylim, 
                                 scale=0.025), #scale parameter here is just the starting scale. It will be tuned.
                    silent = TRUE)
}

MCMC <- buildMCMC(conf)

CompMCMC <- compileNimble(MCMC, project = Rmodel)

nb = 1000
ni = 5000 + nb
nc = 3

start.time2<-Sys.time()
outNim <- runMCMC(CompMCMC, niter = ni , nburnin = nb , nchains =  nc,inits=inits,
                  setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)
end.time<-Sys.time()
end.time-start.time2 # Time difference of 1.639574 mins

# Resultados
summary(outNim[,c('N','D','p0','psi', 'sigma')])

# Inspeccionamos la convergencia de las cadenas de Markov
#xyplot(outNim[,c('N','p0','psi', 'sigma')])
diagPlot(outNim[,c('N','p0','psi', 'sigma')])

gelman.diag(outNim[,c('N','p0','psi', 'sigma')], multivariate = FALSE)

cat("Población que simulamos = ", N, "individuos", "\n")
cat("p0 simulada  = ", p0, "\n")
cat("sigma simulada  = ", sigma, "\n")
cat("Datos recogidos = ", sum(y), "foto-capturas", "\n")


samplesn<-data.matrix(outNim)

# Coeficiente variación
sd(samplesn[,2])/mean(samplesn[,2])

# Histograma
hist(samplesn[,2], col="grey90", main="", xlab="N")

moda <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

modaN<-moda(samplesn[,2])
abline(v=modaN, lty=2, lwd=3, col="red")

################################################################################

# VISUALIZACIÓN ESPACIAL
#========================
mco <- mcmcOutput(outNim)
s1<-mco$s[,,1] ; Sx<-as.matrix(s1)
s2<-mco$s[,,2] ; Sy<-as.matrix(s2)
z<-mco$z       ; z<- as.matrix(z)

delta<-1.5
Xl<-min(X[,1])-delta
Xu<-max(X[,1])+delta
Yl<-min(X[,2])-delta
Yu<-max(X[,2])+delta
obj<-list(Sx=Sx,Sy=Sy,z=z)

par(mfrow=c(1,1))
dev.new(width=8,height=8)
Spat<-SCRdensity(obj, nx=100, ny=100, Xl=Xl, Xu=Xu, Yl=Yl, Yu=Yu)
points(X, pch="+", col="white")
points(S, col="red", pch=16)
