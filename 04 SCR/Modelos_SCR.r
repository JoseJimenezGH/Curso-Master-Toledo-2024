#==============================================================================#
#                                                                              #
#               MODELOS DE CAPTURA-RECAPTURA ESPACIALMENTE EXPLÍCITOS          #
#                     Seguimiento de la Diversidad Biológica                   #
#                          José Jiménez (CSIC-IREC)                            #
#                      UNIVERSIDAD DE CASTILLA-LA MANCHA                       #
#                             27/02/2024 9:01:34                               #
#                                                                              #
#==============================================================================#
setwd('C:/Users/Administrator/OneDrive/66 Master Seguimiento de la Diversidad Biológica/Lab/03 SCR/')
source("SCR_functions.R")
library(scrbook)
library(spatstat)
library(lattice)
library(coda)

## Simulación de datos
N <- 50  # Tamaño de población
K <- 15  # Ocasiones de muestreo
J <- 100 # número de trampas

data <- SimSCR0(N=N, K=K, array3d = TRUE, discard0=TRUE, rnd=2013)

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


library(nimble)
## definimos el modelo
code <- nimbleCode({

  alpha0 ~ dnorm(0,0.1)
  logit(p0) <- alpha0
  alpha1 ~ dnorm(0,0.1)
  sigma <- sqrt(1/(2*alpha1))
  psi ~ dunif(0,1)
  
  for(i in 1:M){
    z[i] ~ dbern(psi) 
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    d2[i,1:J] <- pow(s[i,1]-X[1:J,1],2) + pow(s[i,2]-X[1:J,2],2)
    p[i,1:J] <- p0*exp(-alpha1*d2[i,1:J])*z[i]

    for(j in 1:J){
      y[i,j] ~ dbinom(p[i,j],K)
    }
  }

  N <- sum(z[1:M])
  D <- N/area
})

constants <- list(M = M, K=K, J=J, area=area)

dataN<- list (y=yaug, X=X, xlim=xlim, ylim=ylim)


## Inicios para s
sst <- cbind(runif(M,xlim[1],xlim[2]),runif(M,ylim[1],ylim[2]))
for(i in 1:nind){
  sst[i,1] <- mean( X[y[i,]>0,1] )
  sst[i,2] <- mean( X[y[i,]>0,2] )
}

inits <- list (p0=0.5, alpha1=1, s=sst,
               z=c(rep(1, nind), rbinom((M-nind),1,0.2)))


# Preparamos el modelo para ejecución en Nimble
Rmodel <- nimbleModel(code=code, constants=constants, data=dataN, inits=inits, check=FALSE)
Cmodel <- compileNimble(Rmodel)
# Establecemos los parámetros a monitorizar
mcmcspec<-configureMCMC(Rmodel, monitors=c('N', 'D', 'sigma','psi', 'p0','s','z'))

# Cambiamos el muestrador de z (opcional)
mcmcspec$removeSamplers('z')
for(node in Rmodel$expandNodeNames('z')) mcmcspec$addSampler(target = node, type = 'slice')

mcmcspec$removeSamplers("s")
ACnodes <- paste0("s[", 1:constants$M, ", 1:2]")
for(node in ACnodes) {
    mcmcspec$addSampler(target = node,
                        type = "RW_block",
                        control = list(adaptScaleOnly = TRUE),
                        silent = TRUE)
}

# Contruimos el modelo
scrMCMC <- buildMCMC(mcmcspec)

# Compilamos
CSCRMCMC <- compileNimble(scrMCMC, project = Rmodel)

# Ejecutamos el modelo
nb=1000      # Iteraciones a desechar
ni=5000 +nb  # Iteraciones
nc=3         # Cadenas


start.time2<-Sys.time()
outNim <- runMCMC(CSCRMCMC, niter = ni , nburnin = nb , nchains = nc, inits=inits,
                  setSeed = TRUE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)  
end.time<-Sys.time()
end.time-start.time2 # tiempo de ejecución

# Resultados
summary(outNim[,c('N','D','p0','psi', 'sigma')])

# Inspeccionamos la convergencia de las cadenas de Markov
xyplot(outNim[,c('N','p0','psi', 'sigma')])

gelman.diag(outNim[,c('N','p0','psi', 'sigma')], multivariate = FALSE)

cat("Población que simulamos = ", data$N, "individuos", "\n")
cat("lambda0 simulada  = ", data$p0, "\n")
cat("sigma simulada  = ", data$sigma, "\n")
cat("Datos recogidos = ", sum(y), "foto-capturas", "\n")


samplesn<-data.matrix(outNim)

# Coeficiente variación
sd(samplesn[,2])/mean(samplesn[,2])

# Histograma
hist(samplesn[,2], col="grey90")

################################################################################

# VISUALIZACIÓN ESPACIAL
#========================

s1 <- samplesn[,c(5:154)]   ; Sx<-as.matrix(s1)
s2 <- samplesn[,c(155:304)] ; Sy<-as.matrix(s2)
z <-  samplesn[,c(306:455)] ; z<- as.matrix(z)

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
