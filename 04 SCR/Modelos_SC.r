#==============================================================================#
#                                                                              #
#                       SCR - UNMARKED (CONTEOS ESPACIALES)                    #
#                     Seguimiento de la Diversidad Biológica                   #
#                          José Jiménez (CSIC-IREC)                            #
#                      UNIVERSIDAD DE CASTILLA-LA MANCHA                       #
#                             27/02/2024 9:01:34                               #
#                                                                              #
#==============================================================================#

# Con este codigo se trata la informacion con el script del libro de Royle et
# al. (2013)
setwd('C:/...')
source("SCR_functions.R")
library(nimble)
library(MASS)
library(coda)
library(lattice)

tr<-seq(1.5,8.5, length=10)
X<-cbind(rep(tr,each=length(tr)),rep(tr,times=length(tr))) # 100 coord. trampas

plot(X, xlim=c(0,10), ylim=c(0,10), pch=3, cex=0.75, xlab="X", ylab="Y", asp=TRUE)

set.seed(2225)
xlim <- c(0,10); ylim <- c(0,10)
N <- 20

S <- cbind(runif(N, xlim[1], xlim[2]), runif(N, ylim[1], ylim[2]))
points(S, pch=16, col=2)

J <- nrow(X)
K <- 15 
sigma <- 0.5
lambda0 <- 0.4
r<-3
yy <- array(NA, c(N, J, K))
for(j in 1:J) {
  dist <- sqrt((X[j,1]-S[,1])^2 + (X[j,2]-S[,2])^2)
  lambda <- lambda0*exp(-dist^2/(2*sigma^2))
  for(k in 1:K) {
    yy[,j,k] <- rpois(N, lambda)
  }
}
sum(yy)

n <- apply(yy, c(2,3), sum);
n[is.na(n)] <- 0
n1<-apply(n, 1, sum)
sum(n)
M<-150

# Lo anadimos al ploteado:
tot<-apply(n, 1,sum)
symbols(X, circles=tot/50, inches=F,bg="#00000022", fg=NULL, add=T)
points(X, pch=3, cex=0.75); points(S, pch=16, col=2)
points(S, pch=16, col="red", cex=1.25)

# SC sin informacion adicional es muy poco preciso. Vamos a usar un prior
# informativo para sigma, generandolo con una distribucion gamma. Si queremos
# usar un sigma:
mode = 0.5
sd = 0.1

# Obtenemos los parámetros de ratio y forma:
ra = ( mode + sqrt( mode^2 + 4*sd^2 ) ) / ( 2 * sd^2 )
sh = 1 + mode * ra

show(sh)
show(ra)

# Ploteamos:
X11()
x = seq(0,mode+5*sd,len=1001)
plot( x , dgamma( x , shape=sh , rate=ra ) , type="l" ,
      main=paste("dgamma, mode=",mode,", sd=",sd,sep="") ,
      ylab=paste("dgamma( shape=",signif(sh,3)," , rate=",signif(ra,3)," )",
                 sep="") )
abline( v=mode , lty="dotted" )


library(nimble)
## define the model
code <- nimbleCode({

  sigma ~ dgamma(sh,ra)
  psi ~ dunif(0,1)
  lam0 ~ dunif(0,5)
  r ~ dgamma(.1,.1)
  
  for(i in 1:M) {
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    d2[i,1:J] <- (s[i,1]-X[1:J,1])^2 + (s[i,2]-X[1:J,2])^2
    lam[i,1:J] <- lam0*exp(-d2[i,1:J]/(2*sigma^2))*z[i]*K
  }

  for(j in 1:J){
    bigLambda[j] <- sum(lam[1:M,j])
    n[j] ~ dpois(bigLambda[j])
  }
  N <- sum(z[1:M])

})



constants <- list(M = M, K=K, J=J, sh=sh, ra=ra)
data<-list(n=n1, X=X, xlim=xlim, ylim=ylim)

s.start <- cbind(runif(M, 0, 10), runif(M, 0, 10))
# Funciones de inicio
d <- e2dist(s.start[1:M,], X)
lam <- 0.2 * exp( -(d^2)/(2 * 0.5^2))

yi <- array(0, c(M, J, K)) # resighting array
for (j in 1:J) {
  for (k in 1:K) {
    if (n[j, k] > 0) {
      probs <- lam[ ,j]
      probs <- probs / sum(probs)
      latent.id <- sample(1:M, n[j,k], prob = probs, replace = FALSE)
      yi[latent.id , j, k] <- 1
    }
  } # end of k
}   # end of j

yis <- apply(yi, c(1,2), sum)

ytot<- apply(yis,c(1,2),sum)
ssarr<- array(NA,dim=c(M,2))
  for(i in 1:M){
    if(sum(ytot[i,])==0) next
    s.start[i,1]<- mean( X[ytot[i,]>0,1] )
    s.start[i,2]<- mean( X[ytot[i,]>0,2] )
  }
  ssarr[,]<- s.start


z<-rep(1,M)
inits <- list (sigma=5, lam0=runif(1,0,1), lam=yis, s=s.start, z=z)

Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits)
Cmodel <- compileNimble(Rmodel)

mcmcspec<-configureMCMC(Rmodel, monitors=c('N','lam0','psi', 'sigma','r'), thin=1)
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


pumpMCMC <- buildMCMC(mcmcspec)

CpumpMCMC <- compileNimble(pumpMCMC, project = Rmodel)

# Ejecutamos el modelo
nb=1000      # Iteraciones a desechar
ni=5000 +nb  # Iteraciones
nc=3         # Cadenas


start.time2<-Sys.time()
outNim <- runMCMC(CpumpMCMC, niter = ni , nburnin = nb , nchains = nc, inits=inits,
                  setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)
end.time<-Sys.time()
end.time-start.time2 # tiempo de ejecucion

summary(outNim[,c('N','lam0','psi', 'sigma','r')])
xyplot(outNim[,c('N','lam0','psi', 'sigma','r')])

gelman.diag(outNim[,c('N','lam0','psi', 'sigma','r')], multivariate = FALSE)

cat("Poblacion que simulamos = ", N, "individuos", "\n")
cat("Fotografias (todos no identificados)", sum(n), "\n")
