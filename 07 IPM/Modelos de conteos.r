#==============================================================================#
#                                                                              #
#                        MODELOS DE ESPACIO DE ESTADOS                         #
#                          José Jiménez (CSIC-IREC)                            #
#                                     UCLM                                     #
#                              15/03/2023 9:32:40                              #
#       Based on the book "Bayesian population analysis using WinBUGS          #
#                       a hierarchical perspective"                            #
#                     by Marc Kéry & Michael Schaub                            #
#                         (2012, Academic Press)                               #
#                                                                              #
#==============================================================================#

library(nimble)
library(mcmcOutput)
library(MCMCvis)

# Ejemplo real: conteos de parejas reproductoras de buitre negro en Cabañeros
# Especificamos el modelo NIMBLE

code <- nimbleCode({
  # Priors
  logN.est[1] ~ dunif(1,1000)   # Prior para el tamaño inicial de la población 
  mean.r ~ dnorm(0, 0.001)      # Media a priori para el ratio de crecimiento
  tau.proc ~ dinvgamma(2,1)
  tau.obs ~ dinvgamma(2,1)
  
  # Verosimilitud
  # Proceso de estados
  for (t in 1:(T-1)){
    r[t] ~ dnorm(mean.r, tau.proc)
    logN.est[t+1] <- logN.est[t] + r[t]
  }
  # Proceso de observación
  for (t in 1:T) {
    y[t] ~ dnorm(logN.est[t], tau.obs)
  }
  
  # Tamaño de población en escala real
  for (t in 1:T) {
    N.est[t] <- exp(logN.est[t])
  }
})


pyears <- 10 # Estima hasta el 2022
hm <- c(63,82,95,83,112,81,82,80,97,102,88,91,118,128,123,111,106,111,
        113,104,140,144,156,162, rep(NA, pyears))
year <- 1989:(2012 + pyears)

# Preparamos los datos
data <- list(y = log(hm))

constants<-list(T = length(year))

# Valores de inicio
inits <- list(tau.proc = runif(1, 0, 1), 
              mean.r = rnorm(1), 
              tau.obs = runif(1, 0, 1), 
              logN.est = c(rnorm(1, 5.6, 0.1), rep(NA, (length(year)-1))))

# Parámetros a monitorizar
params <- c("r", "mean.r", "tau.obs", "tau.proc", "N.est")

Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits,
  check=FALSE, calculate=FALSE)
Rmodel$initializeInfo()
Cmodel <- compileNimble(Rmodel)
conf<-configureMCMC(Rmodel, monitors=params, onlySlice=TRUE, thin=5)

MCMC <- buildMCMC(conf)
CompMCMC<- compileNimble(MCMC, project = Rmodel)

outNim <- runMCMC(CompMCMC, niter = 25000 , nburnin = 1500, nchains = 3, 
                  inits=inits, setSeed = TRUE, progressBar = TRUE, 
                  samplesAsCodaMCMC = TRUE)

options(scipen=999)
MCMCsummary(outNim)
diagPlot(outNim)

hm.ssm<-as.matrix(outNim)
# Producimos gráfico
fitted <- lower <- upper <- numeric()
year <- 1989:2022
n.years <- length(hm)
fitted <- apply(hm.ssm[,1:34],2,mean)
lower <- apply(hm.ssm[,1:34],2,quantile, 0.025)
upper <- apply(hm.ssm[,1:34],2,quantile, 0.975)
m1<-min(lower)
m2<-max(upper)
par(mar = c(4.5, 4, 1, 1))
plot(0, 0, ylim = c(m1, m2), xlim = c(1, n.years), ylab = "Tamaño de población", 
  xlab = "Año", col = "black", type = "l", lwd = 2, axes = FALSE, frame = FALSE)
axis(2, las = 1)
axis(1, at = 1:n.years, labels = year)
polygon(x = c(1:n.years, n.years:1), y = c(lower, upper[n.years:1]), 
  col = "gray90", border = "gray90")
points(hm, type = "l", col = "black", lwd = 2)
points(fitted, type = "l", col = "blue", lwd = 2)
legend(x = 2, y = 550, legend = c("Parejas observadas", "Parejas estimadas"), lty = c(1, 1), 
  lwd = c(2, 2), col = c("black", "blue"), bty = "n", cex = 1)
# Probabilidad de que N(2022) > 200 parejas
mean(hm.ssm[,34] >200)


#... y cual ha sido la realidad?
points(34,281,pch=16, col="red",cex=2)




