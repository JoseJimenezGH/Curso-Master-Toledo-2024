#==============================================================================#
#                                                                              #
#                         MODELOS CORMACK-JOLLY-SEBER                          #
#                        Práctica con el lince ibérico
#                     Seguimiento de la Diversidad Biológica                   #
#                          José Jiménez (CSIC-IREC)                            #
#                      UNIVERSIDAD DE CASTILLA-LA MANCHA                       #
#                            21/02/2024 12:11:16                               #
#                                                                              #
#==============================================================================#



library(jagsUI)


setwd('C:/Users/Usuario/OneDrive/01 IPM Lince/01 data/dataDEF')
CH <- data.matrix(read.table("CJS.txt", header=TRUE))

# Creamos un vector con la ocasión de marcado
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Especificamos el modelo en BUGS
sink("cjs.jags")
cat("
model {
  # Priori y restricciones
  for (i in 1:nind){
    for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- mean.phi
      p[i,t] <- mean.p
    } #t
  } #i
  mean.phi ~ dunif(0, 1)         # Priori para la supervivencia media
  mean.p ~ dunif(0, 1)           # Priori para la recaptura media
  # Definimos la probabilidad
  for (i in 1:nind){
    # Definimos el estado latente en la primera captura
    z[i,f[i]] <- 1
    for (t in (f[i]+1):n.occasions){
      # Proceso de estado
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Proceso de observación
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
    } #t                            
  } #i
}
",fill = TRUE)
sink()

# Función para crear una matriz con información sobre el estado latente 
# conocido z
known.state.cjs <- function(ch){
  state <- ch
  for (i in 1:dim(ch)[1]){
    n1 <- min(which(ch[i,]==1))
    n2 <- max(which(ch[i,]==1))
    state[i,n1:n2] <- 1
    state[i,n1] <- NA
  }
  state[state==0] <- NA
  return(state)
}

# Preparamos datos
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], 
  z = known.state.cjs(CH))

# Función para crear una matriz de valores iniciales para el estado latente z
cjs.init.z <- function(ch,f){
  for (i in 1:dim(ch)[1]){
    if (sum(ch[i,])==1) next
    n2 <- max(which(ch[i,]==1))
    ch[i,f[i]:n2] <- NA
  }
  for (i in 1:dim(ch)[1]){
  ch[i,1:f[i]] <- NA
  }
  return(ch)
}

# Valores iniciales
inits <- function(){list(z = cjs.init.z(CH, f), mean.phi = runif(1, 0, 1),
  mean.p = runif(1, 0, 1))}

# Parametros a monitorizar
parameters <- c("mean.phi", "mean.p", "phi.t")

# Configuración MCMC
ni <- 100000
nt <- 6
nb <- 50000
nc <- 3

library(jagsUI)
# Llamamos a JAGS desde R
cjs <- jags(jags.data, inits, parameters, "cjs.jags", n.chains = nc, 
  n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE)

# Resumen de resultados
print(cjs, digits = 3)

par(mfrow = c(1, 2), las = 1)
hist(cjs$sims.list$mean.phi, nclass = 30, col = "gray", main = "", xlab = "Lynx survival", ylab = "Frequency")
abline(v = cjs$mean$mean.phi, col = "red", lwd = 2)
hist(cjs$sims.list$mean.p, nclass = 30, col = "gray", main = "", xlab = "Recapture", ylab = "")
abline(v = cjs$mean$mean.p, col = "red", lwd = 2)

