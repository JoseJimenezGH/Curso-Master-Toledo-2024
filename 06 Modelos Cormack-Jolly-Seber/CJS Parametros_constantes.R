#==============================================================================#
#                                                                              #
#                         MODELOS CORMACK-JOLLY-SEBER                          #
#                     Seguimiento de la Diversidad Biológica                   #
#                          José Jiménez (CSIC-IREC)                            #
#                      UNIVERSIDAD DE CASTILLA-LA MANCHA                       #
#                            21/02/2024 12:11:16                               #
#                                                                              #
#==============================================================================#

# Modelo con parámetros constantes
# Definimos los valores de parámetro
set.seed(1960)
n.occasions <- 10                  # Número de ocasiones captura
marked <- rep(50, n.occasions-1)  # Número anual de nuevos individuos marcados
phi <- rep(0.65, n.occasions-1)
p <- rep(0.5, n.occasions-1)

# Definimos matrices con probabilidades de supervivencia y recaptura
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked))
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))

# Definimos funciones para simular la matriz de historias de captura (CH)
simul.cjs <- function(PHI, P, marked){
  n.occasions <- dim(PHI)[2] + 1
  CH <- matrix(0, ncol = n.occasions, nrow = sum(marked))
  # Definimos un vector con la ocasión de marcado
  mark.occ <- rep(1:length(marked), marked[1:length(marked)])
  # Llenamos la matriz CH
  for (i in 1:sum(marked)){
    CH[i, mark.occ[i]] <- 1       # Escribimos un 1 en la ocasión de suelta
    if (mark.occ[i]==n.occasions) next
    for (t in (mark.occ[i]+1):n.occasions){
      # Bernoulli: ¿El individuo sobrevive en la ocasión?
      sur <- rbinom(1, 1, PHI[i,t-1])
      if (sur==0) break		# Si muere, nos vamos al siguiente individuo
      # Bernoulli: ¿Se recaptura el individuo?
      rp <- rbinom(1, 1, P[i,t-1])
      if (rp==1) CH[i,t] <- 1
    } #t
  } #i
  return(CH)
}

# Ejecutamos la función
CH <- simul.cjs(PHI, P, marked)

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
parameters <- c("mean.phi", "mean.p")

# Configuración MCMC
ni <- 10000
nt <- 6
nb <- 5000
nc <- 3

library(jagsUI)
# Llamamos a JAGS desde R
cjs <- jags(jags.data, inits, parameters, "cjs.jags", n.chains = nc, 
  n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE)

# Resumen de resultados
print(cjs, digits = 3)

par(mfrow = c(1, 2), las = 1)
hist(cjs$sims.list$mean.phi, nclass = 30, col = "gray", main = "", xlab = "Survival", ylab = "Frequency")
abline(v = phi, col = "red", lwd = 2)
hist(cjs$sims.list$mean.p, nclass = 30, col = "gray", main = "", xlab = "Recapture", ylab = "")
abline(v = p[1], col = "red", lwd = 2)

