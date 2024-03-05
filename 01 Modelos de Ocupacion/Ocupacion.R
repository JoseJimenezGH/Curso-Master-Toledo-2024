#==============================================================================#
#                                                                              #
#                            MODELOS DE OCUPACIÓN                              #
#                     Seguimiento de la Diversidad Biológica                   #
#                          José Jiménez (CSIC-IREC)                            #
#                      UNIVERSIDAD DE CASTILLA-LA MANCHA                       #
#                            21/02/2024 12:11:16                               #
#                                                                              #
#==============================================================================#

setwd("C:/Users/Administrator/OneDrive/66 Master Seguimiento de la Diversidad Biológica/Lab/01 Ocupacion")

# 1. MODELO SIN COVARIABLES
################################################################################
# Simulación de datos
set.seed(24)
M <- 100                             # Número de sitios
J <- 4                               # repeticiones temporales
y <- matrix(NA, nrow = M, ncol = J)  # Matriz para contener los datos 
                                     # observados
# Valores de los parámetros
psi <- 0.6         # Probabilidad de ocupación o presencia
p <- 0.4           # Probabilidad de detección


# Proceso ecológico
#~~~~~~~~~~~~~~~~~~~~
# Lo que hay. Generamos datos de presencia/ausencia de la especie objetivo.
z <- rbinom(n = M, size = 1, prob = psi)  # R no tiene Bernoulli


# PROCESO DE OBSERVACIÓN
#~~~~~~~~~~~~~~~~~~~~~~~~
# Lo que vemos. Generamos datos de detección/no detección.
for(j in 1:J){
   y[,j] <- rbinom(n = M, size = 1, prob = z*p)
}

# Veamos los datos
sum(z)                  # Sitios realmente ocupados
sum(apply(y, 1, max))   # Sitios que se observa la ocupación
head(cbind(z=z, y))     # Realidad y observación para los sitios 1:6



# MODELO EN UNA APROXIMACIÓN BAYESIANA (USANDO JAGS)
#=====================================================
str( win.data <- list(y = y, M = nrow(y), J = ncol(y)) )

cat(file = "model.txt",
"
model {
   # Información a priori
   psi ~ dunif(0, 1)
   p ~ dunif(0, 1)
   # Probabilidad
   for (i in 1:M) {              # Bucle sobre los sitios
      z[i] ~ dbern(psi)          # Proceso ecológico
      for (j in 1:J) {           # Bucle sobre los muestreos replicados
         y[i,j] ~ dbern(z[i]*p)  # Proceso de observación
      }
   }
}
")



zst <- apply(y, 1, max)   # Para evitar conflictos 
                          # entre datos/modelo/inicio
inits <- function(){list(z = zst)}

params <- c("psi", "p")

ni <- 5000   ;   nt <- 1   ;   nb <- 1000   ;   nc <- 3
library(jagsUI)
fm2 <- jags(win.data, inits, params, "model.txt", n.chains = nc,
   n.thin = nt, n.iter = ni, n.burnin = nb)
print(fm2, dig = 3)


# MODELO EN UNA APROXIMACIÓN FRECUENTISTA (USANDO UNMARKED)
#============================================================
library(unmarked)                # Cargamos la librería
umf <- unmarkedFrameOccu(y = y)  # preparamos datos
summary(umf)                     # Resumen de datos

(fm <- occu(~1 ~1, umf))         # Ejecutamos el modelo
# En escala real
backTransform(fm, type="state")
backTransform(fm, type="det")



# 2. MODELO CON COVARIABLES
################################################################################
set.seed(1)
M <- 100             # Número de sitios
J <- 3               # Número de presencias/ausencias
y <- matrix(NA, nrow = M, ncol = J) # para contener los datos 
                                    # de observaciones

# Creamos una covariable para la ocupación de 'altura de la vegetación' y la 
# llamamos vegHt:
vegHt <- sort(runif(M, -1, 1)) # lo ordenamos porque nos conviene 
                               # para los gráficos

# Elegimos valores de parámetros para el modelo                               
beta0 <- 0           # Intercepto en escala logit
beta1 <- 3           # Pendiente en escala logit para vegHt
psi <- plogis(beta0 + beta1 * vegHt) # Probabilidad de ocupación


# Simulamos para cada sitio para obtener las verdaderas presencias/ausencias
z <- rbinom(M, 1, psi)        # Verdadera presencia/ausencia

# Vemos los datos reales que hemos creado
table(z)
# Ploteamos el verdadero estado (relación entre la covariable y la ocupación)
plot(vegHt, z, xlab="Altura de la vegetación", ylab="ocupación")
plot(function(x) plogis(beta0 + beta1*x), -1, 1, add=T, lwd=3, col = "red")

# Vamos ahora a simular la detección, creando una covariable para ella. Será 
# el 'viento', para el que suponemos que va a disminuir la detección con la 7
# velocidad del viento.
viento <- array(runif(M * J, -1, 1), dim = c(M, J))
alpha0 <- -2                          # intercepto en escala logit
alpha1 <- -3                          # Pendiente en escala logit para viento
p <- plogis(alpha0 + alpha1 * viento) # Probabilidad de detección

# Simulamos ahora la detección con esta covariable. Realizamos J = 3 réplicas 
# del muestreo en cada sitio
for(j in 1:J) {
    y[,j] <- rbinom(M, z, p[,j])
}
sum(apply(y, 1, max))      # Número de sitios con presencias observadas

# Ploteamos los datos observados y efectos del viento en la probabilidad de 
# detección
plot(viento, y, xlab="Viento", ylab="Probabilidad de detección", 
  frame = F, cex = 0.75)
plot(function(x) plogis(alpha0 + alpha1*x), -1, 1, add=T, lwd=3, col = "red")

# Inspeccionamos los datos: ocupación real y lo que percibimos (y)
head(cbind(psi=round(psi,2), z=z, y1=y[,1], y2=y[,2], y3=y[,3]))


# MODELO EN UNA APROXIMACIÓN BAYESIANA USANDO NIMBLE
#====================================================
str(data    <- list(y = y, vegHt = vegHt, viento = viento))
str(constants<-list(M = nrow(y), J = ncol(y),
                     XvegHt =  seq(-1, 1, length.out=100),
                     Xviento = seq(-1, 1, length.out=100)))



library(nimble)
code <- nimbleCode({

  # A priori
  mean.p ~ dunif(0, 1)     # prior del intercepto de deteccion
  alpha0 <- logit(mean.p)  # Intercepto de deteccion
  alpha1 ~ dunif(-20, 20)  # Prior de la covariable deteccion **viento**
  mean.psi ~ dunif(0, 1)   # distribución de la ocupacion
  beta0 <- logit(mean.psi) # Prior del intercepto de ocupación
  beta1 ~ dunif(-20, 20)   # Prior para la covariable ocupación **vegHt**

  # Probabilidad
  for (i in 1:M) {
    # Modelo de estado o proceso ecológico
    z[i] ~ dbern(psi[i])      # Verdadera ocupacion (z) en el sitio i
    logit(psi[i]) <- beta0 + beta1 * vegHt[i]
    for (j in 1:J) {
      # Modelo de observacion para la observacion real
      y[i,j] ~ dbern(p.eff[i,j])    # deteccion-no deteccion en i y j
      p.eff[i,j] <- z[i] * p[i,j]
      logit(p[i,j]) <- alpha0 + alpha1 * viento[i,j]
    }
  }

  # Cantidades derivadas
  N.occ <- sum(z[1:M])       # Numero de sitios ocupados
  psi.fs <- N.occ/M          # Proporcion de sitios ocupados
  for(k in 1:100){
    logit(psi.pred[k]) <- beta0 + beta1 * XvegHt[k] # predicciones de psi
    logit(p.pred[k]) <- alpha0 + alpha1 * Xviento[k]  # predicciones de p
  }
})



zst <- apply(y, 1, max)
str(inits <- list(z = zst, mean.p = runif(1), alpha1 = runif(1),
              mean.psi = runif(1), beta1 = runif(1)))

params <- c('alpha0', 'alpha1', 'beta0', 'beta1', 'N.occ', 'psi.fs',
  'psi.pred', 'p.pred', 'z') 


Rmodel <- nimbleModel(code=code, constants=constants, data=data, 
                      inits=inits, check=FALSE, calculate=FALSE)
Cmodel <- compileNimble(Rmodel)
mcmcspec<-configureMCMC(Rmodel, monitors=params, nthin=10)
pumpMCMC <- buildMCMC(mcmcspec)
CpumpMCMC <- compileNimble(pumpMCMC, project = Rmodel)
## Ejecutamos con estas características
nb=10000
ni=25000+nb
nc=3

start.time<-Sys.time()
outNim <- runMCMC(CpumpMCMC, niter = ni , nburnin = nb , 
                  nchains = nc, inits=inits,
                  setSeed = TRUE, progressBar = TRUE, 
                  samplesAsCodaMCMC = TRUE)
end.time<-Sys.time()
end.time-start.time


# Comparamos la realidad con la inferencia bayesiana en la tabla:
realidad <- c('alpha0'=alpha0, 'alpha1'=alpha1, 'beta0'=beta0, 
              'beta1'=beta1, 'N.occ'=sum(z), 'psi.fs'=sum(z)/M)
samples<-rbind(as.matrix(outNim$chain1),
               as.matrix(outNim$chain2),
               as.matrix(outNim$chain2))
ref0<-apply(samples[,2:5],2,mean)
ref1<-sum(apply(samples[,207:306],2,mean))
print(cbind('realidad'=realidad, 'estimado'=c(ref0, ref1, ref1/M)))

# Ploteado
par(mfrow = c(2, 2), mar = c(5, 5, 2, 2), las = 1, 
  cex.lab = 0.75, cex = 0.8)
plot(function(x) plogis(alpha0 + alpha1*x), -1, 1, 
  lwd=3, col = "red")
plot(vegHt, z, xlab="Altura de la vegetacion", 
  ylab="Probabilidad de ocupacion (P. media)",
  las = 1, frame = F)   # Verdadera presencia/ausencia
lines(vegHt, psi, lwd=3, col="red")   # psi real
lines(constants$XvegHt, apply(samples[,107:206],2,mean), 
  col="blue", lwd = 2)
lower.psi<-apply(samples[,107:206],2,quantile, 0.025)
upper.psi<-apply(samples[,107:206],2,quantile, 0.975)
matlines(constants$XvegHt, cbind(lower.psi,upper.psi), 
         col="grey", lty = 1)

plot(viento, y, xlab="Viento", 
     ylab="Probabilidad de deteccion", frame = F)
plot(function(x) plogis(alpha0 + alpha1*x), -1, 1, 
     add=T, lwd=3, col = "red")
lines(constants$Xviento, apply(samples[,6:105],2,mean), 
      col="blue", lwd = 2)
lower.p<-apply(samples[,6:105],2,quantile, 0.025)
upper.p<-apply(samples[,6:105],2,quantile, 0.975)
matlines(constants$Xviento, cbind(lower.p,upper.p), 
         col="grey", lty = 1)

# Ploteado de posteriores del numero de sitios ocupados
hist(samples[,1], col = "grey", breaks = 60, xlim = c(20, 100),
main = "", freq = F)
abline(v = mean(samples[,1]), col = "blue", lwd = 3) 
lower.Nocc<-quantile(samples[,1],0.025)
upper.Nocc<-quantile(samples[,1],0.975)
abline(v = c(lower.Nocc,upper.Nocc), col = "grey", lwd = 3) 
abline(v = sum(apply(y, 1, max)), lty = 2, lwd = 3)
abline(v = sum(z), col = "red", lwd = 2)

