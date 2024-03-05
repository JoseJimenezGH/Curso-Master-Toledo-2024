#==============================================================================#
#                                                                              #
#                      MODELOS DE POBLACIÓN INTEGRADOS (IPM)                   #
#                      Seguimiento de la Diversidad Biológica                  #
#                           José Jiménez (CSIC-IREC)                           #
#                       UNIVERSIDAD DE CASTILLA-LA MANCHA                      #
#                               27/02/2024 9:01:34                             #
#                                                                              #
#==============================================================================#

library(nimble)
library(popbio)
library(lattice)
library(coda)

# Usamos un modelo pre-reproducción, y asumimos que los individuos se empiezan a 
# reproducir al llegar a la edad adulta, y vamos a distinguir a efectos del 
# modelo dos clases de edad: de 1 año (juv), y de más de 1 año (ad).


# Vamos a ver ahora los datos que usaremos en el modelo. El primer conjunto 
# de datos es el de conteo de machos localizados con reclamo (asumimos, por el 
# conocimiento de la biología de la especie, que hay idéntico número de 
# hembras reproductoras):
# Conteos de machos adultos (desde el año 1 al 10)
y <- c(42, 36, 35, 33, 36, 37, 39, 40, 43, 39)


# El segundo conjunto de datos son los de productividad. Se se han tomado en 
# el seguimiento de nidos con puesta. Son el número de pollos que vuelan (J) 
# y el número de nidos (R):
# Datos de productividad (desde el año 1 al 9)
J <- c(4, 12, 20, 12, 16, 18,  4, 26, 26, 10)   # Pollos volados
R <- c(37, 38, 36, 32, 38, 35, 36, 48, 41, 35)  # Nidos localizados


# El tercer conjunto de datos son de captura-recaptura:
# Datos de captura-recaptura (en formato m-array, desde el año 1 a 10)
m<-matrix(c(1,  0,  0,  0,  0,  0,  0,  0,  0,  2,
            0,  3,  0,  0,  0,  0,  0,  0,  0,  2,
            0,  0,  5,  0,  0,  0,  0,  0,  0,  5,
            0,  0,  0,  1,  0,  0,  0,  0,  0,  2,
            0,  0,  0,  0,  2,  0,  0,  0,  0,  3,
            0,  0,  0,  0,  0,  1,  0,  0,  0,  5,
            0,  0,  0,  0,  0,  0,  2,  1,  0,  3,
            0,  0,  0,  0,  0,  0,  0,  1,  0,  5,
            0,  0,  0,  0,  0,  0,  0,  0,  2,  6,
           17,  7,  1,  0,  0,  0,  0,  0,  0,  2,
            0, 18,  1,  0,  0,  0,  0,  0,  0,  4,
            0,  0, 25,  2,  1,  0,  0,  0,  0,  1,
            0,  0,  0, 24,  4,  1,  1,  0,  0,  2,
            0,  0,  0,  0, 20,  1,  1,  0,  0,  5,
            0,  0,  0,  0,  0, 19,  4,  1,  1,  2,
            0,  0,  0,  0,  0,  0, 17,  1,  0,  4,
            0,  0,  0,  0,  0,  0,  0, 18,  6,  1,
            0,  0,  0,  0,  0,  0,  0,  0, 19,  3), ncol = 10, byrow = TRUE)


# Modelo en BUGS utilizando Nimble
code <- nimbleCode({

#-------------------------------------------------
# MODELO DE POBLACIÓN INTEGRADO 
#  - Modelo estructurado por edades con 2 clases:
#    joven (hasta 1 año) y adulto (al menos 2 años)  
#  - Edad de primera reproducción: 1 año
#  - Conteos pre-reproductores (machos territoriales)
#  - Todos los ratios vitales se asumen constantes
#-------------------------------------------------

  #----------------------------------------------------------
  # 1. Definimos la información a priori para los parámetros
  #----------------------------------------------------------
  # Error de observación
  tauy <- pow(sigma.y, -2)
  sigma.y ~ dunif(0, 50)
  sigma2.y <- pow(sigma.y, 2)

  # Tamaños iniciales de población
  n1 ~ T(dnorm(3, tauy),0,1000)      # 1-año
  nad ~ T(dnorm(35, tauy),0,1000)    # Adultos
  N1[1] <- round(n1)
  Nad[1] <- round(nad)

  # Supervivencia, productividad y probabilidad de recaptura
  for (t in 1:(nyears-1)){
    sjuv[t] <- mean.sjuv
    sad[t] <- mean.sad
    p[t] <- mean.p
    f[t] <- mean.fec
  }

  mean.sjuv ~ dunif(0, 1)
  mean.sad ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  mean.fec ~ dunif(0, 20)

  #-------------------------------------------------
  # 2. Parámetros derivados
  #-------------------------------------------------
  # Ratio de crecimiento poblacional
  for (t in 1:(nyears-1)){
    lambda[t] <- Ntot[t+1] / Ntot[t]
  }

  #-------------------------------------------------
  # 3. Probabilidad de los conjuntos de datos
  #-------------------------------------------------
  # 3.1. Probabilidad de los datos de conteos
     # 3.1.1 Proceso de sistema (realidad biológica)
    for (t in 2:nyears){
      mean1[t] <- f[t-1] / 2 * sjuv[t-1] * Ntot[t-1]
      N1[t] ~ dpois(mean1[t])
      Nad[t] ~ dbin(sad[t-1], Ntot[t-1])
    }
    for (t in 1:nyears){
      Ntot[t] <- Nad[t] + N1[t]
    }

    # 3.1.2 Proceso de observación
    for (t in 1:nyears){
      y[t] ~ dnorm(Ntot[t], tauy)
    }

  # 3.2 Probabilidad de los datos de captura-recaptura:    
  # modelo CJS con dos clases de edad
  
  # Probabilidad multinomial
  for (t in 1:2*(nyears-1)){
    m[t,1:nyears] ~ dmulti(pr[t,1:nyears], r[t])
  }

  # Probabilidades del m-array para juveniles
  for (t in 1:(nyears-1)){
    # Diagonal principal
    q[t] <- 1-p[t]
    pr[t,t] <- sjuv[t] * p[t]
    # Por encima de la diagonal principal
    for (j in (t+1):(nyears-1)){
      pr[t,j] <- sjuv[t]*prod(sad[(t+1):j])*prod(q[t:(j-1)])*p[j]
    } #j
    # Bajo la diagonal principal
    for (j in 1:(t-1)){
      pr[t,j] <- 0
    } #j
    # Última columna: probabilidad de no-recaptura
    pr[t,nyears] <- 1-sum(pr[t,1:(nyears-1)])
  } #t

  # Probabilidades del m-array para adultos
  for (t in 1:(nyears-1)){
    # Diagonal principal
    pr[t+nyears-1,t] <- sad[t] * p[t]
    # Por encima de la diagonal principal
    for (j in (t+1):(nyears-1)){
      pr[t+nyears-1,j] <- prod(sad[t:j])*prod(q[t:(j-1)])*p[j]
    } #j
    # Bajo la diagonal principal
    for (j in 1:(t-1)){
      pr[t+nyears-1,j] <- 0
    } #j
    # Última columna: probabilidad de no-recaptura
    pr[t+nyears-1,nyears] <- 1 - sum(pr[t+nyears-1,1:(nyears-1)])
  } #t

  # 3.3. Probabilidad para datos de productividad: regresión de Poisson
  for (t in 1:(nyears-1)){
    J[t] ~ dpois(rho[t])
    rho[t] <- R[t]*f[t]
  }
})

str (data  <-  list(m = m, 
                    y = y,
                    J = J))

str(constants<-list(R = R, 
                    nyears = dim(m)[2], 
                    r = rowSums(m)))

set.seed(1960)
N1s<-rep(2, 10)
Nads<-rep(35, 10)
Ntots<-N1s+Nads

str(inits  <-  list(mean.sjuv = runif(1, 0, 1),
                    mean.sad = runif(1, 0, 1),
                    mean.p = runif(1, 0, 1),
                    mean.fec = runif(1, 0, 10),
                    sigma.y = runif(1, 0, 1),
                    n1 = rpois(1, 30),
                    nad = rpois(1, 30),
                    N1=N1s,
                    Nad=Nads,
                    Ntot=Ntots,
                    mean1=rep(10,10),
                    lambda=runif(9,.5,1)))

params <- c("mean.sjuv", "mean.sad", "mean.p", "mean.fec", 
  "N1", "Nad", "Ntot","lambda")


# Preparamos el modelo para ejecución en Nimble
Rmodel <- nimbleModel(code=code, constants=constants, data=data, 
                      inits=inits, check=FALSE)
Rmodel$initializeInfo()
Cmodel <- compileNimble(Rmodel)
# Establecemos los parámetros a monitorizar
mcmcspec<-configureMCMC(Rmodel, monitors=params)

# Contruimos el modelo
IPM_MCMC <- buildMCMC(mcmcspec)

# Compilamos
CIPM_MCMC <- compileNimble(IPM_MCMC, project = Rmodel)

# Ejecutamos el modelo
nb=1000      # Iteraciones a desechar
ni=5000 +nb  # Iteraciones
nc=3         # Cadenas


start.time2<-Sys.time()
outNim <- runMCMC(CIPM_MCMC, niter = ni , nburnin = nb , nchains = nc, 
                  inits=inits, setSeed = TRUE, progressBar = TRUE, 
                  samplesAsCodaMCMC = TRUE)
end.time<-Sys.time()
end.time-start.time2 # tiempo de ejecución

summary(outNim)
xyplot(outNim[,1:10])

ipm<-as.data.frame(rbind(outNim$chain1,outNim$chain2,outNim$chain2))
Ntot<-ipm[,21:30]
head(Ntot)
lambda<-ipm[,31:39]
head(lambda)


# Vamos a plotear la evolución 
par(cex = 1.2)
meanNtot<-apply(Ntot, 2, mean)
lowerNt <- apply(Ntot, 2, quantile, p=0.025)
upperNt <- apply(Ntot, 2, quantile, p=0.975)

plot(meanNtot, type = "b", ylim = c(10, 55),
  ylab = "Tamaño de población",
  xlab = "Año",
  las = 1, pch = 16, col = "blue", frame = F, cex = 1.5)
segments(1:10, lowerNt, 1:10, upperNt, col = "blue")
points(y, type = "b", col = "black", pch = 16, lty = 2, cex = 1.5)
legend(x = 1, y = 20, legend = c("Conteos", "Predicciones"),
  pch = c(16, 16), col = c("black", "blue"), lty = c(2, 1), bty = "n")

meanlambda<-apply(lambda, 2, mean)
lowerlambda <- apply(lambda, 2, quantile, p=0.025)
upperlambda <- apply(lambda, 2, quantile, p=0.975)

plot(1:9, meanlambda, type = "l", ylim = c(0.5, 1.4),
  ylab = "Tasa de crecimiento de la población",
  xlab = "Año", las = 1,
  pch = 16, col = "black",
  frame = F, cex.lab=1.25)
polygon(x = c(1:9, 9:1), y = c(lowerlambda, upperlambda[9:1]),
  col = "gray90", border = "gray90")
points(1:9, meanlambda, type = "l", pch = 16, col = "black", cex = 1, lwd=2)
abline(h=1, lty=2)


sjuv=0.36; sad=0.93; f=0.48
A  <- matrix(c(
               0,     sad*f/2,
            sjuv,       sad), nrow = 2, byrow = TRUE)
lambda(A)


vr<-list(sjuv=sjuv, sad=sad, f=f)
# Creamos la matriz del modelo de población pre-reproducción:
sp.A <- expression(
               0,      sad*f/2,
            sjuv,       sad)


# Vemos la sensibilidad y elasticidad
x <- vitalsens(sp.A, vr)
x
barplot(t(x[,2:3]), beside=TRUE, legend=TRUE, las=1, 
  xlab="Ratios vitales",
  main="Sensibilidad y elasticidad de los ratios vitales")
abline(h=0)


# Vamos a ver que podría pasar con esa población, añadiendo la propia
# variación que ya presenta lambda, y suponiendo que esa variación sea
# aleatoria
lam<-c(0.931,0.9647,0.9771,1.041,1.0344,1.0424,1.0319,1.0352,0.9644)
meanlog<-mean(log(lam))
sdlog<-sd(lam)
hist(exp(rnorm(1000,meanlog,sdlog)))

# Múltiples simulaciones
set.seed(3)
sims <- 100
years<- 50


outmat<-matrix(NA,years,sims)       
outmat[1,]<-50

umbral<- 30
set.seed(1960)
for (i in 1:sims){
  for (t in 2:years){
    lam<-exp(rnorm(1,meanlog,sdlog))
    outmat[t,i]<-outmat[(t-1),i]*lam
    if(outmat[t,i]<=umbral) break
  }
}

matplot(1:years, outmat, type = "l", las = 1, ylab = "Tamaño de población",
    xlab = "Años")
abline(h = umbral, lty = 2)

time<-NULL # creamos un vector contenedor
for (i in 1:sims){ # bucle sobre las simulaciones
  t<-max(which(outmat[,i]>0))
  time<-c(time,t)
}
time.under<-time[which(time<50)]
(length(time.under)*100/sims)
