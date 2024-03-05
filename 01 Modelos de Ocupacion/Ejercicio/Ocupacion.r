#==============================================================================#
#                                                                              #
#                        Modelos de Ocupación                                  #
#                       Jose Jimenez. CSIC-IREC                                #
#                        09/03/2023 17:54:46                                   #
#                          MASTER UCLM 2023                                    #
#                                                                              #
#                                                                              #
#==============================================================================#
setwd('C:/Users/Usuario/OneDrive/00 Master Seguimiento de la Diversidad Biológica 2023/AulaVirtual/01 Modelos de Ocupacion/Ejercicio')

y<-read.table("DatosPresencia.txt", header=TRUE)


# Estima naïve
sum(apply(y,1,max,na.rm=TRUE))/dim(y)[1]

### Definicion del modelo
# Preparamos los datos para el analisis
str( win.data <- list(y = y, M = nrow(y), J = ncol(y)) )
# [1] 0.57

#  Especificamos el modelo en BUGS
cat(file = "model.txt",
"
model {
  # Informacion a priori
  psi ~ dunif(0, 1)
  p ~ dunif(0, 1)
  # Probabilidad
  for (i in 1:M) {             # Bucle sobre los sitios
    z[i] ~ dbern(psi)          # Modelo de estado
    for (j in 1:J) {           # Bucle sobre los muestreos replicados
      y[i,j] ~ dbern(z[i]*p)   # Modelo de observacion
    }
  }
}

")


#  Valores de inicio
zst <- apply(y, 1, max, na.rm=TRUE)     # Para evitar conflictos
                                        # entre datos/modelo/inicio
zst[zst=="-Inf"]<-0
inits <- function(){list(z = zst)}

#  Parametros a monitorizar
params <- c("psi", "p")

#  Configuracion MCMC
ni <- 15000   ;   nt <- 1   ;   nb <- 1000   ;   nc <- 3

### Ajuste del modelo

#  Ejecutamos JAGS
library(jagsUI)
start.time2<-Sys.time()
fm2 <- jags(win.data, inits, params, "model.txt", n.chains = nc,
   n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE)
end.time<-Sys.time()
end.time-start.time2 # tiempo de ejecución

#  Vemos resultados
print(fm2, dig = 3)

library(lattice)
xyplot(fm2$samples)
densityplot(fm2$samples)


# UTILIZANDO UNMARKED:
#=====================
library(unmarked)

umf <- unmarkedFrameOccu(y=y, siteCovs=NULL, 
    obsCovs=NULL)       # organizamos los datos
umf
plot(umf, panels=1)
summary(umf)            # resumen de los datos      
fm <- occu(~1 ~1, umf)  # ajustamos el modelo

# Para extraer los valores a escala real, con sus IC:
(btlp<-backTransform(fm, "det"))
confint(btlp, level = 0.95)
(btls<-backTransform(fm, "state"))
confint(btls, level = 0.95)


# Test de bondad del ajuste
fitstats <- function(fm) {
    observed <- getY(fm@data)
    expected <- fitted(fm)
    resids <- residuals(fm)
    sse <- sum(resids^2)
    chisq <- sum((observed - expected)^2 / expected)
    freeTuke <- sum((sqrt(observed) - sqrt(expected))^2)
    out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
    return(out)
    }

(pb <- parboot(fm, fitstats, nsim=1000, report=1))
plot(pb, main="")
