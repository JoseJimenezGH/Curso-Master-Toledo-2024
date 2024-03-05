#==============================================================================#
#                                                                              #
#                              PVA utilizando IPM                              #
#                   Medidas de actuación y planes de recuperación              #
#                          José Jiménez (CSIC-IREC)                            #
#                         UNIVERSIDAD DE CASTILLA-LA MANCHA                    #
#                             27/02/2024 9:01:34                               #
#                                                                              #
#==============================================================================#
# El conjunto de datos del IPM que vamos a ver se ha simulado con 10 años de 
# datos, con una población inicial de 400 jóvenes y 200 adultos, con una 
# supervivencia anual de 0.36 y 0.93, y una productividad de 0.24 hembras/hembra 
# reproductora a partir del segundo año de vida. Se ha usado una probabilidad 
# de captura del 0.95 para jóvenes y 0.8 para adultos, y una proabilidad de 
# encontrar una pareja reproductora, y registrar el éxito reproductor, del 0.8.
# Para simular los datos se ha usado el código de Abadi et al. (2010) adaptado 
# por Schaub.

setwd('C:/Users/Administrator/OneDrive/66 Master Seguimiento de la Diversidad Biológica/Lab/07 PVA-IPM/')

library(nimble)
library(popbio)
library(lattice)
library(coda)

set.seed(1960)
# Simulamos la supervivencia juvenil, pero añadiendo un ruido gaussiano 
# estocástico alrededor de la media de la supervivencia juvenil
mean.sjuv<-0.36
mean.logit.sjuv <- logit(mean.sjuv)
eps.sjuv<-rnorm(9,0,0.5)
sjuv <- expit(mean.logit.sjuv + eps.sjuv)

sjuv<-c(0.6050644,0.5740094,0.5180878,0.3774764,0.4163225,0.4685812,0.2780022,
        0.2332290,0.5274793)
        
# Hacemos lo mismo con la supervivencia adulta
mean.sad<-0.93
mean.logit.sad <- logit(mean.sad)
eps.sad<-rnorm(9,0,0.5)
sad <- expit(mean.logit.sad + eps.sad)

sad<-c(0.9585836,0.9092798,0.9734925,0.9504287,0.9341231,0.9641143,0.9360352,
       0.9417192,0.9077353)

f11<-0
f12<-0.24 

vr<-list(mean.sjuv=0.36, mean.sad=0.93, f12=0.48)
stages <- c("juvenil","adulto")
post <- expression( matrix2(c(
                  0,     mean.sad*f12/2, # Aqui se usan hembras (50%)
          mean.sjuv,           mean.sad), stages ))
(A <- eval(post, vr))
lambda(A)

# Conteos de machos adultos (desde el año 1 al 10)
y <- c(28,32,25,38,40,43,36,40,44,49)

# Datos de productividad (desde el año 1 al 9)
J <- c(2,22,16,16,20,24,36,16,18,32)   # Pollos volados
R <- c(42,42,39,44,50,48,48,51,52,55)  # Nidos localizados


# Datos de captura-recaptura (en formato m-array, desde el año 1 a 10)
m<-matrix(c(3,  1,  0,  0,  0,  0,  0,  0,  0,  5,
            0,  7,  3,  0,  0,  0,  0,  0,  0,  9,
            0,  0,  8,  2,  1,  1,  0,  0,  0,  5,
            0,  0,  0,  6,  0,  2,  0,  0,  0, 22,
            0,  0,  0,  0,  4,  0,  1,  0,  0, 12,
            0,  0,  0,  0,  0, 10,  0,  0,  0, 11,
            0,  0,  0,  0,  0,  0,  1,  2,  0, 17,
            0,  0,  0,  0,  0,  0,  0,  5,  0, 22,
            0,  0,  0,  0,  0,  0,  0,  0, 12, 20,
           44, 10,  1,  0,  0,  0,  0,  0,  0,  3,
            0, 39, 11,  2,  0,  0,  0,  0,  0,  8,
            0,  0, 49,  4,  1,  1,  0,  0,  0,  4,
            0,  0,  0, 57, 11,  2,  0,  0,  0,  4,
            0,  0,  0,  0, 55, 10,  1,  0,  0,  6,
            0,  0,  0,  0,  0, 44, 15,  5,  1,  7,
            0,  0,  0,  0,  0,  0, 52,  8,  2,  8,
            0,  0,  0,  0,  0,  0,  0, 50, 13,  7,
            0,  0,  0,  0,  0,  0,  0,  0, 50, 20), 
            ncol = 10, byrow = TRUE)

code <- nimbleCode({

#-------------------------------------------------
# MODELO DE POBLACIÓN INTEGRADO 
#  - Modelo estructurado por edades con 2 clases:
#    joven (hasta 1 año) y adulto (al menos 2 años)  
#  - Edad de primera reproducción: 1 año
#  - Conteos pre-reproductores (machos territoriales)
#  - Todos los ratios vitales se asumen constantes
#-------------------------------------------------

  ##############################################################################
  #                           1. PRIORS
  ##############################################################################
  # Error de observación
  tauy <- pow(sigma.y, -2)
  sigma.y ~ dunif(0, 50)
  sigma2.y <- pow(sigma.y, 2)

  # Tamaños iniciales de población
  n1 ~ T(dnorm(5, tauy),0,1000)      # 1-año
  n2 ~ T(dnorm(50, tauy),0,1000)     # Adultos
  N[1,1,1] <- round(n1)
  N[2,1,1] <- round(n2)

  # Supervivencia, productividad y probabilidad de recaptura
  # Priors
 	mean.logit.sjuv <- logit(mean.sjuv)
	mean.sjuv ~ dunif(0, 1) 
	mean.logit.sad <- logit(mean.sad)
	mean.sad ~ dunif(0, 1)
  
 	sigma.sjuv ~ dunif(0, 10)
	tau.sjuv <- pow(sigma.sjuv, -2)
	sigma.sad ~ dunif(0, 10)
	tau.sad <- pow(sigma.sad, -2)
  
 	mean.log.f <- log(mean.f)
	mean.f ~ dunif(0, 10)
	sigma.f ~ dunif(0, 10)
	tau.f <- pow(sigma.f, -2)

  mean.p ~ dunif(0, 1) 
 	
  for (t in 1:(nyears-1)){
		p[t] <- mean.p
	}

  ## 1. CONTROL: SIN GESTIÓN (OPCION NULA)
  #================================================ 
  for (t in 1:(nyears-1+K)){
		# Supervivencia juvenil
    logit.sjuv[t] <- mean.logit.sjuv + eps.sjuv[t]
		eps.sjuv[t] ~ dnorm(0, tau.sjuv)
		sjuv[t] <- expit(logit.sjuv[t])
		
    logit.sad[t,1] <- mean.logit.sad + eps.sad[t,1]
		eps.sad[t,1] ~ dnorm(0, tau.sad)
		sad[t,1] <- expit(logit.sad[t,1])
  }
  for (t in 1:(nyears+K)){
		log.f[t,1] <- mean.log.f + eps.f[t,1]
		eps.f[t,1] ~ dnorm(0, tau.f)
		f[t,1] <- exp(log.f[t,1])
  }
  
  # 2: INCREMENTAR LA PRODUCTIVIDAD EN UN 25%
  #================================================
	# Hasta inicio de la gestión
	for (t in 1:nyears){      
		log.f[t,2] <- log.f[t,1]
		eps.f[t,2] <- eps.f[t,1]
		f[t,2] <- f[t,1]
	}
	
	# Futura: incremento de un 25%
	for (t in (nyears+1):(nyears+K)){
		log.f[t,2] <- mean.log.f + log(1.25) + eps.f[t,2]
		eps.f[t,2] ~ dnorm(0, tau.f)
		f[t,2] <- exp(log.f[t,2])
	}

  # 3. REDUCCION VARIABILIDAD SUPERVIVENCIA ADULTA
  #================================================
 	# Hasta inicio de la gestión
 	for (t in 1:(nyears-1)){   
    logit.sad[t,2] <- logit.sad[t,1]
		eps.sad[t,2] <- eps.sad[t,1]
		sad[t,2] <- sad[t,1]
 	}
 	
 	# Proyección a futuro
 	for (t in nyears:(nyears-1+K)){    
    logit.sad[t,2] <- mean.logit.sad + eps.sad[t,2]
		eps.sad[t,2] ~ dnorm(0, tau.sad*2)   # Incrementamos la precisión
		sad[t,2] <- ilogit(logit.sad[t,2])
 	}
  
  # 4. SUELTA DE 5 HEMBRAS CADA AÑO DE GESTIÓN
  #=================================================
  # No afectamos a la supervivencia ni a la fecundidad, sino al
  # tamaño poblacional (ver debajo)




  ##############################################################################
  #                           2. PROBABILIDAD
  ##############################################################################

  # 2.1. Probabilidad de los datos de conteos
    # 2.1.1 Proceso de sistema (realidad biológica)
    # 1. CONTROL (SIN GESTIÓN)
    #=================================================
    for (t in 1:nyears-1+K){
      N[1,t+1,1] ~ dpois(f[t,1] * sjuv[t] * (N[1,t,1] + N[2,t,1]))
      N[2,t+1,1] ~ dbin(sad[t,1], (N[1,t,1] + N[2,t,1]))
    }
    
    # 2. CON INCREMENTO DE PRODUCTIVIDAD
    #=================================================
  	# Hasta inicio de la gestión
    for (t in 1:nyears){
  		N[1,t,2] <- N[1,t,1]
  		N[2,t,2] <- N[2,t,1]
  	}
  	# Proyección a futuro
  	for (t in nyears:(nyears-1+K)){
  		N[1,t+1,2] ~ dpois(f[t,2] * sjuv[t] * (N[1,t,2] + N[2,t,2]))
  		N[2,t+1,2] ~ dbin(sad[t,1], (N[1,t,2] + N[2,t,2]))
  	}
  	
 		# 3. REDUCCION VARIABILIDAD SUPERVIVENCIA ADULTA
 		#=================================================
  	# Hasta inicio de la gestión
  	for (t in 1:nyears){
  		N[1,t,3] <- N[1,t,1]
  		N[2,t,3] <- N[2,t,1]
  	}
  	# Proyección a futuro
  	for (t in nyears:(nyears-1+K)){
  		N[1,t+1,3] ~ dpois(f[t,1] * sjuv[t] * (N[1,t,3] + N[2,t,3]))
  		N[2,t+1,3] ~ dbin(sad[t,2], (N[1,t,3] + N[2,t,3]))
  	}
  	
 		# 4. SUELTA DE 5 HEMBRAS CADA AÑO DE GESTIÓN
 		#=================================================
  	# Hasta inicio de la gestión
  	for (t in 1:nyears){
  		N[1,t,4] <- N[1,t,1]
  		N[2,t,4] <- N[2,t,1]
  	}
  	# Proyección a futuro
  	for (t in nyears:(nyears-1+K)){
  		N[1,t+1,4] ~ dpois(f[t,1] * sjuv[t] * (N[1,t,4] + N[2,t,4] + 5))
  		N[2,t+1,4] ~ dbin(sad[t,1], (N[1,t,4] + N[2,t,4] + 5))
  	}

    # 3.1.2 Proceso de observación
    for (t in 1:nyears){
      y[t] ~ dnorm(N[1,t,1] + N[2,t,1], tauy)
    }

  # 2.2 Probabilidad de los datos de captura-recaptura:    
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
      pr[t,j] <- sjuv[t]*prod(sad[(t+1):j,1])*prod(q[t:(j-1)])*p[j]
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
    pr[t+nyears-1,t] <- sad[t,1] * p[t]
    # Por encima de la diagonal principal
    for (j in (t+1):(nyears-1)){
      pr[t+nyears-1,j] <- prod(sad[t:j,1])*prod(q[t:(j-1)])*p[j]
    } #j
    # Bajo la diagonal principal
    for (j in 1:(t-1)){
      pr[t+nyears-1,j] <- 0
    } #j
    # Última columna: probabilidad de no-recaptura
    pr[t+nyears-1,nyears] <- 1 - sum(pr[t+nyears-1,1:(nyears-1)])
  } #t

  # 2.3. Probabilidad para datos de productividad: regresión de Poisson
  for (t in 1:(nyears-1)){
    J[t] ~ dpois(rho[t])
    rho[t] <- R[t]*f[t,1]
  }
 	# Derived parameters
	for (t in 1:(nyears+K)){
		Ntot[t,1] <- N[1,t,1] + N[2,t,1] # Sin gestión
		Ntot[t,2] <- N[1,t,2] + N[2,t,2] # Gestión de productividad
		Ntot[t,3] <- N[1,t,3] + N[2,t,3] # Disminuir variabilidad sad
		Ntot[t,4] <- N[1,t,4] + N[2,t,4] # Suelta de 5 hembras anuales
	}

})

str (data  <-  list(m = m, 
                    y = y,
                    J = J))

str(constants<-list(R = R, 
                    nyears = dim(m)[2], 
                    r = rowSums(m),
                    K = 5))


nyears<-10; K<-5

Ns<- array(NA,c(2,15,4))
for(i in 1:4){
  Ns[,,i]<-matrix(rbind(rep(5,15),rep(30,15)))
}

eps.ads<-matrix(runif(15*2*2,0,1),ncol=2)
eps.sjuvs<-c(runif(15*1*2,0,1))
eps.fs<-cbind(runif((nyears+K),0.3,1),runif((nyears+K),0.3,1))

str(inits  <-  list(mean.sjuv = runif(1, 0.5, 1),
                    sigma.sjuv=0.1,
                    tau.sjuv=0.1,
                    mean.sad = runif(1, 0.5, 1),
                    sigma.sad=0.1,
                    tau.sad=0.1,
                    sigma.f=0.1,
                    tau.f=0.1,
                    eps.sjuv=eps.sjuvs,
                    eps.sad=eps.ads,
                    eps.f=eps.fs,
                    mean.p = runif(1, 0.3, 1),
                    mean.f = runif(1, 0.3, 5),
                    sigma.y = runif(1, 0.2, 1),
                    n1 = rpois(1, 5),
                    n2 = rpois(1, 30),
                    N=round(Ns,0),
                    Ntot=apply(Ns,c(2,3),sum)))

params <- c('mean.sjuv','sjuv','mean.sad','sad','mean.p','mean.f', 'f',
            'N','Ntot')


# Preparamos el modelo para ejecución en Nimble
Rmodel <- nimbleModel(code=code, constants=constants, data=data, 
                      inits=inits, check=FALSE)
Rmodel$initializeInfo()

Cmodel <- compileNimble(Rmodel)
# Establecemos los parámetros a monitorizar
mcmcspec<-configureMCMC(Rmodel, monitors=params, nthin=10)

# Contruimos el modelo
IPM_MCMC <- buildMCMC(mcmcspec)

# Compilamos
CIPM_MCMC <- compileNimble(IPM_MCMC, project = Rmodel)

# Ejecutamos el modelo
nb=5000       # Iteraciones a desechar
ni=50000 +nb  # Iteraciones
nc=3          # Cadenas


start.time2<-Sys.time()
outNim <- runMCMC(CIPM_MCMC, niter = ni , nburnin = nb , nchains = nc, 
                  inits=inits, setSeed = TRUE, progressBar = TRUE, 
                  samplesAsCodaMCMC = TRUE)
end.time<-Sys.time()
end.time-start.time2 # tiempo de ejecución

summary(outNim)
xyplot(outNim[,10:21])

ipm<-as.data.frame(rbind(outNim$chain1,outNim$chain2,outNim$chain2))
names(ipm) <- gsub('\"', "", names(ipm))

# Probabilidad de que incrementar la productividad sea mejor que el control
round(mean(ipm$'Ntot[15, 2]' >  ipm$'Ntot[15, 1]'), 3)
round(sd(ipm$'Ntot[15, 2]' >  ipm$'Ntot[15, 1]'), 3)

# Probabilidad de reducir la variabilidad sea mejor que el control
round(mean(ipm$'Ntot[15, 3]' >  ipm$'Ntot[15, 1]'), 3)
round(sd(ipm$'Ntot[15, 3]' >  ipm$'Ntot[15, 1]'), 3)

# Probabilidad de que gestion con sueltas sea mejor que el control
round(mean(ipm$'Ntot[15, 4]' >  ipm$'Ntot[15, 1]'), 3)
round(sd(ipm$'Ntot[15, 4]' >  ipm$'Ntot[15, 1]'), 3)

# Es evidente que la cuarta opción es la que va a tener un mejor resultado en 
# el incremento de la población, tanto considerando el valor medio como la
# variación