#==============================================================================#
#                                                                              #
#                        MUESTREO DE DISTANCIAS JERÁRQUICO                     #
#                     Seguimiento de la Diversidad Biológica                   #
#                          José Jiménez (CSIC-IREC)                            #
#                      UNIVERSIDAD DE CASTILLA-LA MANCHA                       #
#                             27/02/2024 9:01:34                               #
#                                                                              #
#==============================================================================#

# Premisas del muestreo de distancias
#====================================== 
# - Los individuos en la línea de observación son detectados con probabilidad 1 
# - Los animales deben ser detectados antes de efectuar algún movimiento de 
#   reacción frente al observador
# - Cada detección debe ser independiente
# - No debe existir error de medición
# 
# 
# Simulador para crear datos
#=============================
# La mejor forma de ver si un código en R funciona, es usar simulaciones. 
# Creamos una población de parámetros conocidos, y vemos si la ejecución del 
# código nos devuelve la población de partida

set.seed(1975)
dist.sim<-function(Dha=1, sigma=250, wmax=500, L=1000, nL=10){
  Dm2<-Dha/10000
  lambda<-Dm2*L*wmax*2
  n<-rpois(nL,lambda)
  # Función de detección seminormal
  p.detect<-function(x,sigma){
    det.p<-exp(-(x^2/(2*sigma^2)))
  }
  dist.data<-data.frame(line=numeric(0),length=numeric(0),distance=numeric(0))
  for(i in 1:nL) {
    all.dist<-runif(n[i],0,wmax)
    detections<-rbinom(length(all.dist),1,p.detect(all.dist,sigma=sigma))
    obs.dist<-all.dist[detections==1]
    line<-rep(i,length(obs.dist))
    length<-rep(L,length(obs.dist))
    new.data<-data.frame(line=line,length=L,distance=obs.dist)
    dist.data<-rbind(dist.data,new.data)
  }
  return<-list(dist.data,n)
}


# Creación de datos
#====================
# Vamos a seleccionar una densidad determinada, el valor de sigma 
# (parámetro de la distribución seminormal, que nos informa de como decrece la 
# detectabilidad con la distancia), longitud de los transectos (todos serán 
# iguales) y número de transectos a realizar.
set.seed(1995) 

D<-0.15
sim.data<-dist.sim(Dha=D, sigma=250, wmax=500, L=1000, nL=10) 
# Esto equivale a S=500*2*1000*10 (1000 ha; 0.15*1000=150 animales)
dat<-sim.data[[1]]                                                               

# Disponemos la información para tratarla en unmarked
breaks<-seq(0,500,50)
n<-length(breaks)-1
labels<-as.numeric(1:n)

# Agrupamos los datos de distancia por intervalos y calculamos las frecuencias
# para cada linea
sim.data[[1]]$binned.x<-cut(sim.data[[1]]$distance,breaks=breaks,labels=labels)
y<-with(sim.data[[1]],table(line,binned.x))

# Convertimos en matriz
class(y) <- "matrix"
# Vamos a ver la matriz
y

# Lista de longitudes de transecto
tlength<-with(sim.data[[1]],aggregate(length~line,FUN="mean"))$length
tlength


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                              #
#                     APROXIMACIÓN FRECUENTISTA CON UNMARKED                   #
#                                                                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Cargamos unmarked
library(unmarked)

# Preparamos los datos para tratarlos en _unmarked_
sim.data.input<-unmarkedFrameDS(y=y, dist.breaks=breaks, tlength=tlength, 
  survey="line",unitsIn="m")
# Vemos el marco de datos creado:
head(sim.data.input)


# Ajuste de los modelos competidores
mod1<-distsamp(~1 ~1, data=sim.data.input, keyfun="halfnorm",output="density",
  unitsOut="ha")
mod2<-distsamp(~1 ~1, data=sim.data.input, keyfun="hazard", output="density",
  unitsOut="ha")
mod3<-distsamp(~1 ~1, data=sim.data.input, keyfun="exp", output="density",
  unitsOut="ha")
mod4<-distsamp(~1 ~1, data=sim.data.input, keyfun="uniform", output="density",
  unitsOut="ha")  

# Seleccionamos el modelo de mejor AIC
fmList <- fitList(Halfnorm=mod1, Hazard=mod2, Exponential=mod3, Uniform=mod4)
modSel(fmList)

# El mejor modelo es mod1
summary(mod1)

# Veamos el histograma
hist(mod1, xlab="Distancia (m)", col="lightgrey", lwd=4)
x<-sim.data[[1]]$distance
rug(x,side=1)

# Transformamos a escala real
backTransform(mod1, type="state")   # animales / ha
backTransform(mod1, type="det")

re1 <- ranef(mod1, K=100)

# Función auxiliar
postPredN <- function(ranefOut, nSims=100) {
    post <- ranefOut@post
    postDim <- dim(post)
    nSites <- postDim[1]
    nProbs <- postDim[2]
    possibleN <- 0:(nProbs-1)
    nYears <- postDim[3]
    NsitePostSamples <- array(NA_integer_,
                              c(nSites, nYears, nSims))

    for(i in 1:nSites) {
        for(t in 1:nYears) {
            NsitePostSamples[i,t,] <- sample(
                possibleN, size=nSims,
                replace=TRUE, prob=post[i,,t])
        }
    }
    NpostSamples <- apply(NsitePostSamples, c(2,3), sum)
    return(NpostSamples)
}

# Obtenemos media, mediana, SD y CI95%
Npost <- postPredN(re1, nSims=1000)
rowMeans(Npost)             ## Posterior mean of annual abundance
apply(Npost, 1, median)     ## Posterior median
apply(Npost, 1, sd)
apply(Npost, 1, quantile,
      prob=c(0.025, 0.975)) ## 95% CI



# Veamos la bondad del ajuste (GOF)
(fm <- distsamp(~1 ~1, keyfun="halfnorm", data=sim.data.input))

# Bondad del ajuste

# Vamos a ver el GOF con *parboot* (boostrapping)
# Esta función nos devuelve tres estadísticos.
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

(pb <- parboot(mod1, fitstats, nsim=1000, report=1))
# chisquare test statistic
(c.hat <- pb@t0[2] / mean(pb@t.star[,2])) # c-hat as ratio of observed/expected

par(mfrow=c(2,2))
hist(pb@t.star[,1], xlab='SSE', col="grey90", breaks=15, main="")
abline(v=pb@t0[1],lty=2, lwd=3)
hist(pb@t.star[,2], xlab='Chisq', col="grey90", breaks=15, main="")
abline(v=pb@t0[2],lty=2, lwd=3)
hist(pb@t.star[,3], xlab='freemanTukey', col="grey90", breaks=15, main="")
abline(v=pb@t0[3],lty=2, lwd=3)

# Podemos estimamos la abundancia en cada transecto
re <- ranef(mod1, K=50)
# Best Unbiased Predictors
media<-bup(re, stat="mean", real)
ci<-confint(re, level=0.9) # 90% CI
(comp<-data.frame('media'=media,'lower'=ci[,1],'upper'=ci[,2],'real'=sim.data[[2]]))

# Real
cat("Poblacion simulada = ", sum(comp[,4]), "individuos", "\n")
# Estimado
cat("Poblacion estimada = ", round(sum(comp[,1]),0),"(",sum(comp[,2]),"-",sum(comp[,3]),") individuos\n") 



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                              #
#                           APROXIMACIÓN BAYESIANA                             #
#                                                                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Vamos a escalar los datos de distancia
dat$distance<-dat$distance/100
dat$y<-1
dat<-data.frame(site=dat$line,y=dat$y,d=dat$distance)
B<-5

# Vamos a obtener el número de individuos detectados por sitio
# ncap = 1 más el número de individuos detectados por sitio
ncap <- table(dat[,1])             # ncap = 1 si no hay observaciones
sites0 <- dat[is.na(dat[,2]),][,1] # sitios sin detecciones
ncap[as.character(sites0)] <- 0    # Rellenamos con 0 los sitios sin
                                   #  detecciones
ncap <- as.vector(ncap)

# Preparamos los datos
site <- dat[!is.na(dat[,2]),1]    # ID de sitio por observacion
delta <- 0.1                      # agrupamos distancias en franjas
midpt <- seq(delta/2, B, delta)   # puntos medios de las franjas
dclass <- dat[,3] %/% delta + 1   # convertimos distancias en grupos
nD <- length(midpt)               # Número de intervalos de distancia
dclass <- dclass[!is.na(dat[,2])] # Categorias observadas
nind <- length(dclass)            # Número total de individuos observados

# Preparamos los datos
str(data <- list(midpt=midpt, ncap=ncap, dclass=dclass))
str(constants<-list(nsites=10, site=site, nind=nind, B=B, nD=nD, delta=delta))

library(nimble)
## define the model
code <- nimbleCode({
  # Priors
  alpha0 ~ dunif(-10,10)
  beta0 ~ dunif(-10,10)
  
  for(i in 1:nind){
    dclass[i] ~ dcat(fc[site[i],1:nD]) # Parte 1 del HDS
  }
  
  for(s in 1:nsites){
    # Construímos las probabilidades por celdas
    for(g in 1:nD){             # midpt = punto central de cada celda
      log(p[s,g]) <- -midpt[g] * midpt[g] / (2*sigma[s]*sigma[s])
      pi[s,g] <- delta / B      # probabilidad por intervalo
      f[s,g] <- p[s,g] * pi[s,g]
      fc[s,g] <- f[s,g] / pcap[s]
    }
    pcap[s] <- sum(f[s,1:nD])          # Pr(captura): suma de las areas
    
    ncap[s] ~ dbin(pcap[s], N[s])  # Parte 2 del HDS
    N[s] ~ dpois(lambda[s])        # Parte 3 del HDS
    log(lambda[s]) <- beta0
    log(sigma[s])<- alpha0
  }
  # Parámetros derivados
  Ntotal <- sum(N[1:nsites])
  area<- nsites*1*2*B  # Unidad long.== 1, franja = 2xB
  D<- Ntotal/(10*area) # animales por 100 ha
})


# Inicios
Nst <- ncap + 1
str(inits <- list(alpha0=0, beta0=0, N=Nst))

# Parámetros a salvar
params <- c('alpha0', 'beta0', 'Ntotal','D')


# Preparamos el modelo para ejecución en Nimble
Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits, check=FALSE)
Cmodel <- compileNimble(Rmodel)
# Establecemos los parámetros a monitorizar
mcmcspec<-configureMCMC(Rmodel, monitors=params, nthin=5)

# Contruimos el modelo
hdsMCMC <- buildMCMC(mcmcspec)

# Compilamos
ChdsMCMC <- compileNimble(hdsMCMC, project = Rmodel)

# Ejecutamos el modelo
nb=5000      # Iteraciones a desechar
ni=25000 +nb  # Iteraciones
nc=3         # Cadenas


start.time2<-Sys.time()
outNim <- runMCMC(ChdsMCMC, niter = ni , nburnin = nb , 
                  nchains = nc, inits=inits,
                  setSeed = TRUE, progressBar = TRUE, 
                  samplesAsCodaMCMC = TRUE)  
end.time<-Sys.time()
end.time-start.time2 # tiempo de ejecución

# Resultados
summary(outNim)

# Inspeccionamos la convergencia de las cadenas de Markov
xyplot(outNim)

