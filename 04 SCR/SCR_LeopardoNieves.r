#==============================================================================#
#                                                                              #
#                             CAPTURA-RECAPTURA SCR                            #
#                            LEOPARDO DE LAS NIEVES                            #
#                                 Jose Jimenez                                 #
#                             07/02/2021 10:52:32                              #
#                                                                              #
#==============================================================================#

rm(list=ls(all=TRUE))
setwd('C:/Leopard')  # poned aquí vuestro path
source("SCR_functions.R")
library(scrbook)
library(calibrate)
library(aspace)
library(nimble)
library(coda)
library(lattice)

# Si quereis más información sobre los datos, podeis visitar la página
# browseURL("https://osf.io/hr6xm/")

X1<-read.table("traps.txt", header=TRUE)
dim(X1)
# Veamos el contenido
X1

# Vamos a ver la ubicación de la trampas:
xlims2<-c(min(X1[,2]),max(X1[,2]))+c(-1000,1000)
ylims2<-c(min(X1[,3]),max(X1[,3]))+c(-1000,1000)


# Veamos la distancia entre las trampas más cercanas:
leopardo<-  X1[,2:3]
output2 <- matrix(ncol=nrow(leopardo), nrow=nrow(leopardo))
for (i in 1: nrow(leopardo)){
     singlePoint <- c(leopardo[i,1], leopardo[i,2])
     output2[i,] <- c(distances(centre.xy=singlePoint, destmat=leopardo, verbose=FALSE))
     		  }

# Convierto en NA la diagonal
output2[output2=="0"]<-NA
DM<-apply(output2, 2, min, na.rm=TRUE); length(DM)

# Distancia media entre trampas
mean(DM, na.omit=TRUE)

# Rango de distancias
range(DM)

# ¿Es adecuada la distancia entre trampas?... Lo veremos más tarde

# Lo ploteamos usando asp=1 para que ajuste la relación x/y
plot(X1[,2:3],pch="+",col="blue", asp=1, xlim=xlims2, ylim=ylims2)
# Vamos a etiquetar las cámaras trampas:
textxy(X1[,2], X1[,3], seq(1,37), m=c(0.5,0.5), cex=1)


# Escalamos y centramos las coordenadas, ya que en general los programas que
# vamos a utilizar trabajan mejor alrededor del cero
X1 <- X1[,2:3]/1000
names(X1)<-c('X','Y')
# Estandarizamos
newtrapsx<-X1[,1]-mean(X1[,1])
newtrapsy<-X1[,2]-mean(X1[,2])
X<-cbind(newtrapsx,newtrapsy)
n.llx=min(X[,1])
n.upx=max(X[,1])
n.lly=min(X[,2])
n.upy=max(X[,2])
# Vemos el resultado
head(X)

# OPERATIVIDAD DE LAS CÁMARAS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Veamos la operatividad de las cámaras (cuando están en funcionamiento)
MASK<-data.matrix(read.table("MASK.txt", header=FALSE))
J<-dim(MASK)[1] # trampas
K<-dim(MASK)[2] # dias de muestreo

## Tenemos que convertir esto en matriz
x<-as.matrix(MASK)  
image(1:K,1:J,t(x), yaxt = "n", xlab="Occasion", ylab="",
  col=topo.colors(2),, cex.lab=1.25)
mtext(side = 2, "Camera", line = 2.5, cex=1.25)
axis(2, rev(seq(1, 43, by=1)))
box()

# Y apilamos los dias de operatividad de cada cámara
KT<-apply(MASK,1,sum)
colnames(MASK)<-c(1:K)

# Vamos a leer las capturas para ponerlas en un formato utilizable, y 
# organizarlas para su uso en el codigo SCR:
wtraps<-cbind(1:J,X,MASK)
wcaps<-as.matrix(read.table("ID.txt", header=TRUE))
head(wcaps)

# Vamos a ver si hay algún duplicado
wcaps[duplicated(wcaps),]

# Ahora usamos una función que a partir de las trampas y las capturas crea 
# una matriz tridimensional (individuo-trampa-ocasión)
datYknown <-SCR23darray(wcaps, wtraps)

# ¿Cuantos individuos hemos visto?
nind<-dim(datYknown)[1]; nind

yr<-apply(datYknown, c(2,3),sum) # Esto lo uso para comprobar que la máscara es correcta
sum(yr)
sum(yr*MASK) # Está OK!!

## Aumentado de datos
M<-100   # Ponemos más del triple de lo que aproximadamente pensamos 
         # que podría haber.
nind<- dim(datYknown)[1]
Yaug <- array(0, dim = c(M, J, K))
Yaug[1:nind, , ] <- datYknown
y<-apply(Yaug, c(1,2), sum)

# Vamos a plotear las capturas y la ubicación promedio de los diferentes 
# individuos y los extremos donde se ha localizado (spiderplot)

xlim<-c(n.llx,n.upx)
ylim<-c(n.lly,n.upy)
plot(X, xlim=xlim+c(-2,2), ylim=ylim+c(-2,2), pch="+", asp=TRUE)
datn<-apply(datYknown, c(2,3), sum)
tot<-apply(datn, 1,sum)
symbols(X, circles=tot, inches=F, bg="#00000022", fg=NULL, add=T)
points(X, pch="+", cex=0.5)
# Spiderplot
spiderplotJJ4(datYknown, X, buffer=2, lwd=2)


# DEFINICIÓN DEL ESPACIO DE ESTADOS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
buff<- 20             # el buffer debe ser tal que no detectamos ningun animal
                      # más allá del buffer (al menos 2.5*sigma)
xlim<-xlim+c(-buff,buff)
ylim<-ylim+c(-buff,buff)
area<-diff(xlim)*diff(ylim)


# DEFINIMOS EL MODELO EN BUGS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Escribimos el modelo
code <- nimbleCode({

  alpha0 ~ dnorm(0,.1)
  logit(p0) <- alpha0
  alpha1 ~ dnorm(0,.01)
  sigma <- sqrt(1/(2*alpha1))
  psi ~ dunif(0,1)
  
  for(i in 1:M){
    z[i] ~ dbern(psi) 
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    d[i,1:J] <- pow(s[i,1]-X[1:J,1],2) + pow(s[i,2]-X[1:J,2],2)
    p[i,1:J] <- p0*z[i]*exp(- alpha1*d[i,1:J])
    
    for(j in 1:J){   
      y[i,j] ~ dbin(p[i,j],KT[j])
    }
  }
  N <- sum(z[1:M])
  D <- 100*N/area    # individuos/100 km2
})


# CONSTANTES
#~~~~~~~~~~~~~~~
str(constants <- list(M = M, K=K, J=J, area=area))

# DATOS
#~~~~~~~~~~~~~~~
str(data <- list (y=y, X=X, xlim=xlim, ylim=ylim, KT=KT))

# Generamos inicios
z <- c(rep(1,nind),rep(0,M-nind))
sst <- cbind(runif(M,xlim[1],xlim[2]),runif(M,ylim[1],ylim[2]))
for(i in 1:nind){
  sst[i,1] <- mean( X[y[i,]>0,1] )
  sst[i,2] <- mean( X[y[i,]>0,2] )
}

# INICIOS
#~~~~~~~~~~~~~~~
str(inits <- list (p0=0.5, alpha1=1, s=sst, z=z))


# COMPILACIÓN
#~~~~~~~~~~~~~~~
Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits, check=FALSE)

Cmodel <- compileNimble(Rmodel)
params<- c('N', 'D', 'sigma','psi', 'p0', 's','z')
mcmcspec<-configureMCMC(Rmodel, monitors=params)

pumpMCMC <- buildMCMC(mcmcspec)

CpumpMCMC <- compileNimble(pumpMCMC, project = Rmodel)

# Ejecutamos el modelo
nb=1000      # Iteraciones a desechar
ni=5000 +nb  # Iteraciones
nc=3         # Cadenas


start.time2<-Sys.time()
outNim <- runMCMC(CpumpMCMC, niter = ni , nburnin = nb , nchains = nc, inits=inits,
                  setSeed = TRUE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)
end.time<-Sys.time()
end.time-start.time2 # tiempo de ejecución

## OBTENEMOS EL OUTPUT:
## ====================
summary(outNim[,c('N', 'D', 'sigma','psi', 'p0')])

# Inspección de la convergencia
xyplot(outNim[,c('N', 'D', 'sigma','psi', 'p0')])

# Veamos el histograma de densidad poblacional, a ver si el aumentado de datos 
# era suficiente
hist(as.matrix(outNim[,'D']), main="", xlab="D", ylab="Frecuencia")

# PLOTEADO DE CENTROS DE ACTIVIDAD
samplesn<-as.matrix(outNim)
s1 <- samplesn[,c(5:104)] ; Sx<-as.matrix(s1)
s2 <- samplesn[,c(105:204)] ; Sy<-as.matrix(s2)
z <- samplesn[,c(206:305)] ; z<- as.matrix(z)
delta<-10
Xl<-min(X[,1])-delta
Xu<-max(X[,1])+delta
obj<-list(Sx=Sx,Sy=Sy,z=z)
JJ<-SCRdensity(obj, nx=50, ny=50, scalein=1, scaleout = 100,
  Xl=xlim[1], Xu=xlim[2], Yl=ylim[1], Yu=ylim[2])

# Añadimos las trampas (escaladas)
points(X, pch="+", col=1, lwd=1, cex=1.5)