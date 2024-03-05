#==============================================================================#
#                                                                              #
#                        MODELOS DE POBLACIÓN MATRICIALES                      #
#                     SEGUIMIENTO DE LA DIVERSIDAD BIOLÓGICA                   #
#                          José Jiménez (CSIC-IREC)                            #
#                       UNIVERSIDAD DE CASTILLA-LA MANCHA                      #
#                             21/02/2023 18:46:07                              #
#                                                                              #
#==============================================================================#

# 1. BUHO MOTEADO DEL NORTE
#-------------------------------------------------------------------------------
# Vamos a crear la matriz de proyección. Se trata de una especie con tres 
# estados: juvenil (primer año), subadulto (segundo) y adulto (tres o más)-
# Vamos a usar el paquete popbio.

library(popbio)
library(diagram)
vr <- list(s0=0.108, s1=0.71, s2=0.942, f2=0.48)
stages <- c("juvenil", "subadulto", "adulto")
post <- expression( matrix2(c(
             0,       0,       s2*f2*0.5,
             s0,      0,       0,
             0,       s1,      s2), stages))
A <- eval(post, vr)

# Vemos la matriz
A


# CICLO VITAL
#==============
# Antes de seguir vamos a visualizar el ciclo vital de la especie:
A.label <- c("Juv", "Sub", "Ad") # Etiquetamos los estados
plotmat(round(A,2), pos=3.5,curve=0.5,name=A.label,lwd=2,arr.len=0.6,
        arr.width=0.25,my=-0.1, box.cex=1.00, relsize = 1.175,
        box.size=0.05,arr.type="triangle",dtext= 0.95,
        main="")

# PROYECCIÓN POBLACIONAL
# ========================
# Fijamos los valores de abundancia inicial y preparamos para hacer una 
# proyección a 50 años:
years <- 0:50
nYears <- length(years)
# Creamos una matriz contenedora con NA
n <- matrix(NA, nYears, 3)
n[1,] <- c(100, 50, 25)    ## Abundancia inicial

# Hacemos la proyección y lo visualizamos
for(t in 2:nYears) {
  n[t,] <- A %*% n[t-1,]   ## Multiplicación de matrices
}

# Vemos como en el plazo de 50 años, la población se reajusta
matplot(years, n, type="o", pch=16, xlab="Tiempo", ylab="Abundancia")
legend(35, 100, c("Juvenil", "Subadulto", "Adulto"),
       bty='n', col=1:3, pch=16, lty=1)


# DISTRIBUCIÓN ESTABLE DE EDADES
# ===============================
# Con el tiempo, las diferentes clases de edad se estabilizan. Obtenemos con 
# esta proyección la distribución estable de estados (o edades).
N <- rowSums(n) ## Abundancia total cada año
c <- n/N
matplot(years, c, type="o", pch=16, ylim=c(0, 1), xlab="Tiempo",
        ylab="Proporción de cada clase de edad")
legend(35, 1, c("Juvenil", "Subadulto", "Adulto"),
       bty='n', col=1:3, pch=16, lty=1)


# Calculamos la distribución estable de edades con la matriz de proyección
SAD.proj <- c[nYears,]
SAD.proj

# Ahora calculamos la distribución estable de edades usando los valores 
# propios:
eA <- eigen(A)
SADu <- Re(eA$vectors[,1])
SAD <- SADu/sum(SADu)
SAD

# Lo mismo podemos hacer usando popbio
stable.stage(A)


# TASA ANUAL DE CRECIMIENTO (LAMBDA)
# ===================================
lambda.it <- n[-1,]/n[-nYears,] ## Dividimos cada registro por el anterior
matplot(years[-1],## Hay nYears-1 *intervalos* entre años
       lambda.it, type="o", pch=16, xlab="Tiempo", ylab="Tasa de crecimiento")
abline(h=1, lty=2, col="grey") ## LÃ­nea horizontal
legend(35, 2.35, c("Juvenil", "Subadulto", "Adulto"),
       bty='n', col=1:3, pch=16, lty=1)



lambda.it[nYears-1,] ## Deben ser el mismo valor
lambda.proj <- lambda.it[nYears-1,1]
lambda.proj ## Ratio de crecimiento asintótico

# Ahora vamos a calcular el ratio de crecimiento usando valores propios
lambda <- Re(eA$values[1])
lambda  ## ratio asintótico de crecimiento

# Lo mismo podemos hacer usando popbio
lambda(A)

# VALOR REPRODUCTIVO POR CLASES DE EDAD
#========================================
# Ahora calculamos el valor reproductivo (RV) para cada clase de edad.
eAT <- eigen(t(A))
RVu <- Re(eAT$vectors[,1])
RV <- RVu/RVu[1]
RV

reproductive.value(A)

# ANÁLISIS DE SENSIBILIDAD Y ELASTICIDAD
#========================================
# Lo haremos directamente con popbio:
sensitivity(A)
elasticity(A)

elm <- expression( 
            0,       0,  s2*f2*0.5,
           s0,       0,          0,
            0,      s1,         s2)

x<-vitalsens(elm, vr)
x

# Vamos a verlo gráficamente
barplot(t(x[,2:3]), beside=TRUE, legend=FALSE, las=1, xlab="Ratios vitales",
  main="")
legend(1,0.9, legend=c("Sensitividad", "Elasticidad"),
       fill=c("black","grey"),
       cex=1.25,
       bty='n')
abline(h=0)

# 2. PANTERA DE FLORIDA
#-------------------------------------------------------------------------------
# Este es otro caso de modelo matricial estructurado en cuatro estados, basado 
# en otro caso real, el puma (Puma concolor coryi) o pantera de Florida.
vr <- list(s0=0.62, s1=0.70, s2=0.70, s3=0.75, f=1.5)
stages <- c("Cria", "Juv", "Subad", "Ad")

post <- expression( matrix2(c(
            0,     0,      s2*f*0.5,   s3*f*0.5,
            s0,    0,      0,          0,
            0,     s1,     0,          0,
            0,     0,      s2,         s3), stages ))

A <- eval(post, vr)
# Veamos la matriz
A
lambda(A)  # Vemos directamente la tasa de crecimiento poblacional


# CICLO VITAL
#==============
A.label <- c("Cria", "Juv", "Subad", "Ad")

plotmat(round(A,2), pos=4, curve=0.5,name=A.label,lwd=2,arr.len=0.6,
  arr.width=0.25,my=-0.2, box.cex=1.00,
  box.size=0.05,arr.type="triangle",dtext= 0.95,
  main="")


# PROYECCIÓN POBLACIONAL
# ========================
years <- 0:50
nYears <- length(years)
n <- matrix(NA, nYears, 4)
n[1,] <- c(50, 20, 10, 5)    ## Abundancia inicial

for(t in 2:nYears) {
  n[t,] <- A %*% n[t-1,]     ## Multiplicación de matrices
}

matplot(years, n, type="o", pch=16, xlab="Tiempo", ylab="Abundancia")
legend(40,50, legend=c("Cria", "Juv", "Subad", "Ad"),
   bty='n', col=c(4,1,2,3), pch=16, lty=1)


# ANÁLISIS DE SENSIBILIDAD Y ELASTICIDAD
#========================================
elm <- expression(
            0,     0,  s2*f*0.5,   s3*f*0.5,
           s0,     0,         0,          0,
            0,    s1,         0,          0,
            0,     0,        s2,         s3)

x<-vitalsens(elm, vr)
x


barplot(t(x[,2:3]), beside=TRUE, legend=FALSE, las=1, xlab="Ratios vitales",
  main="")
legend(1,0.9, legend=c("Sensitividad", "Elasticidad"),
       fill=c("black","grey"),
       cex=1.25,
       bty='n')
abline(h=0)


# Otra forma de evaluar las tasas vitales vs tasa de crecimiento poblacional
n<-length(vr)
vr<-seq(0.7,1,.05)  # Secuencia de valores de los ratios vitales
# hacemos una matriz de almacenamiento de los valores de lambda
vrsen<-matrix(NA, nrow=length(vr), ncol=n, dimnames=list(vr, names(vr)))

# Este bucle funciona con un ratio cada vez
for (h in 1:n){
  vr2<-list(s0=0.62, s1=0.70, s2=0.70, s3=0.75, f=1.5)
  for (i in 1:length(vr)){
    vr2[[h]]<-vr[i]
    A<-matrix(sapply(elm, eval,vr2 , NULL), nrow=sqrt(length(elm)), byrow=TRUE)
    vrsen[i,h] <- max(Re(eigen(A)$values))
  }
}

matplot(rownames(vrsen), vrsen, type='l', lty=1, lwd=2, las=1, col=1:5,
 ylab="Crecimiento de la población de la pantera de Florida", 
 xlab="Valores de los ratios vitales",
 main="Efectos del cambio de los ratios vitales")


# Construimos la leyenda
vrn<-expression(s0, s1, s2, s3, f)
legend(0.7, 1.165, vrn, lty=1, lwd=2, col=1:5, cex=1.2, bty='n')


