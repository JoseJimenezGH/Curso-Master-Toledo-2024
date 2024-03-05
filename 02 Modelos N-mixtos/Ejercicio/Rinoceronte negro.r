#==============================================================================#
#                                                                              #
#                             Modelos N-mixtos                                 #
#                         Jose Jimenez. CSIC-IREC                              #
#                          09/03/2023 17:54:46                                 #
#                            MASTER UCLM 2023                                  #
#                                                                              #
#                                                                              #
#==============================================================================#

# El caso de estudio es un modelo N-mixto abierto (plurianual) con los datos
# del artículo:

# Kidwai, Z., J. Jiménez, C. J. Louw, H. P. Nel, and J. P. Marshal. 2019.
# Using N-mixture models to estimate abundance and temporal trends of black
# rhinoceros (Diceros bicornis L.) populations from aerial counts. Global
# Ecology and Conservation 19:e00687.

# La selección de modelos que se presenta es una versión simplificada del
# original, ya que solo se hace selección de covariables en la detección para
# abreviar los cálculos,.

setwd('C:/Users/Usuario/OneDrive/00 Master Seguimiento de la Diversidad Biológica 2023/AulaVirtual/02 Modelos N-mixtos/Ejercicio')
library(unmarked)

load('yy.RData')

# Estudiemos los datos:
yy

# Preparamos los datos en formato utilizable por unmarked
ymat<-matrix(yy,nrow=2, ncol=17*3)

# Estima naïve: vamos a utilizar el conteo más alto de cada año (es una
# simplificación que se usa con frecuencia, pero no es correcta, ya que no
# resuelve el problema de la detección). Extraemos el máximo anual para cada
# sitio:
y<-apply(yy,c(1,3),max)

plot(1999:2015,y[1,],pch=1, ylim=c(0,120), cex=1.25, type='b',
  xlab="Year", ylab="Naive population estimate", col="red")
points(1999:2015, y[2,], cex=1.25, type='b', col="blue")
legend(2009.5,18, lty=1, lwd=c(1,1), pt.cex = 1.25, pch=c(1,1), col=c('red','blue'),
  cex=1.25, legend=c('Madikwe','Pilanesberg'), bty='n')

# Preparamos los datos para su uso en unmarked:
Ft<-data.matrix(read.table('flyTimeRhino.txt', header=FALSE))
Jd<-data.matrix(read.table('julianDate.txt', header=FALSE))
Rm<-data.matrix(read.table('removal.txt', header=FALSE))

# Covariable categórica de sitio
sc <- data.frame(Site = c('A','B'))  # A: Madikwe; B: Pilanesberg

# Covariable específica de la observacion:
oc <- list(
    Obs = matrix(1:51, nrow=2, ncol=51, byrow=TRUE),
    Jd = Jd)       # Julian date
# Covariables de sitio y año
ysc <- list(
    Year = matrix(1:17, nrow=2, ncol=17, byrow=TRUE),
    Ft = Ft,       # Hora de inicio del vuelo
    Rm=Rm)         # Extracciones

umf <- unmarkedFramePCO(y=ymat,
                        siteCovs=sc,
                        obsCovs=oc,
                        numPrimary=17,
                        yearlySiteCovs=ysc)

# Estandarizamos covariables
yearlySiteCovs(umf) <- scale(yearlySiteCovs(umf))
obsCovs(umf) <- scale(obsCovs(umf))

# Vemos los datos
umf
summary(umf)

# Ajustamos el modelo (versión simplificada del original del artículo)
#=======================================================================
# lambda: siteCovs only
# omega/gamma: siteCovs o yearlySiteCovs
# p: siteCovs, yearlySiteCovs, o obsCovs

# Ajustamos un modelo N-mixto abierto: lambda (= abundancia),
#                                      gamma (= ganancia o recruitment),
#                                      omega (= supervivencia),
#                                      p (= detección)

fm.0 <- pcountOpen(~1, ~1, ~1, ~1, data=umf,
  mixture='P', K=130, control=list(trace=TRUE), start=c(4,3.1,1,1)); fm.0

(lam <- exp(coef(fm.0, type='lambda')))
(gam <- exp(coef(fm.0, type='gamma')))
(om <- plogis(coef(fm.0, type='omega')))
(p <- plogis(coef(fm.0, type='det')))

fm.1 <- pcountOpen(~1, ~1, ~1, ~Site, umf,
  start=c(4,3.5,0.2,-0.3,0.8), K=130)
fm.2 <- pcountOpen(~1, ~1, ~1, ~Site+Jd+I(Jd^2), umf,
  start=c(4,3.5,0.2,-0.3,0.8,0,0), K=130)
fm.3 <- pcountOpen(~1, ~1, ~1, ~Site+Jd+I(Jd^2)+Ft, umf,
  start=c(4,3.5,0.2,-0.3,0.8,0,0,0), K=130)
fm.4 <- pcountOpen(~1, ~1, ~1, ~Site+Jd+I(Jd^2)+Ft+Year, umf,
  start=c(4,3.5,0.2,-0.3,0.8,0,0,0,0), K=130)
fm.5 <- pcountOpen(~1, ~1, ~1, ~Site+Jd+I(Jd^2)+Ft+Year+Obs, umf,
  start=c(4,3.5,0.2,-0.3,0.8,0,0,0,0,0), K=130)

brm1<- fitList('lam(.)gam(.)om(.)p(.)'                             = fm.0,
               'lam(.)gam(.)om(.)p(Site)'                          = fm.1,
               'lam(.)gam(.)om(.)p(Site+Jd+I(Jd^2))'               = fm.2,
               'lam(.)gam(.)om(.)p(Site+Jd+I(Jd^2)+Ft)'            = fm.3,
               'lam(.)gam(.)om(.)p(Site+Jd+I(Jd^2)+Ft+Year)'       = fm.4,
               'lam(.)gam(.)om(.)p(Site+Jd+I(Jd^2)+Ft+Year+Obs)'   = fm.5)
modSel(brm1)
##                                                 nPars    AIC delta   AICwt cumltvWt
## lam(.)gam(.)om(.)p(Site+Jd+I(Jd^2))                 7 655.01  0.00 4.9e-01     0.49
## lam(.)gam(.)om(.)p(Site+Jd+I(Jd^2)+Ft)              8 656.05  1.04 2.9e-01     0.78
## lam(.)gam(.)om(.)p(Site+Jd+I(Jd^2)+Ft+Year)         9 657.78  2.77 1.2e-01     0.91
## lam(.)gam(.)om(.)p(Site+Jd+I(Jd^2)+Ft+Year+Obs)    10 658.37  3.36 9.2e-02     1.00
## lam(.)gam(.)om(.)p(Site)                            5 668.08 13.07 7.1e-04     1.00
## lam(.)gam(.)om(.)p(.)                               4 690.74 35.73 8.6e-09     1.00

(lam <- exp(coef(fm.2, type='lambda')))     # lam es la abundancia inicial
(gam <- exp(coef(fm.2, type='gamma')))      # gam es el reclutamiento
(om <- plogis(coef(fm.2, type='omega')))    # om es la supervivencia aparente
(p <- plogis(coef(fm.2, type='det')))       # p es la probabilidad de detección

black.rhino<-ranef(fm.2)  # Para obtener las estimas por sitio y año
bup(black.rhino)

# Obtenemos los intervalos de confianza
low.Mad<-confint(black.rhino)[1,1,]
hig.Mad<-confint(black.rhino)[1,2,]
low.Pil<-confint(black.rhino)[2,1,]
hig.Pil<-confint(black.rhino)[2,2,]

# Ploteamos los resultados
date<-c(1999:2015)
points(date, bup(black.rhino)[1,], col='red', type='p', pch=16, cex=1.25)
segments(date,low.Mad, date, hig.Mad, col='red', lwd=1)
points(date+0.1, bup(black.rhino)[2,], col='blue', type='p', pch=16, cex=1.25)
segments(date+0.1,low.Pil, date+0.1, hig.Pil, col='blue', lwd=1)

# Predicción de abundancia inicial
newdat <- data.frame(Site =c('A','B'), Jd=0)
(pop<-predict(fm.2, type='det', newdata=newdat))

# Probabilidad de detección de cada sitio
par(mar=c(5,5,3,3)-0.1, las=1 )
plot(c(0.5,1.5), pop$Predicted, pch=16, cex=3,
  xlim=c(0,2), ylim=c(0,1),
  xaxt='n',
  ylab='Detection probabilty',
  xlab='Reserve', cex.lab=1.25)
segments(c(0.5,1.5),pop$lower,c(0.5,1.5),pop$upper, lwd=3)
v1<-c(0.5, 1.5)
v2<-c('Madikwe','Pilanesberg')
axis(side = 1, at = v1, labels = v2, las=1)

# Predecimos la detección con la fecha juliana
Jd.ss <- umf@obsCovs$Jd
# Pilanesberg
newdat2 <- data.frame(Jd = seq(min(Jd.ss),max(Jd.ss),,100), Site="B")
pred2 <- predict(fm.2, type='det', newdata=newdat2, append = T)
Jd.new<-seq(min(Jd),max(Jd),,100)
plot(Jd.new, pred2$Predicted,
  main="Pilanesberg",
  type='l', ylim=c(0,1), lwd=3,
  xlab='Julian date',
  ylab='Detection')
lines(Jd.new, pred2$lower, lty=2)
lines(Jd.new, pred2$upper, lty=2)

