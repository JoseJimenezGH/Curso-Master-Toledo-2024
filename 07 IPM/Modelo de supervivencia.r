#==============================================================================#
#                                                                              #
#                        MODELOS DE SUPERVIVENCIA MULTIESTADO                  #
#                          José Jiménez (CSIC-IREC)                            #
#                                     UCLM                                     #
#                              15/03/2023 9:32:40                              #
#       Based on the book "Bayesian population analysis using WinBUGS          #
#                       a hierarchical perspective"                            #
#                     by Marc Kéry & Michael Schaub                            #
#                         (2012, Academic Press)                               #
#                                                                              #
#==============================================================================#

rm(list=ls())
library(tidyverse)
library(coda)
library(mcmcOutput)
library(nimble)
library(diagram)

set.seed(1982)
# Definir la función para simular los datos de captura-recaptura de 
# varios estados
simul.ms <- function(PSI.STATE, PSI.OBS, marked, unobservable = NA){
   # Inobservable: número que define el estado no observable
   n.occasions <- dim(PSI.STATE)[4] + 1
   CH <- CH.TRUE <- matrix(NA, ncol = n.occasions, nrow = sum(marked))
   # Definir un vector con la ocasión de marcaje
   mark.occ <- matrix(0, ncol = dim(PSI.STATE)[1], nrow = sum(marked))
   g <- colSums(marked)
   for (s in 1:dim(PSI.STATE)[1]){
      if (g[s]==0) next  # Para evitar el mensaje de error si no hay nada 
                         # que sustituir
      mark.occ[(cumsum(g[1:s])-g[s]+1)[s]:cumsum(g[1:s])[s],s] <-
      rep(1:n.occasions, marked[1:n.occasions,s])
      } #s
   for (i in 1:sum(marked)){
      for (s in 1:dim(PSI.STATE)[1]){
         if (mark.occ[i,s]==0) next
         first <- mark.occ[i,s]
         CH[i,first] <- s
         CH.TRUE[i,first] <- s
         } #s
      for (t in (first+1):n.occasions){
         # Procesos multinomiales para las transiciones de estado
         if (first==n.occasions) next
         state <- which(rmultinom(1, 1, PSI.STATE[CH.TRUE[i,t-1],,i,t-1])==1)
         CH.TRUE[i,t] <- state
         # Procesos multinomiales para el proceso de observación
         event <- which(rmultinom(1, 1, PSI.OBS[CH.TRUE[i,t],,i,t-1])==1)
         CH[i,t] <- event
         } #t
      } #i
   # Sustituir NA y el número de estado más alto (muerto) por 0
   CH[is.na(CH)] <- 0
   CH[CH==dim(PSI.STATE)[1]] <- 0
   CH[CH==unobservable] <- 0
   id <- numeric(0)
   for (i in 1:dim(CH)[1]){
      z <- min(which(CH[i,]!=0))
      ifelse(z==dim(CH)[2], id <- c(id,i), id <- c(id))
      }
   return(list(CH=CH[-id,], CH.TRUE=CH.TRUE[-id,]))
    # CH: historiales de captura que se utilizarán
    # CH.TRUE: historiales de captura con observación perfecta
   }


# Definimos los parámetros de la simulación
s <- 0.8
F <- 0.6
r <- 0.5
p <- 0.9
n.occasions <- 10   
n.states <- 4      # Numero de estados verdaderos simulados
n.obs <- 3         # Número de estados que observamos
marked <- matrix(0, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(40, n.occasions)	# individuos "liberados" en la zona
                                    # de estudio


# 1. Proceso de estados
totrel <- sum(marked)*(n.occasions-1)
PSI.STATE <- array(NA, dim=c(n.states, n.states, totrel, n.occasions-1))
for (i in 1:totrel){
   for (t in 1:(n.occasions-1)){
      PSI.STATE[,,i,t] <- matrix(c(
            s*F, s*(1-F), (1-s)*r,   (1-s)*(1-r),
            0,     s,     (1-s)*r,   (1-s)*(1-r),
            0,     0,           0,             1,
            0,     0,           0,             1), nrow = n.states, byrow = TRUE)
      } #t
   } #i



######################### DIAGRAMA DE ESTADOS ##################################

# Usamos la misma matriz:
A<-  matrix(c(
            s*F, s*(1-F), (1-s)*r,   (1-s)*(1-r),
            0,     s,     (1-s)*r,   (1-s)*(1-r),
            0,     0,           0,             1,
            0,     0,           0,             1), nrow = n.states, byrow = TRUE)

A.label <- (c('Vivo\nen el area\nnatal', 
             'Vivo\nfuera del area\nnatal', 
             'Recuperado\nmuerto', 
             'Muerto y no\n recuperado')) # Etiquetamos los estados

A<-round(A,3)
plotmat(t(A), 
        name=A.label,
        lwd = 1, 
        box.lwd = 2, 
        cex.txt = 0,  # 0.7
        box.cex = 0.75,
        box.size = 0.11, 
        box.col = c('light blue','khaki'),
        arr.length=.3,
        arr.width=.2,
        self.cex = .6,
        self.shiftx = -0.1, 
        relsize = 0.875,
        main="")

################################################################################

# 2.Proceso de observación
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
   for (t in 1:(n.occasions-1)){
      PSI.OBS[,,i,t] <- matrix(c(
            p,    0,    1-p,
            0,    0,      1,
            0,    1,      0,
            0,    0,      1), nrow = n.states, byrow = TRUE)
      } #t
   } #i


# Ejecutamos la simulación
sim <- simul.ms(PSI.STATE, PSI.OBS, marked)
CH <- sim$CH
zsimul <- sim$CH.TRUE

# Calculamos la fecha de la primera captura
get.first <- function(x) min(which(x!=0))
first <- apply(CH, 1, get.first)

# Recodificamos la matriz de observaciones ("0" no permitido)
rCH <- CH  CH re-codificada
rCH[rCH==0] <- 3


code <- nimbleCode({
  # -------------------------------------------------
  # Parámetros:
  # s: probabilidad de supervivencia
  # F: probabilidad de fidelidad a la zona de cría
  # r: probabilidad de recuperación
  # p: probabilidad de recaptura
  # -------------------------------------------------
  # Estados (S):
  # 1 vivo en el área de estudio
  # 2 vivo fuera de la zona de estudio
  # 3 recientemente muerto y recuperado
  # 4 muerto recientemente, pero no recuperado, o muerto (absorbente)
  # Observaciones (O):
  # 1 visto vivo
  # 2 recuperado muerto
  # 3 ni visto ni recuperado
  # -------------------------------------------------

  # Priors
  s ~ dunif(0, 1)     # supervivencia media
  F ~ dunif(0, 1)     # fidelidad media a la zona
  r ~ dunif(0, 1)     # probabilidad de recuperación
  p ~ dunif(0, 1)     # probabilidad de recaptura

  # Definir las matrices de transición de estado y de observación
 # Definir las probabilidades del estado S(t+1) dado S(t)
  ps[1,1] <- s*F
  ps[1,2] <- s*(1-F)
  ps[1,3] <- (1-s)*r
  ps[1,4] <- (1-s)*(1-r)
  
  ps[2,1] <- 0
  ps[2,2] <- s
  ps[2,3] <- (1-s)*r
  ps[2,4] <- (1-s)*(1-r)
  
  ps[3,1] <- 0
  ps[3,2] <- 0
  ps[3,3] <- 0
  ps[3,4] <- 1
  
  ps[4,1] <- 0
  ps[4,2] <- 0
  ps[4,3] <- 0
  ps[4,4] <- 1

  # Definir las probabilidades de O(t) dado S(t)
  po[1,1] <- p
  po[1,2] <- 0
  po[1,3] <- 1-p
  
  po[2,1] <- 0
  po[2,2] <- 0
  po[2,3] <- 1
  
  po[3,1] <- 0
  po[3,2] <- 1
  po[3,3] <- 0
  
  po[4,1] <- 0
  po[4,2] <- 0
  po[4,3] <- 1

  # Probabilidades
  for (i in 1:nind){
    # Definir el estado latente en la primera captura
    z[i,first[i]] <- y[i,first[i]]
    for (t in (first[i]+1):n.occasions){
      z[i,t] ~ dcat(ps[z[i,t-1], 1:4])
      # Proceso de observación: obtener O(t) dado S(t)
      y[i,t] ~ dcat(po[z[i,t], 1:3])
    } #t
  } #i
})


# Datos
str(data   <-    list(y = rCH))

# Constantes
str(constants <- list(first = first,
                      n.occasions = dim(rCH)[2],
                      nind = dim(rCH)[1]))

zinit <- zsimul
for (i in 1:nrow(zinit)){
   zinit[i, first[i]] <- NA
}

# Inicios
str(inits   <-   list(s = runif(1, 0, 1),
                      F = runif(1, 0, 1),
                      p = runif(1, 0, 1),
                      r = runif(1, 0, 1),
                      z = zinit))

# Parámetros a monitorizar
params <- c("s", "F", "r", "p")

Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits,
  check=FALSE, calculate=FALSE)
Rmodel$initializeInfo()
Cmodel <- compileNimble(Rmodel)
conf<-configureMCMC(Rmodel, monitors=params, onlySlice=TRUE, thin=5)

MCMC <- buildMCMC(conf)
CompMCMC<- compileNimble(MCMC, project = Rmodel)

outNim <- runMCMC(CompMCMC, niter = 25000 , nburnin = 1500, nchains = 3, 
                  inits=inits, setSeed = TRUE, progressBar = TRUE, 
                  samplesAsCodaMCMC = TRUE)

options(scipen=999)
MCMCsummary(outNim, digits=2)
diagPlot(outNim)