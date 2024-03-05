#######################################################################
#
# This code create data to be used in integrated population models
# - data are completely independent
# - data sets simulated are survey, capture-recapture, and fecundity
# Original in Abadi et al. 2010, Ecology 91: 7-14 - adapted by M.Schaub
#
#######################################################################
set.seed(1945)
####################
# Input parameters
####################

# Number of years
ti <- 10
ni <- ti-1

# Size of the data sets: number of individuals involved for capture-recapture, survey, and productivity
nd <- c(300, 150, 150)

# Initial population size of each age classes
Ni <- c(400, 200)

# Age specific survival rates (Juv, adult)
phi1 <- 0.36
phi2 <- 0.93

# Fecundity rate (females)
fl1 <- 0          # female productivity of one year old females
fl2 <- 0.24          # female productivity of females older than one year

# Capture and Recapture probabilities
pjuv <- 0.95     # this is the initial capture prob of juveniles
padults <- 0.80  # this is the initial capture as well as the recapture probability of adults

# Prob. to record breeding pair (survey)
psur <- 0.80

# Prob. to record the reproductive success
pprod <- 0.80

#######################################

# Combine the values to vectors

# Simulamos la supervivencia juvenil, pero añadiendo un ruido gaussiano 
# estocástico alrededor de la media de la supervivencia juvenil
mean.sjuv<-0.36
mean.logit.sjuv <- logit(mean.sjuv)
eps.sjuv<-rnorm(ti,0,0.5)
sjuv <- expit(mean.logit.sjuv + eps.sjuv); mean(sjuv)
PHI1 <- sjuv

# Hacemos lo mismo con la supervivencia adulta
mean.sad<-0.93
mean.logit.sad <- logit(mean.sad)
eps.sad<-rnorm(ti,0,0.5)
sad <- expit(mean.logit.sad + eps.sad); mean(sad)
PHI2 <- sad

PJUV <- rep(pjuv, ti)
PADULTS <- rep(padults, ti)
PSUR <- rep(psur, ti)
PPROD <- rep(pprod, ti)
FL1 <- rep(fl1, ti+1)
FL2 <- rep(fl2, ti+1)


###########################
# Population simulation
##########################

# 1. Matrix model to estimate roughly the population size with a Leslie matrix
N <- matrix(data = NA, nrow = length(Ni), ncol = ti+1)
N[,1] <- Ni
for (i in 1:ti){
  les <- matrix(data = c(PHI1[i]*FL1[i], PHI1[i]*FL2[i], PHI2[i], PHI2[i]), byrow=T, nrow = length(Ni), ncol = length(Ni))
  N[,(i+1)]<-les%*%N[,i]
}
no.anim <- sum(N)
no.ani <- round(no.anim*5)

# 2. Define array for each individual
ind <- array(data = NA, dim = c(5, ti+1, no.ani))   # infn about [1ye, adu, dead, rep, juv]

# 3. Simulate the fates of initial individuals

# 1-Year old individuals
for (i in 1:Ni[1]){
  if(Ni[1]==0) break
  ind[1,1,i] <- 1
  ind[4,1,i] <- rpois(1, FL1[1])
  j <- 1
  z1 <- rbinom(1, 1, PHI2[j])           # survived or not
  if(z1==0){ind[3,j+1,i] <- 1
  next}                                 # if dead
  ind[2,j+1,i] <- 1
  ind[4,j+1,i] <- rpois(1, FL2[j+1])    # no.of fledglings produced
  for (j in 2:ti){
    z2 <- rbinom(1, 1, PHI2[j])
    if(z2==0){
    ind[3,j+1,i] <- 1
    break}
    ind[2,j+1,i] <- 1
    ind[4,j+1,i] <- rpois(1, FL2[j+1])
  }
}

# Adults
for (i in 1:Ni[2]){
  if(Ni[2]==0) break
  m <- Ni[1]+i
  ind[2,1,m] <- 1
  ind[4,1,m] <- rpois(1, FL2[1])
  for (j in 1:ti){
    z<-rbinom(1, 1, PHI2[j])
    if(z==0){
      ind[3,j+1,m] <- 1
      break}
    ind[2,j+1,m] <- 1
    ind[4,j+1,m] <- rpois(1, FL2[j+1])
  }
}


# 4. Simulate the fates of newborn
juv <- rep(0,ti)
for (g in 1:ti){
  m <- sum(Ni)+sum(juv)
  juv[g] <- sum(ind[4,g,], na.rm = T)
  if(juv[g]==0) next
  for (i in 1:juv[g]){
    ind[5,g,m+i] <- 1
    # first time step
    z1 <- rbinom(1, 1, PHI1[g])
    if(z1==0){ind[3,g+1,m+i] <- 1
      next}
    ind[1,g+1,m+i] <- 1
    ind[4,g+1,m+i] <- rpois(1, FL1[g+1])
    if(g>=ti) next                       # if cycle too long
    # second time step
    z2 <- rbinom(1, 1, PHI2[g+1])
    if(z2==0){ind[3,g+2,m+i] <- 1
      next}
    ind[2,g+2,m+i] <- 1
    ind[4,g+2,m+i] <- rpois(1, FL2[g+2])
    # third time step and above
    if(g>=(ti-1)) next
    for (k in (g+2):ti){
      z3 <- rbinom(1, 1, PHI2[k])
      if(z3==0){ind[3,k+1,m+i] <- 1
        break}
      ind[2,k+1,m+i] <- 1
      ind[4,k+1,m+i] <- rpois(1, FL2[k+1])
    }
  }
}


# 5. Total number of animals
Ntotal <- sum(Ni)+sum(ind[4,1:ti,], na.rm = T)

# Reshape the array, remove empty cells
IND <- ind[,,1:Ntotal]
rownames(IND) <- c("1-Year", "Adu", "Dead", "Rep", "Juv")

###################################################
#  Create completely independent samples
###################################################

# Three independent samples
Ntot <- dim(IND)[3]
xt1 <- matrix(data = seq(1:Ntot), ncol = 1)
resamp1 <- resamp2 <- resamp3 <- numeric()

# Sample 1
resamp1 <- sample(xt1, nd[1], replace = F)
# Sample 2
resamp2 <- sample(xt1[-resamp1,], nd[2], replace = F)
# Sample 3
resamp3 <- sample(xt1[c(-resamp1, -resamp2),], nd[3], replace = F)

# Individuals selected for capture-recapture, survey and reproductive success data
INDC <- IND[,,resamp1]
INDB <- IND[,,resamp2]
INDR <- IND[,,resamp3]


#################################
# Create Capture-Recapture Data
#################################

# 1. Construct capture histories
CR <- matrix(data = rep(0, nd[1]*(ti+3)), ncol = ti+3, nrow = nd[1])
for ( i in 1:nd[1]){
  for (r in 1:ti){
  # Juveniles
  if(!is.na(INDC[5,r,i])){
    y1 <- rbinom(1, 1, PJUV[r])
    if(y1==1){
      CR[i,r] <- 1
      CR[i,ti+1] <- 1
    }
  }
  # 1-year old
  if(!is.na(INDC[1,r,i])){
    y2 <- rbinom(1, 1, PADULTS[r])
      if(y2==1){CR[i,r] <- 2
      if(CR[i,ti+1]==0) CR[i,ti+2] <- 1
    }
  }
  # Adults
  if(!is.na(INDC[2,r,i])){
    y3 <- rbinom(1, 1, PADULTS[r])
      if(y3==1){CR[i,r] <- 3
      if(CR[i,ti+1]==0&CR[i,ti+2]==0) CR[i,ti+3] <- 1
      }
    }
  }
}
colnames(CR)<-c(seq(1,ti), "Juveniles", "1-year", "Adults")


# 2. Create m-array for the capture-histories

# 2.1. Two functions for data manipulation:
# 2.1.1. Function to remove individuals without histories
clean <- function(crma){
  su <- rep(NA, nrow(crma))
  for (i in 1:nrow(crma)){
    su[i] <- sum(crma[i,])}
  final <- crma[which(su!=0),]
}

# 2.1.2. Function to remove all first captures
rmfirst <- function(crma){
  newma <- matrix(data = 0, ncol = ncol(crma), nrow = nrow(crma))
  for (i in 1:nrow(crma)){
    if (sum(crma[i,])==1) next  # if the individual is captured only once, move to next individual
    k <- min(which(crma[i,]==1))
    newma[i,] <- crma[i,]
    newma[i,k] <- 0
  }
  final <- clean(newma)
  if (is.null(ncol(final))==TRUE) {final <- matrix(final, ncol = length(final), nrow = 1)}
  final <- final
}

# 2.2. Recode capture-recapture matrix
cr2 <- replace(CR,which(CR[,]==2),1)
cr3 <- replace(cr2,which(cr2[,]==3),1)
cr <- clean(cr3)
m.juv <- matrix(data = 0, ncol = ti, nrow = ti)  # juv m-array
m.adu <- matrix(data = 0, ncol = ti, nrow = ti)  # 1y + adults m-array
cr.juv <- cr[which(cr[,(ti+1)]==1),1:ti]         # indivs captured as juv
cr.adu <- cr[c(which(cr[,(ti+2)]==1), which(cr[,(ti+3)]==1)), 1:ti]  # indivs captured as 1y + adults

# 2.3. Calculate m-array for juveniles
for (tt in 1:(ti-1)){
  cr.h <- matrix(data = 0, ncol = ti, nrow = nrow(cr.juv))
  if(sum(cr.juv[,tt])==0)next
  for (i in 1:nrow(cr.juv)){
    if (cr.juv[i,tt]==1){
      cr.h[i,] <- cr.juv[i,]
      cr.juv[i,] <- 0
    }
  }

  m.juv[tt,1] <- sum(cr.h[,tt])

  # Remove all first captures in cr.h
  if(sum(cr.h)==0) next
  # cr.h1<-rmfirst(clean(cr.h))
  cr.h1 <- clean(cr.h)
  cr.h1 <- matrix(cr.h1, ncol = (ti))
  cr.h1 <- rmfirst(cr.h1)
  if(sum(cr.h1)==0) next
  # Determine when the birds were first recaptured
  po<-rep(0, nrow(cr.h1))
  for (j in 1:nrow(cr.h1)){
    po[j] <- min(which(cr.h1[j,]==1))
  }

  k <- as.data.frame(table(po))    #  a table of time first recaptured (po) by #indivs (freq)
  if(length(k$Freq)>=1){
    for (i in 1:length(k$Freq)){
      m.juv[tt,as.numeric(as.vector(k[i,1]))] <- k[i,2]   #give the m-array of juveniles
    }
  }
  cr.adu <- rbind(cr.adu, cr.h1[which(po>=(tt+1)),])  #  merge adults cap hist and juvs recaptured as adults(1y & adults) following first release
}

# 2.4. m-array for adults
for (tt in 1:(ti-1)){
  cr.h <- matrix(data = 0, ncol = ti, nrow = nrow(cr.adu))
  if(sum(cr.adu[,tt])==0) next
  for (i in 1:nrow(cr.adu)){
    if (cr.adu[i,tt]==1){
      cr.h[i,] <- cr.adu[i,]
      cr.adu[i,] <- 0
    }
  }
  m.adu[tt,1] <- sum(cr.h[,tt])
  # Remove all first captures in cr.h
  if(sum(cr.h)==0) next
  # cr.h1 <- rmfirst(clean(cr.h))
  cr.h1 <- clean(cr.h)
  cr.h1 <- matrix(cr.h1, ncol = (ti))
  cr.h1 <- rmfirst(cr.h1)
  if(sum(cr.h1)==0) next
  # Determine when the birds were first recaptured
  po <- rep(0,nrow(cr.h1))
  for (j in 1:nrow(cr.h1)){
    po[j] <- min(which(cr.h1[j,]==1))
  }
  k <- as.data.frame(table(po))
  if(length(k$Freq)>=1){
    for (i in 1:length(k$Freq)){
      m.adu[tt,as.numeric(as.vector(k[i,1]))] <- k[i,2]
    }
  }
  cr.adu <- rbind(cr.adu, cr.h1)
}

# 2.5. Rearrange the m-array
m <- matrix(NA, ncol = ti, nrow = 2*(ti-1))
# Juv
for (i in 1:(ti-1)){
  for (j in 1:(ti-1)){
    m[i,j] <- m.juv[i,j+1]
    m[i,ti] <- m.juv[i,1]-sum(m.juv[i,2:ti])
  }
}

# Adults
for (i in 1:(ti-1)){
  for (j in 1:(ti-1)){
    m[i+ti-1,j] <- m.adu[i,j+1]
    m[i+ti-1,ti] <- m.adu[i,1]-sum(m.adu[i,2:ti])
  }
}

######################################################
# Create population survey data
######################################################
SUR <- rep(0,ti)
for (i in 1:nd[2]){
  for (u in 1:ti){
    if(!is.na(INDB[1,u,i])){
      y3 <- rbinom(1, 1, PSUR[u])
      if(y3==1){
        SUR[u] <- SUR[u]+1
      }
    }
    if(!is.na(INDB[2,u,i])){
      y3 <- rbinom(1, 1, PSUR[u])
      if(y3==1){
        SUR[u] <- SUR[u]+1
      }
    }
  }
}


#############################################
# Create reproductive success data
#############################################

R <- rep(0, ti) # number of pairs whose productivity was observed
nestlings <- rep(0, ti)  # total number of nestlings recoded in a year
for (i in 1:nd[3]){
  for (v in 1:ti){
    if(!is.na(INDR[4,v,i])){
        y <- rbinom(1, 1, PPROD[v])
      if(y==1){
        R[v] <- R[v]+1
        nestlings[v] <- nestlings[v]+INDR[4,v,i]
      }
    }
  }
}

#######################
# Bundle data
#######################

data <- list(m = m, y = SUR, J = nestlings*2, R = R, nyears = dim(m)[2], r = rowSums(m))



# 11.3.2. Analysis of the model
# Specify model in BUGS language
sink("ipm.jags")
cat("
model {
  #-------------------------------------------------
  #  Integrated population model
  #  - Age structured model with 2 age classes:
  #		1-year old and adult (at least 2 years old)
  #  - Age at first breeding = 1 year
  #  - Prebreeding census, female-based
  #  - All vital rates assumed to be constant
  #-------------------------------------------------

  #-------------------------------------------------
  # 1. Define the priors for the parameters
  #-------------------------------------------------
  # Observation error
  tauy <- pow(sigma.y, -2)
  sigma.y ~ dunif(0, 50)
  sigma2.y <- pow(sigma.y, 2)

  # Initial population sizes
  # Note that this part is different than in BUGS:
  # 1. JAGS seems to be very sensitive to the choice of the prior distribution for the initial population sizes and of the choice of the observation model (priors)
  # 2. Since the initial population sizes are used in binomial distributions, the numbers must be integers, otherwise JAGS does not run.
  # The following specification seems to work very well and produces similar results as BUGS:

  n1 ~ dnorm(3, tauy)T(0,)     # 1-year
  nad ~ dnorm(35, tauy)T(0,)    # Adults
  N1[1] <- round(n1)
  Nad[1] <- round(nad)

  # Survival and recapture probabilities, as well as productivity
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
  # 2. Derived parameters
  #-------------------------------------------------
  # Population growth rate
  #for (t in 1:(nyears-1)){
  #  lambda[t] <- Ntot[t+1] / Ntot[t]
  #}

  #-------------------------------------------------
  # 3. The likelihoods of the single data sets
  #-------------------------------------------------
  # 3.1. Likelihood for population population count data (state-space model)
    # 3.1.1 System process
    for (t in 2:nyears){
      mean1[t] <- f[t-1] / 2 * sjuv[t-1] * Ntot[t-1]
      N1[t] ~ dpois(mean1[t])
      Nad[t] ~ dbin(sad[t-1], Ntot[t-1])
    }
    for (t in 1:nyears){
      Ntot[t] <- Nad[t] + N1[t]
    }

    # 3.1.2 Observation process
    for (t in 1:nyears){
      y[t] ~ dnorm(Ntot[t], tauy)
    }

  # 3.2 Likelihood for capture-recapture data: CJS model (2 age classes)
  # Multinomial likelihood
  for (t in 1:2*(nyears-1)){
    m[t,1:nyears] ~ dmulti(pr[t,], r[t])
  }

  # m-array cell probabilities for juveniles
  for (t in 1:(nyears-1)){
    # Main diagonal
    q[t] <- 1-p[t]
    pr[t,t] <- sjuv[t] * p[t]
    # Above main diagonal
    for (j in (t+1):(nyears-1)){
      pr[t,j] <- sjuv[t]*prod(sad[(t+1):j])*prod(q[t:(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
      pr[t,j] <- 0
    } #j
    # Last column: probability of non-recapture
    pr[t,nyears] <- 1-sum(pr[t,1:(nyears-1)])
  } #t

  # m-array cell probabilities for adults
  for (t in 1:(nyears-1)){
    # Main diagonal
    pr[t+nyears-1,t] <- sad[t] * p[t]
    # Above main diagonal
    for (j in (t+1):(nyears-1)){
      pr[t+nyears-1,j] <- prod(sad[t:j])*prod(q[t:(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
      pr[t+nyears-1,j] <- 0
    } #j
    # Last column
    pr[t+nyears-1,nyears] <- 1 - sum(pr[t+nyears-1,1:(nyears-1)])
  } #t

  # 3.3. Likelihood for productivity data: Poisson regression
  for (t in 1:(nyears-1)){
    J[t] ~ dpois(rho[t])
    rho[t] <- R[t]*f[t]
  }
}
",fill = TRUE)
sink()


# Initial values
inits <- function(){list(mean.sjuv = runif(1, 0, 1), mean.sad = runif(1, 0, 1), mean.p = runif(1, 0, 1), mean.fec = runif(1, 0, 10), sigma.y = runif(1, 0, 1), n1 = rpois(1, 30), nad = rpois(1, 30))}

# Parameters monitored
parameters <- c("mean.sjuv", "mean.sad", "mean.p", "mean.fec", "N1", "Nad",
  "Ntot", "sigma2.y", "lambda")

# MCMC settings
ni <- 20000
nt <- 6
nb <- 5000
nc <- 3

library(jagsUI)
# Call JAGS from R (BRT 2 min)
ipm <- jags(data, inits, parameters, "ipm.jags", n.chains = nc,
  n.thin = nt, n.iter = ni, n.burnin = nb)

print(ipm, digits = 3)

# Note that sometimes the model does not run, because of the estimation of
# the population growth rate (lambda). Presumably, awkward numbers such as 0
# are produced during the initialising phase and then the ratio cannot be
# computed. There are two options: 1) just run the model again until it works,
# or 2) delete the calculation of lambda within the JAGS code and compute it
# afterwards based on the posterior samples of Ntot.


# Produce Fig. 11-4
par(cex = 1.2)
lower <- upper <- numeric()
for (i in 1:10){
   lower[i] <- quantile(ipm$sims.list$Ntot[,i], 0.025)
   upper[i] <- quantile(ipm$sims.list$Ntot[,i], 0.975)
   }
plot(ipm$mean$Ntot, type = "b", ylim = c(0, 100),
  ylab = "Population size", xlab = "Year", las = 1, pch = 16, col = "blue",
  frame = F, cex = 1.5)
segments(1:10, lower, 1:10, upper, col = "blue")
points(SUR, type = "b", col = "black", pch = 16, lty = 2, cex = 1.5)
legend(x = 1, y = 65, legend = c("Counts", "Estimates"), pch = c(16, 16),
  col = c("black", "blue"), lty = c(2, 1), bty = "n")

