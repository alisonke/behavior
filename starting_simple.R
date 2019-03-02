####Starting simple again

library("boot")
library("rjags")

#Generate data
set.seed(432104)
#######abundance process############################

#set number of sites, visits, behaviors
site<-10 #number of sites
visit<-3 #number of visits
behav<-4 #number of behaviors

#### simulate data ####
#simulate forest covers
forest<-c(rep(0,site)) #empty forest vector
for (j in 1:site){
  forest[j]<-runif(10,0,1.2)} #forest cover for each site

#make up a beta for relationship between forest and abundance
beta<-4

#simulate actual abundance data for single species
lambda<-c(rep(0,site)) #empty lambda vector
N<-c(rep(0,site)) #empty N vector
for (j in 1:site){
  lambda[j]<-exp(beta*forest[j]) #lambda for each site
  N[j]<-rpois(1,lambda[j]) }#N for each site

### behavior process
#simulate actual animals doing behaviors. each matrix is a different visit.
##create underlying behavior probabilities (intercepts)

behav.prob.int <- c(.2,.1,.6,.1)

#each behavior is controlled in some way by forest cover. for example, let's pretend it's eating, breeding, calling, aggression

behav.prob.beta <- c(4,2,3,1)

behav.subprob <- behav.prob.int + behav.prob.beta # + species level covariate?

behav.prob <- behav.subprob/sum(behav.subprob)

D<-array(c(rep(0,site*behav*visit)),dim=c(site,behav,visit)) #behavior matrix

for (j in 1:site){
  for (k in 1:visit){
    D[j,1:behav,k] <- rmultinom(1,N[j],behav.prob)}}

### detection process ###

#make up a beta for relationship between forest and detection

Bd = -2

#make up an alpha intercept for detection of each behavior that changes at each visit depending on Mu.v and Tau.v

Mu.v <- c(1,2,3,4)  # made up Mu.v values
Tau.v <- c(.5,.2,.4,.3) #made up Tau.v values

alpha0 <- matrix(c(rep(0,behav*visit)),nrow=visit) #empty alpha matrix

for(k in 1:visit){
  for(b in 1:behav){
    alpha0[k,b]<-rnorm(1,Mu.v[b],Tau.v[b])
  }
}

#simulate detection probabilities

detect <-array(c(rep(0,site*behav*visit)),dim=c(site,behav,visit)) #empty detection matrix

for(k in 1:visit){
  for(b in 1:behav){
    for(j in 1:site){
      detect[j,b,k]<-inv.logit(alpha0[k,b]+Bd*forest[j])
    }
  }
}

#simulate observations
Y <-array(c(rep(0,site*behav*visit)),dim=c(site,behav,visit)) #empty observation matrix

for(k in 1:visit){
  for(b in 1:behav){
    for(j in 1:site){
      Y[j,b,k]<-rbinom(1,D[j,b,k],detect[j,b,k])
    }
  }
}


####### now model abundance with detection parameter #############
### compile data for JAGS model

Ya<-matrix(c(rowSums(Y[,,1]),rowSums(Y[,,2]),rowSums(Y[,,2])),nrow=3,byrow=T)

data <- list(nSites = site, nVisits = visit, forest = forest, nBehav = behav, Ya=Ya, Y=Y)

model.string <-"
model {
# Priors
beta ~ dnorm(0,0.001)
Bd ~ dnorm(0,0.0001)
b.beta[1] <- 0; # zero contrast for baseline behavior
for (j in 2 : nBehav) { b.beta[j] ~ dnorm(0, 0.0001)} # vague priors for other behaviors

# Likelihood, detection plus a random error intercept

for(i in 1:nSites){
logit(detect[i]) = Bd*forest[i]
}

# Likelihood, abundance

for (i in 1:nSites){
for (k in 1:nVisits){
Ya[k,i] ~ dbin(detect[i],N[i])
}}

for (i in 1:nSites){
lambda[i]<-exp(beta*forest[i])
N[i] ~ dpois(lambda[i])}

for (i in 1:nSites){
for(k in 1:nVisits){
Y[i,1:4,k] ~ dmulti(p[i,1:4] , Ya[k,i])}}

for (i in 1:nSites){
for(j in 1:nBehav){
p[i,j] <- phi[i,j]/sum(phi[i,])
log(phi[i,j]) <- b.beta[j]*forest[i]
}
}

}"
modeljags<-textConnection(model.string)

#inits function
inits <- function(){list(N=c(rep(4*(max(Y)+1),10)), Nb=matrix(c(rep(max(Y)+1,40)),nrow=10),beta= rnorm(1), Bd = rnorm(1), Mu.v =rep(rnorm(1),3), Tau.v= rep(rlnorm(1),3))}
# Parameters to estimate
params <- c("beta","Bd", "Mu.v", "Tau.v","lambda","N","detect","alpha0","b.beta")

# MCMC settings
nc =3 ; ni= 1200 ; nb =200 ; nt= 1

# Start Gibbs sampler
sat.jags <- jags.model(modeljags,data=data,n.chains=3,n.adapt =1000,inits=inits)
samps.coda <- coda.samples(sat.jags, params, ni, nt,n.burnin=nb)
samps.jags <- jags.samples(sat.jags, params, ni, nt,n.burnin=nb)
print(samps.jags)
summary(samps.coda)
