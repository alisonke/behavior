###goal: to get an estimated coefficient for effect of forest cover on abundance of each behavior for a single species

library("boot")
library("rjags")

#Generate data
set.seed(4321)
#######abundance process############################

#set number of sites, visits, behaviors
site<-10 #number of sites
visit<-3 #number of visits
behav<-4 #number of behaviors

#### simulate data ####
#simulate forest covers
forest<-runif(site,0,1.2) #forest cover for each site

#make up a beta for relationship between forest and abundance
beta<-4

#simulate actual abundance data for a single species
lambda<-c(rep(0,site)) #empty lambda vector
N<-c(rep(0,site)) #empty N vector
for (i in 1:site){
  lambda[i]<-exp(beta*forest[i]) #lambda for each site
  N[i]<-rpois(1,lambda[i]) }#N for each site

### behavior process
#simulate animal doing behaviors. each matrix is a different visit.
##create underlying behavior probabilities (intercepts)

behav.prob.int <- c(.2,.1,.6,.1)

#betas: each behavior is controlled in some way by forest cover.

behav.prob.beta <- c(2,4,3,1)

D<-array(c(rep(0,site*behav*visit)),dim=c(site,behav,visit)) #behavior matrix

for (i in 1:site){
    for (k in 1:visit){
      behav.subprob<-behav.prob.int + behav.prob.beta*forest[i] #+species covariate?
      behav.prob <- behav.subprob/sum(behav.subprob)
      D[i,1:behav,k] <- rmultinom(1,N[i],behav.prob)}}

### detection process ###

#make up a beta for relationship between forest and detection

Bd <- -2

#make up an alpha intercept for detection of each behavior that changes at each visit depending on Mu.v and Tau.v

Mu.v <- c(1,2,3,4)  # made up Mu.v values
Tau.v <- c(.5,.2,.4,.3) #made up Tau.v values

alpha0 <- matrix(c(rep(0,behav*visit)),nrow=visit) #empty alpha matrix

for(k in 1:visit){
  for(j in 1:behav){
    alpha0[k,j]<-rnorm(1,Mu.v[j],Tau.v[j])
  }
}

#simulate detection probabilities

detect <-array(c(rep(0,site*behav*visit)),dim=c(site,behav,visit)) #empty detection matrix

for(k in 1:visit){
  for(j in 1:behav){
    for(i in 1:site){
      detect[i,j,k]<-inv.logit(alpha0[k,j]+Bd*forest[i])
    }}}

#simulate observations
Y <-array(c(rep(0,site*behav*visit)),dim=c(site,behav,visit)) #empty observation matrix

for(k in 1:visit){
  for(j in 1:behav){
    for(i in 1:site){
      Y[i,j,k]<-rbinom(1,D[i,j,k],detect[i,j,k])
    }
  }
}
####### now model abundance with detection parameter #############
### compile data for JAGS model
Ya<-matrix(c(rowSums(Y[,,1]),rowSums(Y[,,2]),rowSums(Y[,,2])),nrow=3,byrow=T)

data <- list(nSites = site, nVisits = visit, forest = forest, nBehav = behav, Y=Y)

model.string <-"
model {
#### Priors 

#forest effect on abundance
beta ~ dnorm(0,0.001)

#forest effect on detection
Betadet ~ dnorm(0,0.0001)

#forest effect on behavior
b.beta[1] <- 0; # zero contrast for baseline behavior
for (j in 2 : nBehav) { b.beta[j] ~ dnorm(0, 0.0001)} # vague priors for other behaviors

#### Likelihood, abundance

for (i in 1:nSites){
lambda[i]<-exp(beta*forest[i])
N[i] ~ dpois(lambda[i])}

##out of the abundances, there are 4 behaviors that come from a multinomial distribution
#Nb is the actual abundance of each behavior for each visit

for (i in 1:nSites){
for(k in 1:nVisits){
Nb[i,1:4,k] ~ dmulti(p[i,1:4] , N[i])}}

##the behavior probabilities relate to forest cover (b.beta). They are constant across visits.

for(i in 1:nSites){
for(j in 1:nBehav){
p[i,j] <- phi[i,j]/sum(phi[i,])
log(phi[i,j]) <- b.beta[j]*forest[i]
}
}

###Detection process. detection probability is constant across behaviors and visits for now (not true in simulation)

for(i in 1:nSites){
logit(detect[i]) = Betadet*forest[i]
}

for (i in 1:nSites){
for(j in 1:nBehav){
for (k in 1:nVisits){
Y[i,j,k] ~ dbin(detect[i],Nb[i,j,k])
}}}

}"
modeljags<-textConnection(model.string)

#inits function
inits <- function(){list(Betadet=rnorm(1),N=rep(max(Y)*4,10),beta= rnorm(1))}
# Parameters to estimate
params <- c("beta","Betadet","b.beta","lambda","N","Nb","detect")

# MCMC settings
nc =3 ; ni= 5000 ; nb =500 ; nt= 1

# Start Gibbs sampler
sat.jags <- jags.model(modeljags,data=data,n.chains=3,n.adapt =1000,inits=inits)
samps.coda <- coda.samples(sat.jags, params, ni, nt,n.burnin=nb)
samps.jags <- jags.samples(sat.jags, params, ni, nt,n.burnin=nb)
print(samps.jags)
summary(samps.coda)

#have to give N initial values for it to run, but then N estimates are constant
