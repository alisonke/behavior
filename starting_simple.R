####Starting simple again

###goals: to get an estimated abundance of behavior for each site
#         to get a coefficient for effect of forest cover on abundance of each behavior

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

behav.prob.beta <- c(2,4,3,1)

behav.subprob <- behav.prob.int + behav.prob.beta # + species level covariate?

behav.prob <- behav.subprob/sum(behav.subprob)

D<-array(c(rep(0,site*behav*visit)),dim=c(site,behav,visit)) #behavior matrix

for (j in 1:site){
  for (k in 1:visit){
    D[j,1:behav,k] <- rmultinom(1,N[j],behav.prob)}}

### detection process ###

#make up a beta for relationship between forest and detection

Bd <- -2

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
    }}}

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
Betadet ~ dnorm(0,.0001)
b.beta[1] <- 0; # zero contrast for baseline behavior
for (j in 2 : nBehav) { b.beta[j] ~ dnorm(0, 0.0001)} # vague priors for other behaviors. this is effect of forest on each behavior j

# Likelihood, abundance
##Ya is total abundance observed at each visit.
#for each site, total abundance is a poisson distribution. forest cover affects total abundance through beta.

for (i in 1:nSites){
lambda[i]<-exp(beta*forest[i])
N[i] ~ dpois(lambda[i])}

##out of the abundances, there are 4 behaviors that come from a multinomial distribution
#Nb is the actual abundance of each behavior 

for (i in 1:nSites){
Nb[i,1:4] ~ dmulti(p[i,1:4] , N[i])}

##the behavior probabilities relate to forest cover (b.beta). They are constant across visits.

for(i in 1:nSites){
for(j in 1:nBehav){
p[i,j] <- phi[i,j]/sum(phi[i,])
log(phi[i,j]) <- b.beta[j]*forest[i]
}
}

###Detection process. detection probability is constant across behaviors and visits for now.

#this is effect of forest cover on detection. Here, it is the same for each behavior (not true in simulation)

for(i in 1:nSites){
logit(detect[i]) = -2*forest[i]
}
for (i in 1:nSites){
for(j in 1:nBehav){
for (k in 1:nVisits){
Y[i,j,k] ~ dbin(detect[i],Nb[i,j])
}}}

##the values I want: Nb, b.beta


}"
modeljags<-textConnection(model.string)

#inits function
inits <- function(){list(Betadet=-2, N=c(rep(4*(max(Y)+1),10)), Nb=matrix(c(rep(max(Y)+1,40)),nrow=10),beta= rnorm(1), Mu.v =rep(rnorm(1),3), Tau.v= rep(rlnorm(1),3))}
# Parameters to estimate
params <- c("beta","Betadet", "lambda","N","Nb","detect","b.beta")

# MCMC settings
nc =3 ; ni= 10000 ; nb =1000 ; nt= 1

# Start Gibbs sampler
sat.jags <- jags.model(modeljags,data=data,n.chains=3,n.adapt =1000,inits=inits)
samps.coda <- coda.samples(sat.jags, params, ni, nt,n.burnin=nb)
samps.jags <- jags.samples(sat.jags, params, ni, nt,n.burnin=nb)
print(samps.jags)
summary(samps.coda)

##this runs at least. I had to change Betadet to a number and the N[i] abundance estimates are constant, which is a problem.
