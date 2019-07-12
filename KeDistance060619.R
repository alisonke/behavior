library(rjags)
library(extraDistr)
set.seed(2019)

n.behav=10
# Half-normal detection function
g <- function(x, sig) exp(-x^2/(2*sig^2))

nSites <- 30	# number of point count surveys
strip.width = 5 #largest distance category
int.w <- 1	# width of distance categories (v)
dist.breaks <- seq(0, strip.width, by = int.w)	# distance break points
n.int <- length(dist.breaks) - 1			# number of distance categories

############################################################################################################
### detection component: sigma with intercept for each behavior and one covariate (fixed effect)

### detection component: sigma with one covariate (fixed effect)
mu.alpha.d<-log(3)					# mean of species-level random effect on intercept of sigma
sig.alpha.d<- 0.25					# standard deviation of species-level random effect on intercept of sigma
beta.d<- -0.2					# fixed effect of observation covariate on sigma

alpha.d<-rnorm(n.behav, mu.alpha.d, sig.alpha.d)

#alpha.d <- log(4)      # behavior detection intercepts for this species
#beta.d<- -0.4				# fixed effect of observation covariate on detection

###look at distribution of sigma intercepts 
hist(exp(rnorm(1000, mu.s, sig.s)))

#smallest intercept and largest distance
g(10,6)

noise<-runif(nSites)			# binary observation covariate

###makes a behavior by site matrix for Scale parameter of half-normal detection function 
sigma <- exp(matrix(alpha.d, nrow=n.behav, ncol=nSites) + 
               matrix(c(beta.d*noise), nrow=n.behav, ncol=nSites, byrow=T))  


############################################################################################################
###abundance component, with one covariate

landuse<-rnorm(nSites, 0, 1) ##spatial covariate


###abundance component, with one covariate (random species level random effect)
#Intercept
mu.alpha.a<-log(5)				# mean of species-level random effect on intercept of log(expected abundance)
sig.alpha.a<-1				# SD of species-level random effect on intercept of log(expected abundance)
mu.beta.a<-0					# mean of species-level random effect on coefficient of log(expected abundance)
sig.beta.a<-0.5					# SD of species-level random effect on coefficient of log(expected abundance)

#alpha.a<-5			# intercept of abundance for behaviors
#beta.a<-c(0 ,-0.3,  0.7,  1.2, -0.4)					# coefficient of land use on each behavior

alpha.a <- rnorm(n.behav, mu.alpha.a, sig.alpha.a)
beta.a<-rnorm(n.behav, mu.beta.a, sig.beta.a)

##Poisson mean (log(expected abundance))
lambda<-exp(matrix(alpha.a, nrow=n.behav, ncol=nSites) + 
              matrix(rep(beta.a,each=nSites)*rep(landuse, times=n.behav), nrow=n.behav, ncol=nSites, byrow=T) ) 

##abundance
N <- matrix(rpois(n.behav*nSites, as.vector(lambda)), nrow=n.behav, ncol=nSites ) 

### total number of individuals in all sampled points for each species
N.tot<-apply(N,1,sum)


#####simulate continuous distance data
# y=number of individuals detected in each distance interval
y <- array(0, c(n.behav, nSites, length(dist.breaks)-1))

for (i in 1:n.behav) {
  for (j in 1:nSites) {
    if (N[i, j] == 0)
      next
    
    # Distance from observer to the individual
    probs<-((dist.breaks[-1]^2-(dist.breaks[-1]-1)^2))/((tail(dist.breaks,1))^2)
    d <-rcat(N[i, j], probs)		# distribution of animals, probabilities proportional to areas of distance categories
    
    p <- g(x = d, sig = sigma[i, j]) # Detection probability of each distance category
    seen <- rbinom(N[i, j], 1, p)  # Which individuals are detected
    if (all(seen == 0))
      next
    d1 <- d[seen == 1] 				# The distance data for seen individuals
    counts <- table(cut(d1, dist.breaks, include.lowest = TRUE))
    y[i, j,] <-
      counts 	# The number of detections in each distance interval
  }
}

y.sum<-apply(y,1:2, sum)

##skip data sets with unobserved species or species with abundance=0
if (any(apply(y.sum,1,sum)==0) | any(N.tot==0) ) next

##### if data passes both criteria, continue on ################################################

##### convert data to JAGS format

nind<-sum(y)

beh<-site<-dclass<-NULL

for (i in 1:n.behav){
  for(j in 1:nSites){
    for (k in 1:n.int){
      if (y[i,j,k]==0) next
      beh<-c(beh, rep(i, y[i,j,k]))		#behavior index
      site<-c(site, rep(j, y[i,j,k]))		#site index
      dclass<-c(dclass, rep(k, y[i,j,k]))	#distance category index
      
    }}}


###write data to .R file for post-processing
dat<-list(N=N, y=y, beta.a=beta.a, alpha.a=alpha.a, noise=noise, landuse=landuse)


################## run JAGS model ##############################################################

### compile data for JAGS model
data <- list(
  n.behav = n.behav,
  n.int = n.int,
  db = dist.breaks,
  int.w = int.w,
  pi = probs,
  nSites = nSites,
  noise = noise,
  y = t(y.sum),
  nind = nind,
  dclass = dclass,
  behavior = beh,
  site = site,
  landuse = landuse
)



modeltext="model{

###abundance and detection parameters priors

for (b in 1:n.behav){
alpha.d[b]~dnorm(mu.alpha.d, tau.alpha.d)
beta.a[b]~dnorm(mu.beta.a,tau.beta.a)
alpha.a[b]~dnorm(mu.alpha.a,tau.alpha.a)
}
mu.alpha.d~dnorm(0,0.01)
tau.alpha.d<-1/(tau_s*tau_s)
tau_s~dunif(0,500)

mu.beta.a~dnorm(0,0.01)
tau.beta.a<-1/sqrt(tau_b1)
tau_b1~dgamma(0.1,0.1)

mu.alpha.a~dnorm(0,0.01)
tau.alpha.a<-1/sqrt(tau_a)
tau_a~dgamma(0.1,0.1)


### fixed observation coefficient

beta.d~dnorm(0,0.01)


### detection process

for (b in 1:n.behav){
for (j in 1:nSites){
sigma[b,j]<-exp(alpha.d[b] + beta.d*noise[j]) #sigma for detection probability
f.0[b,j] <- 2 * dnorm(0,0, 1/sigma[b,j]^2) #2 times density of normal distribution at zero

for(k in 1:n.int){
### integral over distance intervals

up[b,j,k]<-pnorm(db[k+1], 0, 1/sigma[b,j]^2) #area of upper interval
low[b,j,k]<-pnorm(db[k], 0, 1/sigma[b,j]^2) #area of lower interval
p[b,j,k]<- 2 * (up[b,j,k] - low[b,j,k]) #area of interval

f[b,j,k]<- p[b,j,k]/f.0[b,j]/int.w  #area of interval divided by full area (scaled) divided by interval width (1 in this case). probability of an individual being in this interval

fc[b,j,k]<- f[b,j,k] * pi[k]  #detection prob is prob of individual being detected there times prob of individual being there

fct[b,j,k]<-fc[b,j,k]/sum(fc[b,j,1:n.int])  #scale all of these to sum to 1
}

pcap[b,j]<-sum(fc[b,j,1:n.int])    # overall detection probability (note it is fc not fct)

###abundance process
lambda[j,b]<- exp(alpha.a[b] + beta.a[b]*landuse[j]) 

### for a flexible number of covariates, use:
### lambda[j,b]<- exp(alpha[b] + inprod(beta1[b,]*v1[j,]))
### see seabird application for an example

y[j,b]~ dbin(pcap[b,j],N[j,b]) #overall detection probability with N
N[j,b]~dpois(lambda[j,b])
}}
}"

modeltext<-textConnection(modeltext)

### initial values for N
N.in<-t(y.sum)+1

inits<-function(){list(N=N.in, mu.alpha.d=runif(1,0,1), tau_s=runif(1,0,1), mu.beta.a=runif(1), tau_b1=runif(1), mu.alpha.a = runif(1, 0.5, 1.5), tau_a=runif(1) )}

### set parameters to monitor
params1<-c('mu.alpha.d', 'mu.alpha.a', 'mu.beta.a','beta.d','tau_s','tau_b1','tau_a')

### compile and adapt JAGS model, then generate posterior samples (adjust n.iter for n=5 to 20000)
mod<-jags.model(modeltext, data, inits, n.chain=3, n.adapt=500)
out<-coda.samples(mod, params1, n.iter=8000, thin=8)

summary(out)
gelman.diag(out)
max(gelman.diag(out)$psrf[,1])

traceplot(out)
