##next step: add random effect for site due to multiple visits
library(rjags)

set.seed(2019)

#number of species in the community
n.spec=5

#number of behaviors. seen or heard for example
n.behav<-2 

# Half-normal detection function
g <- function(x, sig) exp(-x^2/(2*sig^2))

n.sites <- 50					# number of point count surveys
strip.width <- 5 				# width of each interval
int.w<-1					# width of distance categories (v)
dist.breaks <- seq(0, strip.width, by=int.w)	# distance break points
nG<-length(dist.breaks)-1			# number of distance categories


############################################################################################################
### detection component: sigma with one covariate. depends on behavior and not species. In the future, depend on strata too? Everything starting with s = SIGMA
# Each behavior has an underlying mean detection probability and fixed relationship wiht landuse
s.mu.alpha<-c(log(2.5), log(4))		#.92		# mean of behavior-level random effect on intercept of sigma 
s.sig.alpha<- c(0.25,.3)				# standard deviation of behavior-level random effect on intercept of sigma 
s.beta<- c(-0.2,.4)					# fixed effect of observation covariate on sigma. differs by behavior. For example, noise will affect heard more than seen observations

###look at distribution of sigma intercepts 
hist(exp(rnorm(1000, s.mu.alpha[2], s.sig.alpha[2])))

##detection prob at farthest distance interval for largest sigma
g(5,10)

#simulate site x behavior level detection parameter sigma

##simulate site-level binary covariate  
noise<-rbinom(n.sites, 1, 0.6)			# noise as example

#intercept depends on behavior
s.alpha<-c(rep(NA,n.behav))
for(i in 1:n.behav){
  s.alpha[i]<-rnorm(1,s.mu.alpha[i], s.sig.alpha[i])		# detection intercept
}

###makes a behavior by site matrix for scale parameter of half-normal detection function 
sigma<-matrix(NA,nrow=n.behav,ncol=n.sites)
for(i in 1:n.behav){
  for(j in 1:n.sites){
    sigma[i,j]<- exp(s.alpha[i]+s.beta[i]*noise[j])
  }
}

############################################################################################################
###abundance component, with one covariate, random effect species intercept and species x behavior level fixed effect. Everything starting with l = LAMBDA
#Intercept for random effect for species intercepts
l.mu.alpha<-log(2.5)	#.4	# mean of species-level random effect on intercept of log(expected abundance)
l.sig.alpha<-1				# SD of species-level random effect on intercept of log(expected abundance)

##sp-specific intercept
l.alpha <- rnorm(n.spec, l.mu.alpha, l.sig.alpha)

##species x behavior fixed coefficient for covariate
l.beta<-matrix(rnorm(n.spec*n.behav),nrow=n.spec)

##spatial covariate
landuse<-rnorm(n.sites, 0, 1)

##Poisson mean (log(expected abundance))
lambda<-array(NA,dim=c(n.spec,n.behav,n.sites))
for(i in 1:n.spec){
  for(j in 1:n.behav){
    for(k in 1:n.sites){
lambda[i,j,k]<-exp(l.alpha[i] + l.beta[i,j]*landuse[k])
}}}

##abundance
N<-array(NA,dim=c(n.spec,n.behav,n.sites))
for(i in 1:n.spec){
  for(j in 1:n.behav){
    for(k in 1:n.sites){
      N[i,j,k]<-rpois(1,lambda[i,j,k])
    }}}

### total number of individuals in all sampled transects for each species
N.tot<-apply(N,1,sum)

#####simulate continuous distance data
# y=number of individuals detected in each distance interval
y <- array(0, c(n.spec, n.behav, n.sites, length(dist.breaks)-1))

for(i in 1:n.spec){
  for(j in 1:n.behav){
    for(k in 1:n.sites){
    if(N[i,j,k] == 0)
      next
      #probabilties of occurring in each circular band
    probs<-((dist.breaks[-1]^2-(dist.breaks[-1]-1)^2))/((tail(dist.breaks,1))^2)
    # Distance from observer to the individual
    d <- rcat(N[i, j, k], probs)		# uniform distribution of animals
    p <- g(x=d, sig=sigma[j,k])   		# Detection probability
    seen <- rbinom(N[i,j,k], 1, p) 		# Which individuals are detected
    if(all(seen == 0))
      next
    d1 <- d[seen==1] 				# The distance data for seen individuals
    counts <- table(cut(d1, dist.breaks, include.lowest=TRUE))
    y[i,j,k,] <- counts 				# The number of detections in each distance interval
  }
}
}

y.sum<-apply(y,1:3, sum) #combining all the distance bands together

##skip data sets with unobserved species or species with abundance=0
if (any(apply(y.sum,1,sum)==0) | any(N.tot==0) ) next

##### if data passes both criteria, continue on ################################################

##### convert data to JAGS format

nind<-sum(y)

spp<-sst<-dclass<-NULL

for (i in 1:n.spec){
  for(j in 1:n.behav){
    for(k in 1:n.sites){
      for (l in 1:nG){
        if (y[i,j,k,l]==0) next
        sp.ind<-c(spp, rep(i, y[i,j,k,l]))		#species index
        b.ind<-c(beh, rep(j, y[i,j,k,l]))		#behavior index
        st.ind<-c(sst, rep(k, y[i,j,k,l]))		#site index
        d.ind<-c(dclass, rep(l, y[i,j,k,l]))	#distance category index
      
    }}}}


################## run JAGS model #################################################

### compile data for JAGS model
data <- list(
    n.spec = n.spec,
    n.behav = n.behav,
    nG = nG,
    db = dist.breaks,
    int.w = int.w,
    pi = probs,
    n.sites = n.sites,
    landuse = landuse,
    y = y.sum,
    #sp.ind = sp.ind,
    #b.ind = b.ind,
    #st.ind = st.ind,
    #d.ind = d.ind,
    noise = noise
  )

modeltext <- "model{

#behavior specific parameters
#Sigma depends on behavior
for(j in 1:n.behav){
s.alpha[j]~dnorm(s.mu.alpha[j], s.tau.alpha[j])
s.beta[j]~dnorm(0,0.01)
}

#species specific parameters
for (i in 1:n.spec){
l.alpha[i]~dnorm(l.mu.alpha,l.tau.alpha)
}

###species and behavior specific parameters

for (i in 1:n.spec){
for(j in 1:n.behav){
l.beta[i,j]~dnorm(0,.01)
}}

###hyperparameters of species or behavior level random effects
for(j in 1:n.behav){
s.mu.alpha[j]~dnorm(0,0.01)
s.tau.alpha[j]<-1/(s.sig.alpha[j]*s.sig.alpha[j])
s.sig.alpha[j]~dunif(0,500)}

l.mu.alpha~dnorm(0,0.01)
l.sig.alpha<-1/sqrt(l.tau.alpha)
l.tau.alpha~dgamma(0.1,0.1)

for(j in 1:n.behav){
for (k in 1:n.sites){
sigma[j,k]<-exp(s.alpha[j] + s.beta[j]*noise[k])
f.0[j,k] <- 2 * dnorm(0,0, 1/sigma[j,k]^2)

for(l in 1:nG){

### actual integral over distance categories
up[j,k,l]<-pnorm(db[l+1], 0, 1/sigma[j,k]^2) 
low[j,k,l]<-pnorm(db[l], 0, 1/sigma[j,k]^2) 
p[j,k,l]<- 2 * ( up[j,k,l] - low[j,k,l])
f[j,k,l]<- p[j,k,l]/f.0[j,k]/int.w                          
fc[j,k,l]<- f[j,k,l] * pi[l]       
fct[j,k,l]<-fc[j,k,l]/sum(fc[j,k,1:nG])  
}

pcap[j,k]<-sum(fc[j,k,1:nG])    # overall detection probability

for(i in 1:n.spec){
lambda[i,j,k]<- exp(l.alpha[i] + l.beta[i,j]*landuse[k])

y[i,j,k]~ dbin(pcap[j,k],N[i,j,k])
N[i,j,k]~dpois(lambda[i,j,k])

###create replicate abundances for Bayesian p-value on abundance component

N.new[i,j,k]~dpois(lambda[i,j,k])

 ### residuals for 'observed' and new abundances 
resid[i,j,k]<-pow(sqrt(N[i,j,k])-sqrt(lambda[i,j,k]),2)
resid.new[i,j,k]<-pow(sqrt(N.new[i,j,k])-sqrt(lambda[i,j,k]),2)

}}}

for(i in 1:n.spec){for(j in 1:n.behav){
sum.resid[i,j]<-sum(resid[i,j,1:n.sites])
sum.resid.new[i,j]<-sum(resid.new[i,j,1:n.sites])
}}

# Bayesian p-value
BP.N<-sum(sum.resid.new)>sum(resid.new)

}
"

modelFile1<-textConnection(modeltext)


### create initial values, the ones for N are important!
N.in<-y.sum+1

inits<-function(){list(N=N.in, l.mu.alpha=runif(1,0,1), l.tau.alpha=runif(1,0,1), 
                        s.mu.alpha = runif(n.behav, 0.5, 1.5), s.sig.alpha=runif(n.behav),s.beta=runif(2) )}

### set parameters to monitor
params1<-c('s.mu.alpha', 's.sig.alpha', 'l.mu.alpha', 'l.sig.alpha', 'l.beta','l.alpha', 'Nspec', 's.alpha', 's.beta','BP.N')


mod<-jags.model(modelFile1, data, inits, n.chain=3, n.adapt=500)
out<-coda.samples(mod, params1, n.iter=8000, thin=8)

summary(out)
gelman.diag(out)



