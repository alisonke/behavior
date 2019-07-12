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
mu.alpha.d<-log(2.5)					# mean of species-level random effect on intercept of sigma
sig.alpha.d<- 0.25					# standard deviation of species-level random effect on intercept of sigma
beta.d<- -0.2					# fixed effect of observation covariate on sigma

alpha.d<-rnorm(n.behav, mu.alpha.d, sig.alpha.d)

#alpha.d <- log(4)      # behavior detection intercepts for this species
#beta.d<- -0.4				# fixed effect of observation covariate on detection

###look at distribution of sigma intercepts 
hist(exp(rnorm(1000, mu.alpha.d, sig.alpha.d)))

#smallest intercept and largest distance
g(3,3)

noise<-rbinom(nSites,1,.6)			# binary observation covariate

###makes a behavior by site matrix for Scale parameter of half-normal detection function 
sigma <- exp(matrix(alpha.d, nrow=n.behav, ncol=nSites) + 
               matrix(beta.d*noise, nrow=n.behav, ncol=nSites, byrow=T))  

############################################################################################################
###abundance component, with one covariate

landuse<-rnorm(nSites, 0, 1) ##spatial covariate


###abundance component, with one covariate (random species level random effect)
#Intercept
mu.alpha.a<-log(1.5)				# mean of species-level random effect on intercept of log(expected abundance)
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
  
  ###species specific parameters
  
  for (b in 1:n.behav){
    alpha.sigma[b]~dnorm(mu_sigma, tau_sigma)
    beta.lambda[b]~dnorm(mu_beta.lambda,tau_beta.lambda)
    alpha.lambda[b]~dnorm(mu_alpha.lambda,tau_alpha.lambda)
  }
  
  ###hyperparameters of species level random effects
  mu_sigma~dnorm(0,0.01)
  tau_sigma<-1/(sig_sigma*sig_sigma)
  sig_sigma~dunif(0,500)
  
  mu_alpha.lambda~dnorm(0,0.01)
  sig_alpha.lambda<-1/sqrt(tau_alpha.lambda)
  tau_alpha.lambda~dgamma(0.1,0.1)
  
  mu_beta.lambda~dnorm(0,0.01)
  sig_beta.lambda<-1/sqrt(tau_beta.lambda)
  tau_beta.lambda~dgamma(0.1,0.1)
  
  
  ### fixed observation coefficient
  
  beta.sigma~dnorm(0,0.01)
  
  
  for (s in 1:n.behav){
    
    for (j in 1:nSites){
      
      sigma[s,j]<-exp(alpha.sigma[s] + beta.sigma*noise[j])
      
      f.0[s,j] <- 2 * dnorm(0,0, 1/sigma[s,j]^2)
      
      for(k in 1:n.int){
        
        ### actual integral over distance categories
        
        up[s,j,k]<-pnorm(db[k+1], 0, 1/sigma[s,j]^2) 
        low[s,j,k]<-pnorm(db[k], 0, 1/sigma[s,j]^2) 
        p[s,j,k]<- 2 * ( up[s,j,k] - low[s,j,k])
        f[s,j,k]<- p[s,j,k]/f.0[s,j]                         
        fc[s,j,k]<- f[s,j,k] * pi[k]       
        fct[s,j,k]<-fc[s,j,k]/sum(fc[s,j,1:n.int])  
        
      }
      
      
      pcap[s,j]<-sum(fc[s,j,1:n.int])    # overall detection probability
      
      lambda[j,s]<- exp(alpha.lambda[s] + beta.lambda[s]*landuse[j])
      
      y[j,s]~ dbin(pcap[s,j],N[j,s])
      N[j,s]~dpois(lambda[j,s])
      
      
      ###create replicate abundances for Bayesian p-value on abundance component
      
      Nnew[j,s]~dpois(lambda[j,s])
      
      ### residuals for 'observed' and new abundances 
      FT1[j,s]<-pow(sqrt(N[j,s])-sqrt(lambda[j,s]),2)
      FT1new[j,s]<-pow(sqrt(Nnew[j,s])-sqrt(lambda[j,s]),2)
    }
    
    T1p[s]<-sum(FT1[1:nSites,s])
    T1newp[s]<-sum(FT1new[1:nSites,s])
  }
  
  # Bayesian p-value
  Bp.N<-sum(T1newp[1:n.behav])>sum(T1p[1:n.behav])
  
}
"


modeltext<-textConnection(modeltext)

### initial values for N
N.in<-t(y.sum)+1

inits<-function(){list(N=N.in, mu_alpha.lambda=runif(1,0,1), tau_alpha.lambda=runif(1,0,1), mu_beta.lambda=runif(1), tau_beta.lambda=runif(1), 
                        mu_sigma = runif(1, 0.5, 1.5), sig_sigma=runif(1) )}

### set parameters to monitor
params1<-c('mu_sigma', 'sig_sigma', 'mu_alpha.lambda', 'sig_alpha.lambda', 'mu_beta.lambda', 'sig_beta.lambda','Bp.N', 
            'beta.lambda','alpha.lambda', 'alpha.sigma', 'beta.sigma')

### compile and adapt JAGS model, then generate posterior samples (adjust n.iter for n=5 to 20000)
mod<-jags.model(modeltext, data, inits, n.chain=3, n.adapt=500)
out<-coda.samples(mod, params1, n.iter=8000, thin=8)

summary(out)
gelman.diag(out)
max(gelman.diag(out)$psrf[,1])

traceplot(out)

