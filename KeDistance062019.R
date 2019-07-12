
library(rjags)

set.seed(2019)

#number of species in the community
n.behav<-2 # 5, 15, 30

# Half-normal detection function
g <- function(x, sig) exp(-x^2/(2*sig^2))

nSites <- 30					# number of line transect surveys
strip.width <- 5 				# width of each interval
int.w<-1					# width of distance categories (v)
dist.breaks <- seq(0, strip.width, by=int.w)	# distance break points
nG<-length(dist.breaks)-1			# number of distance categories


############################################################################################################
### detection component: sigma with one covariate (fixed effect)
d.mu.alpha<-log(2.5)		#.92			# mean of species-level random effect on intercept of sigma
d.sig.alpha<- 0.25					# standard deviation of species-level random effect on intercept of sigma
d.beta<- -0.2					# fixed effect of observation covariate on sigma

###look at distribution of sigma intercepts 
hist(exp(rnorm(1000, d.mu.alpha, d.sig.alpha)))

##detection prob at farthest distance interval for largest sigma
g(10,6)


############################################################################################################
###abundance component, with one covariate (random species level random effect)
#Intercept
l.mu.alpha<-log(2.5)	#.4			# mean of species-level random effect on intercept of log(expected abundance)
l.sig.alpha<-1				# SD of species-level random effect on intercept of log(expected abundance)
l.mu.beta<-0					# mean of species-level random effect on coefficient of log(expected abundance)
l.sig.beta<-0.5					# SD of species-level random effect on coefficient of log(expected abundance)





  
  #####simulate site-specific binary covariate and species and site specific detection parameter sigma 
  
  noise<-rbinom(nSites, 1, 0.6)			# observation covariate
  s.alpha<-rnorm(n.behav, d.mu.alpha, d.sig.alpha)		# detection intercept
  
  ###makes a species by site matrix for Scale parameter of half-normal detection function 
  sigma <- exp(matrix(s.alpha, nrow=n.behav, ncol=nSites) + 
                 matrix(d.beta*noise, nrow=n.behav, ncol=nSites, byrow=T) )  
  
  
  ##### Simulate abundance across sites
  
  ##sp-specific intercept
  l.alpha <- rnorm(n.behav, l.mu.alpha, l.sig.alpha)
  
  ##sp-specific coefficient for covariate
  l.beta<-rnorm(n.behav, l.mu.beta, l.sig.beta)
  
  ##spatial covariate
  landuse<-rnorm(nSites, 0, 1)
  
  ##Poisson mean (log(expected abundance))
  lambda<-exp(matrix(l.alpha, nrow=n.behav, ncol=nSites) + 
                matrix(rep(l.beta,each=nSites )*rep(landuse, times=n.behav), nrow=n.behav, ncol=nSites, byrow=T) ) 
  
  ##abundance
  N <- matrix(rpois(n.behav*nSites, as.vector(lambda)), nrow=n.behav, ncol=nSites ) 
  
  ### total number of individuals in all sampled transects for each species
  N.tot<-apply(N,1,sum)
  
  
  #####simulate continuous distance data
  # y=number of individuals detected in each distance interval
  y <- array(0, c(n.behav, nSites, length(dist.breaks)-1))
  
  for (i in 1:n.behav){
    for(j in 1:nSites) {
      if(N[i,j] == 0)
        next
      probs<-((dist.breaks[-1]^2-(dist.breaks[-1]-1)^2))/((tail(dist.breaks,1))^2)
      # Distance from observer to the individual
      d <- rcat(N[i, j], probs)		# uniform distribution of animals
      
      p <- g(x=d, sig=sigma[i,j])   		# Detection probability
      seen <- rbinom(N[i,j], 1, p) 		# Which individuals are detected
      if(all(seen == 0))
        next
      d1 <- d[seen==1] 				# The distance data for seen individuals
      counts <- table(cut(d1, dist.breaks, include.lowest=TRUE))
      y[i,j,] <- counts 				# The number of detections in each distance interval
    }
  }
  
  
  y.sum<-apply(y,1:2, sum)
  
  ##skip data sets with unobserved species or species with abundance=0
  if (any(apply(y.sum,1,sum)==0) | any(N.tot==0) ) next
  
  
  ##### if data passes both criteria, continue on ################################################
  
  ##### convert data to JAGS format
  
  nind<-sum(y)
  
  spp<-sst<-dclass<-NULL
  
  for (i in 1:n.behav){
    for(j in 1:nSites){
      for (k in 1:nG){
        if (y[i,j,k]==0) next
        spp<-c(spp, rep(i, y[i,j,k]))		#species index
        sst<-c(sst, rep(j, y[i,j,k]))		#site index
        dclass<-c(dclass, rep(k, y[i,j,k]))	#distance category index
        
      }}}
  
  
  ###write data to .R file for post-processing
  dat<-list(N=N, y=y, l.beta=l.beta, l.alpha=l.alpha, s.alpha=s.alpha, noise=noise, landuse=landuse)
  dput(dat, paste('Data_spec',n.behav,'_', iter, '.R', sep=''))
  
  
  ################## run JAGS model ##############################################################
  
  ### compile data for JAGS model
  data1<-list(spec=n.behav, nG=nG, db=dist.breaks, v=int.w, 
              pi=probs, nsites=nSites, v1=landuse,
              y=t(y.sum), nind=nind, dclass=dclass, species=spp, site=sst, OBSVAR=noise)
  
  #xg=dist.breaks[-1]-0.5
  
  ### create initial values, the ones for N are important!
  N.in<-t(y.sum)+1
  
  inits1<-function(){list(N=N.in, mu_a=runif(1,0,1), tau_a=runif(1,0,1), mu_b1=runif(1), tau_b1=runif(1), 
                          mu_s = runif(1, 0.5, 1.5), sig_s=runif(1) )}
  
  ### set parameters to monitor
  params1<-c('mu_s', 'sig_s', 'mu_a', 'sig_a', 'mu_b1', 'sig_b1','Bp.N', 
             'Bp.Obs', 'l.beta','alpha', 'Nspec', 'd.alpha', 'd.beta')
  

modeltext <- "model{
  
  ###species specific parameters
  
  for (s in 1:spec){
    d.alpha[s]~dnorm(mu_s, tau_s)
    l.beta[s]~dnorm(mu_b1,tau_b1)
    l.alpha[s]~dnorm(mu_a,tau_a)
  }
  
  ###hyperparameters of species level random effects
  mu_s~dnorm(0,0.01)
  tau_s<-1/(sig_s*sig_s)
  sig_s~dunif(0,500)
  
  mu_a~dnorm(0,0.01)
  sig_a<-1/sqrt(tau_a)
  tau_a~dgamma(0.1,0.1)
  
  mu_b1~dnorm(0,0.01)
  sig_b1<-1/sqrt(tau_b1)
  tau_b1~dgamma(0.1,0.1)
  
  
  ### fixed observation coefficient
  
  d.beta~dnorm(0,0.01)
  
  
  for (s in 1:spec){
    
    for (j in 1:nsites){
      
      sigma[s,j]<-exp(d.alpha[s] + d.beta*OBSVAR[j])
      
      f.0[s,j] <- 2 * dnorm(0,0, 1/sigma[s,j]^2)
      
      for(k in 1:nG){
        
        ### actual integral over distance categories
        
        up[s,j,k]<-pnorm(db[k+1], 0, 1/sigma[s,j]^2) 
        low[s,j,k]<-pnorm(db[k], 0, 1/sigma[s,j]^2) 
        p[s,j,k]<- 2 * ( up[s,j,k] - low[s,j,k])
        f[s,j,k]<- p[s,j,k]/f.0[s,j]/v                          
        fc[s,j,k]<- f[s,j,k] * pi[k]       
        fct[s,j,k]<-fc[s,j,k]/sum(fc[s,j,1:nG])  
        
      }
      
      
      pcap[s,j]<-sum(fc[s,j,1:nG])    # overall detection probability
      
      lambda[j,s]<- exp(l.alpha[s] + l.beta[s]*v1[j])
      
      y[j,s]~ dbin(pcap[s,j],N[j,s])
      N[j,s]~dpois(lambda[j,s])
      
      
      ###create replicate abundances for Bayesian p-value on abundance component
      
      Nnew[j,s]~dpois(lambda[j,s])
      
      ### residuals for 'observed' and new abundances 
      FT1[j,s]<-pow(sqrt(N[j,s])-sqrt(lambda[j,s]),2)
      FT1new[j,s]<-pow(sqrt(Nnew[j,s])-sqrt(lambda[j,s]),2)
    }
    
    T1p[s]<-sum(FT1[1:nsites,s])
    T1newp[s]<-sum(FT1new[1:nsites,s])
  }
  
  # Bayesian p-value
  Bp.N<-sum(T1newp[1:spec])>sum(T1p[1:spec])
  
  
  for(i in 1:nind){
    dclass[i] ~ dcat(fct[species[i],site[i],1:nG]) 
    
    ###generate new observations, calculate residuals for Bayesian p-value on detection component
    dclassnew[i] ~ dcat(fct[species[i],site[i],1:nG]) 
    Tobsp[i]<- pow(1- sqrt(fct[species[i],site[i],dclass[i]]),2)
    Tobspnew[i]<- pow(1- sqrt(fct[species[i],site[i],dclassnew[i]]),2)
  }
  
  Bp.Obs<-sum(Tobspnew[1:nind])>sum(Tobsp[1:nind])
  
  
  ###monitor total abundance
  for (i in 1:spec){
    Nspec[i]<-sum(N[1:nsites,i])
  }
  
  
}"

modelFile1<-textConnection(modeltext)

mod<-jags.model(modelFile1, data1, inits1, n.chain=3, n.adapt=500)
out<-coda.samples(mod, params1, n.iter=8000, thin=8)

summary(out)
gelman.diag(out)

