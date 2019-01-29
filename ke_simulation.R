library(rjags)

set.seed(1111)

#number of species in the community
n.spec<-2

# Half-normal distance detection function
g <- function(x, sig) exp(-x^2/(2*sig^2))

#number of surveys
nSites <- 30
samp.rad <- 50 #radius of sampling circle
int.w<-10					# width of distance categories (v)
dist.breaks <- seq(0, samp.rad, by=int.w)	# distance break points
nG<-length(dist.breaks)-1	# number of distance categories


############################################################################################################

### detection component = sigma,  with one fixed effect
mu.s<-log(2.5)  # mean of species-level random effect on intercept of sigma
sig.s<- .9	  # standard deviation of species-level random effect on intercept of sigma 
beta.s<- -0.2					# fixed effect of observation covariate (e.g. forest cover) on sigma


###look at distribution of sigma intercepts 
hist(exp(rnorm(1000, mu.s, sig.s)))

##detection prob at farthest distance interval for largest sigma
g(50,30)


############################################################################################################
###abundance component, with one covariate (random species level random effect)
#Intercept
mu.lam.alpha<-log(1.5)				# mean of species-level random effect on intercept of log(expected abundance)
sig.lam.alpha<-1				# SD of species-level random effect on intercept of log(expected abundance)
mu.b1<-0					# mean of species-level random effect on coefficient of log(expected abundance)
sig.b1<-0.5					# SD of species-level random effect on coefficient of log(expected abundance)



###########################################################################################################
### begin iterations ######################################################################################
niter<-2					# number of iterations
iter<-1						# starting iteration

while(iter<=niter){
  
  print (iter)
  
  
  #####simulate site-specific binary covariate and species and site specific detection parameter sigma 
  
  obscov<-rbinom(nSites, 1, 0.6)			# observation covariate
  s.alpha<-rnorm(n.spec, mu.s, sig.s)		# detection intercept
  
  ###makes a species by site matrix for Scale parameter of half-normal detection function 
  sigma <- exp(matrix(s.alpha, nrow=n.spec, ncol=nSites) + 
                 matrix(beta.s*obscov, nrow=n.spec, ncol=nSites, byrow=T) )  
  
  
  ##### Simulate abundance across sites
  
  ##sp-specific intercept
  lam.alpha <- rnorm(n.spec, mu.lam.alpha, sig.lam.alpha)
  
  ##sp-specific coefficient for covariate
  b1<-rnorm(n.spec, mu.b1, sig.b1)
  
  ##spatial covariate
  Ncov<-rnorm(nSites, 0, 1)
  
  ##Poisson mean (log(expected abundance))
  lambda<-exp(matrix(lam.alpha, nrow=n.spec, ncol=nSites) + 
                matrix(rep(b1,each=nSites )*rep(Ncov, times=n.spec), nrow=n.spec, ncol=nSites, byrow=T) ) 
  
  ##abundance
  N <- matrix(rpois(n.spec*nSites, as.vector(lambda)), nrow=n.spec, ncol=nSites ) 
  
  ### total number of individuals in all sampled transects for each species
  N.tot<-apply(N,1,sum)
  
  
  #####simulate continuous distance data
  # y=number of individuals detected in each distance interval
  y <- array(0, c(n.spec, nSites, length(dist.breaks)-1))
  
  for (i in 1:n.spec){
    for(j in 1:nSites) {
      if(N[i,j] == 0)
        next
      
      # Distance from observer to the individual
      d <- runif(N[i,j], 0, strip.width) 		# uniform distribution of animals
      
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
  
  for (i in 1:n.spec){
    for(j in 1:nSites){
      for (k in 1:nG){
        if (y[i,j,k]==0) next
        spp<-c(spp, rep(i, y[i,j,k]))		#species index
        sst<-c(sst, rep(j, y[i,j,k]))		#site index
        dclass<-c(dclass, rep(k, y[i,j,k]))	#distance category index
        
      }}}
  
  
  ###write data to .R file for post-processing
  dat<-list(N=N, y=y, b1=b1, lam.alpha=lam.alpha, s.alpha=s.alpha, obscov=obscov, Ncov=Ncov)
  dput(dat, paste('Data_spec',n.spec,'_', iter, '.R', sep=''))
  
  
  ################## run JAGS model ##############################################################
  
  ### compile data for JAGS model
  data1<-list(spec=n.spec, nG=nG, db=dist.breaks, v=int.w, 
              pi=rep(1/(length(dist.breaks)-1), length(dist.breaks)-1), nsites=nSites, v1=Ncov,
              y=t(y.sum), nind=nind, dclass=dclass, species=spp, site=sst, OBSVAR=obscov)
  
  #xg=dist.breaks[-1]-0.5
  
  ### create initial values, the ones for N are important!
  N.in<-t(y.sum)+1
  
  inits1<-function(){list(N=N.in, mu_a=runif(1,0,1), tau_a=runif(1,0,1), mu_b1=runif(1), tau_b1=runif(1), 
                          mu_s = runif(1, 0.5, 1.5), sig_s=runif(1) )}
  
  ### set parameters to monitor
  params1<-c('mu_s', 'sig_s', 'mu_a', 'sig_a', 'mu_b1', 'sig_b1','Bp.N', 
             'Bp.Obs', 'beta1','alpha', 'Nspec', 'asig', 'bsig')
  
  ### read in JAGS model file
  ### NOTE: JAGS model code below!!
  
  modelFile1='Community_DS_Simulations.txt'
  
  ### compile and adapt JAGS model, then generate posterior samples (adjust n.iter for n=5 to 20000)
  mod<-jags.model(modelFile1, data1, inits1, n.chain=1, n.adapt=500)
  out<-coda.samples(mod, params1, n.iter=8000, thin=8)
  
  ### save model output for post-processing
  dput(out,  paste('Output_spec',n.spec,'_', iter, '.R', sep=''))
  
  
  iter<-iter+1
} ##end iteration loop





###################################################################################################
####################################################################################################


