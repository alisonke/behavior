
##next step: add random effect for site due to multiple visits
library(rjags)

set.seed(200)

#number of species in the community
n.spec=5

#number of behaviors. seen or heard for example
n.behav<-2 

# Half-normal detection function
g <- function(x, sig) exp(-x^2/(2*sig^2))

n.sites <- 30					# number of point count surveys
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



############################################################################################################
###abundance component, with one covariate, random effect species intercept and species x behavior level fixed effect. Everything starting with l = LAMBDA
#Intercept for random effect for species intercepts
l.mu.alpha<-log(4)	#.4	# mean of species-level random effect on intercept of log(expected abundance)
l.sig.alpha<-1				# SD of species-level random effect on intercept of log(expected abundance)


###########################################################################################################
### begin iterations ######################################################################################
niter<-1					# number of iterations
iter<-1						# starting iteration

while(iter<=niter){
  print (iter)

##sp-specific intercept
l.alpha <- rnorm(n.spec, l.mu.alpha, l.sig.alpha)

##species x behavior fixed coefficient for covariate
l.beta<-matrix(rnorm(n.spec*n.behav),nrow=n.spec)

##spatial covariates
landuse<-rnorm(n.sites, 0, 1)
noise<-rbinom(n.sites, 1, 0.6)

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

##sp-specific intercept
s.alpha <- rnorm(n.behav, s.mu.alpha, s.sig.alpha)


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
    sigma <- array(rep(NA,n.behav*n.sites),dim=c(n.behav,n.sites))
    sigma[j,k] <- exp(s.alpha[j]+s.beta[j]*noise[k])  
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

spec.site<-apply(y,c(1,3),sum)

##skip data sets with unobserved species or species with abundance=0
if (any(apply(y.sum,1,sum)==0) | any(N.tot==0) ) next

##### if data passes both criteria, continue on ################################################

##### convert data to JAGS format

n.ind<-sum(y)

sp.ind<-b.ind<-st.ind<-d.ind<-NULL

for (i in 1:n.spec){
  for(j in 1:n.behav){
    for(k in 1:n.sites){
      for (l in 1:nG){
        if (y[i,j,k,l]==0) next
        sp.ind<-c(sp.ind, rep(i, y[i,j,k,l]))		#species index
        b.ind<-c(b.ind, rep(j, y[i,j,k,l]))		#behavior index
        st.ind<-c(st.ind, rep(k, y[i,j,k,l]))		#site index
        d.ind<-c(d.ind, rep(l, y[i,j,k,l]))	#distance category index
      
    }}}}

###write data to .R file for post-processing
dat<-list(N=N, y=y, landuse=landuse, noise=noise, l.alpha=l.alpha, s.alpha=s.alpha, l.beta=l.beta)
setwd("~/Google Drive/Git/behavior/output")
dput(dat, paste('Data_sim',n.spec,'_', iter, '.R', sep=''))

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
    n.ind = n.ind,
    sp.ind = sp.ind,
    b.ind = b.ind,
    st.ind = st.ind,
    d.ind = d.ind,
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

for(i in 1:n.ind){
    d.ind[i] ~ dcat(fct[b.ind[i],st.ind[i],1:nG]) 
    
    ###generate new observations, calculate residuals for Bayesian p-value on detection component
    d.ind.new[i] ~ dcat(fct[b.ind[i],st.ind[i],1:nG]) 
    Tobsp[i]<- pow(1- sqrt(fct[b.ind[i],st.ind[i],d.ind[i]]),2)
    Tobspnew[i]<- pow(1- sqrt(fct[b.ind[i],st.ind[i],d.ind.new[i]]),2)
  }
  BP.d<-sum(Tobspnew[1:n.ind])>sum(Tobsp[1:n.ind])

  for (i in 1:n.spec){
    N.sp[i]<-sum(N[i,,]) #number of each species across all behaviors, sites, distances
for(k in 1:n.sites){
sp.abun[i,k]<-sum(N[i,,k]) #number of each species per site across behaviors
}
} 
}
"

modelFile<-textConnection(modeltext)


### create initial values, the ones for N are important!
N.in<-y.sum+1

inits<-function(){list(N=N.in, l.mu.alpha=runif(1,0,1), l.tau.alpha=runif(1,0,1), 
                        s.mu.alpha = runif(n.behav, 0.5, 1.5), s.sig.alpha=runif(n.behav),s.beta=runif(2) )}

### set parameters to monitor
params<-c('s.mu.alpha', 's.sig.alpha', 'l.mu.alpha', 'l.sig.alpha', 'l.beta','l.alpha', 'N.sp', 's.alpha', 's.beta','BP.N','BP.d','sp.abun')


mod<-jags.model(modelFile, data, inits, n.chain=3, n.adapt=500)
out<-coda.samples(mod, params, n.iter=8000, thin=8)

### save model output for post-processing
setwd("~/Google Drive/Git/behavior/output")
dput(out,  paste('Output_sim',n.spec,'_', iter, '.R', sep=''))

iter<-iter+1
} ##end iteration loop




############### post processing of model output ###################################################

library(coda)


######function to get mode out of output. not yet understood
#for continuous data
mfun<-function(x){
  fx<-density(x)
  md<-fx$x[fx$y==max(fx$y)]
  return(md)
}

#for discrete data
mfund<-function(x){
  fx<-table(x)
  if (length(as.numeric(dimnames(fx)[[1]])[fx == max(fx)]) == 1 ) 
    md <-as.numeric(dimnames(fx)[[1]])[fx == max(fx)] 
  
  if (length(as.numeric(dimnames(fx)[[1]])[fx == max(fx)]) > 1 ) 
    md <-sample(as.numeric(dimnames(fx)[[1]])[fx == max(fx)],1) #random draw if more than 1 max
  return(md)
}

##########################################################################################################


###make tables to hold evals
###this requires reading in one output file before running the iter loop

out<-dget(paste('Output_sim',n.spec,'_', 1, '.R', sep=''))

parms<-dimnames(out[[1]])[[2]][-c(1,2)] ###get parameter names


###separate abundance estimates, species level parameters and community parameters
nin<-c(grep('N.', parms),grep('sp.abun', parms))
Nparms<-parms[nin]
indparms<-parms[c(grep('alpha', parms), grep('beta', parms))]
comparms<-parms[-sort(c(nin,c(grep('alpha', parms), grep('sig', parms), grep('beta', parms)) ))]

Ntab<-array(NA,c(n.spec*n.sites+n.spec, niter, 8))
dimnames(Ntab)[[3]]<-c('Mean', 'Mode', 'True', 'AbsBiasMean', 'RelBiasMean', 'AbsBiasMode', 'RelBiasMode', 'CI')
dimnames(Ntab)[[1]]<-Nparms

indtab<-array(NA,c(n.spec*3+2+n.behav*4, niter, 8))
dimnames(indtab)[[3]]<-c('Mean', 'Mode', 'True', 'AbsBiasMean', 'RelBiasMean', 'AbsBiasMode', 'RelBiasMode', 'CI')
dimnames(indtab)[[1]]<-indparms

#dont understand this part yet
comtab<-array(NA,c(length(comparms), niter, 6))
dimnames(comtab)[[3]]<-c('Mean', 'Mode', 'True', 'RelBiasMean', 'RelBiasMode', 'CI')
dimnames(comtab)[[1]]<-comparms

###array for Bp-values
bptab<-array(NA, c(2, niter,2) )
dimnames(bptab)[[3]]<-c('Mean', 'SD')
dimnames(bptab)[[1]]<-c('N', 'Obs')
dimnames(bptab)[[2]]<-1:niter

###vector to check max(Rhat)
Rhat<-NULL


#for (iter in 1:niter){
  
  dat<-dget(paste('Data_sim',n.spec,'_', iter, '.R', sep=''))
  
  ###get iteration-specific species level input values
  l.beta<-dat$l.beta 					# coefficient for spatial covariate on log(expected abundance)
  l.alpha<-dat$l.alpha			# intercept, log(expected abundance)
  s.alpha<-dat$s.alpha				# intercept, log(sigma)
  N.tot<-c(apply(dat$N,1,sum),spec.site)			# total abundance
  
  ####get model results and summarize using coda package
  out<-dget(paste('Output_sim',n.spec,'_', iter, '.R', sep=''))
  sout<-summary(out)
  
  ###get 2.5th and 97.5th percentiles of posteriors
  sqt<-sout[[2]][-c(1,2),c(1,5)]
  
  ###arrange input values in same order as parameters in model output
  inp.ind<-c(l.alpha,l.mu.alpha,l.sig.alpha,s.alpha,s.mu.alpha,s.sig.alpha,l.beta,s.beta)
  #inp.com<-c(s.beta, l.mu.alpha, l.sig.alpha, l.sig.alpha, s.mu.alpha, s.sig.alpha)
  
  #####get posterior mean
  Ntab[,iter, 1]<-sout[[1]][Nparms,1]
  Ntab[,iter, 3]<-N.tot
  indtab[,iter,1]<-sout[[1]][indparms,1]
  indtab[,iter,3]<-inp.ind
  #comtab[,iter,1]<-sout[[1]][comparms,1]
  #comtab[,iter,3]<-inp.com
  
  ### get mode
  mout<-rbind(out[[1]], out[[2]], out[[3]])
  mout<-mout[,-c(1,2)]
  #not looking at the bayesian p value for this part
  
  ###mode
  Ntab[,iter,2]<-apply(mout[,nin],2,mfund)
  indtab[,iter,2]<-apply(mout[,indparms],2,mfun)
  #comtab[,iter,2]<-apply(mout[,comparms],2,mfun)
  
  ###get  bias mean (abs, rel)
  Ntab[,iter, 4]<-Ntab[,iter, 1]-N.tot
  Ntab[,iter, 5]<-Ntab[,iter, 4]/N.tot * 100
  indtab[,iter,4]<-indtab[,iter,1]-inp.ind
  indtab[,iter,5]<-indtab[,iter,4]/inp.ind *100
  
  #comtab[,iter,4]<-(comtab[,iter,1]-inp.com)/inp.com *100
  
  ##exponentiate mu_b1 (which is 0) for relative bias
 # comtab[comparms=='mu_b1',iter,4]<- (exp(comtab[comparms=='mu_b1',iter,1])-exp(inp.com[comparms=='mu_b1']))/
    #exp(inp.com[comparms=='mu_b1']) *100
  
  ###get  bias mode (abs, rel)
  Ntab[,iter, 6]<-Ntab[,iter, 2]-N.tot
  Ntab[,iter, 7]<-Ntab[,iter, 6]/N.tot * 100
  
  indtab[,iter,6]<-indtab[,iter,2]-inp.ind
  indtab[,iter,7]<-indtab[,iter,6]/inp.ind *100
  
  #comtab[,iter,5]<-(comtab[,iter,2]-inp.com)/inp.com *100
  
  #comtab[comparms=='mu_b1',iter,5]<- (exp(comtab[comparms=='mu_b1',iter,2])-exp(inp.com[comparms=='mu_b1']))/
  #  exp(inp.com[comparms=='mu_b1']) *100
  
  ###confidence interval coverage
  Ntab[,iter, 8]<-as.numeric(N.tot >= sqt[Nparms,1] & N.tot <= sqt[Nparms,2])
  indtab[,iter, 8]<-as.numeric(inp.ind >= sqt[indparms,1] & inp.ind <= sqt[indparms,2])
  #comtab[,iter, 6]<-as.numeric(inp.com >= sqt[comparms,1] & inp.com <= sqt[comparms,2])
  
  ##########################################################################
  ## get Bayesian pvalues ##################################################
  
  bptab[,iter,1:2]<-sout[[1]][1:2,1:2]
  
  ##########################################################################
  ### check convergence 
  Rhat[iter]<-max(gelman.diag(out,multivariate=F)$psrf[,1])
  
#} #end iteration loop




