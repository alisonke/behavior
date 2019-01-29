
##########################################################################################################
#### Sollmann, R., Gardner, B., Williams, K.A., Gilbert, A.T. and Veit, R.R. A hierarchical distance #####
#### sampling model to estimate abundance and covariate associations of species and communities	     #####	
####												     #####
#### Appendix S2: R and JAGS code to implement simulation study and case study			     #####
##########################################################################################################


##########################################################################################################
#### Simulation study:										     #####
#### simulate abundance for a community of species, then generate distance sampling data from that   #####
#### and analyze with data-generating model; summarize results across iterations                     #####
##########################################################################################################


library(rjags)

set.seed(2013)

#number of species in the community
n.spec<-30 # 5, 15, 30

# Half-normal detection function
g <- function(x, sig) exp(-x^2/(2*sig^2))

nSites <- 50					# number of line transect surveys
strip.width <- 10 				# strip half-width, w (in this example only one side of the line transect is surveyed)
int.w<-1					# width of distance categories (v)
dist.breaks <- seq(0, strip.width, by=int.w)	# distance break points
nG<-length(dist.breaks)-1			# number of distance categories


############################################################################################################
### detection component: sigma with one covariate (fixed effect)
mu.s<-log(2.5)					# mean of species-level random effect on intercept of sigma
sig.s<- 0.25					# standard deviation of species-level random effect on intercept of sigma
beta.s<- -0.2					# fixed effect of observation covariate on sigma

###look at distribution of sigma intercepts 
hist(exp(rnorm(1000, mu.s, sig.s)))

##detection prob at farthest distance interval for largest sigma
g(10,6)


############################################################################################################
###abundance component, with one covariate (random species level random effect)
#Intercept
mu.lam.alpha<-log(1.5)				# mean of species-level random effect on intercept of log(expected abundance)
sig.lam.alpha<-1				# SD of species-level random effect on intercept of log(expected abundance)
mu.b1<-0					# mean of species-level random effect on coefficient of log(expected abundance)
sig.b1<-0.5					# SD of species-level random effect on coefficient of log(expected abundance)




###########################################################################################################
### begin iterations ######################################################################################
niter<-100					# number of iterations
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

############### JAGS model code ###################################################

#### this code needs to be written to a .txt file called 'Community_DS_Simulations.txt'
#### if given a different name, adjust accordingly in the code above


model{
  
  ###species specific parameters
  
  for (s in 1:spec){
    asig[s]~dnorm(mu_s, tau_s)
    beta1[s]~dnorm(mu_b1,tau_b1)
    alpha[s]~dnorm(mu_a,tau_a)
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
  
  bsig~dnorm(0,0.01)
  
  
  for (s in 1:spec){
    
    for (j in 1:nsites){
      
      sigma[s,j]<-exp(asig[s] + bsig*OBSVAR[j])
      
      f.0[s,j] <- 2 * dnorm(0,0, 1/sigma[s,j]^2)
      
      for(k in 1:nG){
        
        ### approximation to integral of distance function over distance categories
        ### by using mid-point of distance categories
        ### works for point surveys (with appropriate values for pi[k]) 
        ### if used, delete line above starting with 'f.0'
        
        ### p[s,j,k]<- exp( -xg[k]*xg[k]/(2*sigma[s,j]*sigma[s,j]) ) #
        ### fc[s,j,k]<- p[s,j,k]*pi[k]                         
        ### fsc[s,j,k]<- fc[s,j,k]/pcap[s,j]       
        ### fct[s,j,k]<-fsc[s,j,k]/sum(fsc[s,j,1:nG])  
        
        
        ### actual integral over distance categories
        
        up[s,j,k]<-pnorm(db[k+1], 0, 1/sigma[s,j]^2) 
        low[s,j,k]<-pnorm(db[k], 0, 1/sigma[s,j]^2) 
        p[s,j,k]<- 2 * ( up[s,j,k] - low[s,j,k])
        f[s,j,k]<- p[s,j,k]/f.0[s,j]/v                          
        fc[s,j,k]<- f[s,j,k] * pi[k]       
        fct[s,j,k]<-fc[s,j,k]/sum(fc[s,j,1:nG])  
        
      }
      
      
      pcap[s,j]<-sum(fc[s,j,1:nG])    # overall detection probability
      
      lambda[j,s]<- exp(alpha[s] + beta1[s]*v1[j])
      
      ### for a flexible number of covariates, use:
      ### lambda[j,s]<- exp(alpha[s] + inprod(beta1[s,]*v1[j,]))
      ### see seabird application for an example
      
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
  
  
}



###################################################################################################
####################################################################################################

############### post processing of model output ###################################################

library(coda)


######function to get mode out of output
mfun<-function(x){
  fx<-density(x)
  md<-fx$x[fx$y==max(fx$y)]
  return(md)
}

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
out<-dget(paste('Output_spec',n.spec,'_', 1, '.R', sep=''))

parms<-dimnames(out[[1]])[[2]][-c(1,2)] ###get parameter names


###separate abundance estimates, species level parameters and community parameters
nin<-grep('N', parms)
Nparms<-parms[nin]
indparms<-parms[c(grep('alpha', parms), grep('asig', parms), grep('beta1', parms))]
comparms<-parms[-sort(c(nin,c(grep('alpha', parms), grep('asig', parms), grep('beta1', parms)) ))]

Ntab<-array(NA,c(n.spec, niter, 8))
dimnames(Ntab)[[3]]<-c('Mean', 'Mode', 'True', 'AbsBiasMean', 'RelBiasMean', 'AbsBiasMode', 'RelBiasMode', 'CI')
dimnames(Ntab)[[1]]<-Nparms

indtab<-array(NA,c(n.spec*3, niter, 8))
dimnames(indtab)[[3]]<-c('Mean', 'Mode', 'True', 'AbsBiasMean', 'RelBiasMean', 'AbsBiasMode', 'RelBiasMode', 'CI')
dimnames(indtab)[[1]]<-indparms

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


for (iter in 1:niter){
  
  dat<-dget(paste('Data_spec',n.spec,'_', iter, '.R', sep=''))
  
  ###get iteration-specific species level input values
  b1<-dat$b1 					# coefficient for spatial covariate on log(expected abundance)
  lam.alpha<-dat$lam.alpha			# intercept, log(expected abundance)
  s.alpha<-dat$s.alpha				# intercept, log(sigma)
  N.tot<-apply(dat$N,1,sum)			# total abundance
  
  ####get model results and summarize using coda package
  out<-dget(paste('Output_spec',n.spec,'_', iter, '.R', sep=''))
  sout<-summary(out)
  
  ###get 2.5th and 97.5th percentiles of posteriors
  sqt<-sout[[2]][-c(1,2),c(1,5)]
  
  ###arrange input values in same order as parameters in model output
  inp.ind<-c(lam.alpha, s.alpha, b1)
  inp.com<-c(beta.s, mu.lam.alpha, mu.b1, mu.s, sig.lam.alpha, sig.b1, sig.s)
  
  
  #####get posterior mean
  Ntab[,iter, 1]<-sout[[1]][Nparms,1]
  Ntab[,iter, 3]<-N.tot
  indtab[,iter,1]<-sout[[1]][indparms,1]
  indtab[,iter,3]<-inp.ind
  comtab[,iter,1]<-sout[[1]][comparms,1]
  comtab[,iter,3]<-inp.com
  
  ### get mode
  mout<-rbind(out[[1]], out[[2]], out[[3]])
  mout<-mout[,-c(1,2)]
  
  ###mode
  Ntab[,iter,2]<-apply(mout[,nin],2,mfund)
  indtab[,iter,2]<-apply(mout[,indparms],2,mfun)
  comtab[,iter,2]<-apply(mout[,comparms],2,mfun)
  
  ###get  bias mean (abs, rel)
  Ntab[,iter, 4]<-Ntab[,iter, 1]-N.tot
  Ntab[,iter, 5]<-Ntab[,iter, 4]/N.tot * 100
  indtab[,iter,4]<-indtab[,iter,1]-inp.ind
  indtab[,iter,5]<-indtab[,iter,4]/inp.ind *100
  
  comtab[,iter,4]<-(comtab[,iter,1]-inp.com)/inp.com *100
  
  ##exponentiate mu_b1 (which is 0) for relative bias
  comtab[comparms=='mu_b1',iter,4]<- (exp(comtab[comparms=='mu_b1',iter,1])-exp(inp.com[comparms=='mu_b1']))/
    exp(inp.com[comparms=='mu_b1']) *100
  
  ###get  bias mode (abs, rel)
  Ntab[,iter, 6]<-Ntab[,iter, 2]-N.tot
  Ntab[,iter, 7]<-Ntab[,iter, 6]/N.tot * 100
  
  indtab[,iter,6]<-indtab[,iter,2]-inp.ind
  indtab[,iter,7]<-indtab[,iter,6]/inp.ind *100
  
  comtab[,iter,5]<-(comtab[,iter,2]-inp.com)/inp.com *100
  
  comtab[comparms=='mu_b1',iter,5]<- (exp(comtab[comparms=='mu_b1',iter,2])-exp(inp.com[comparms=='mu_b1']))/
    exp(inp.com[comparms=='mu_b1']) *100
  
  ###confidence interval coverage
  Ntab[,iter, 8]<-as.numeric(N.tot >= sqt[Nparms,1] & N.tot <= sqt[Nparms,2])
  indtab[,iter, 8]<-as.numeric(inp.ind >= sqt[indparms,1] & inp.ind <= sqt[indparms,2])
  comtab[,iter, 6]<-as.numeric(inp.com >= sqt[comparms,1] & inp.com <= sqt[comparms,2])
  
  ##########################################################################
  ## get Bayesian pvalues ##################################################
  
  bptab[,iter,1:2]<-sout[[1]][1:2,1:2]
  
  ##########################################################################
  ### check convergence 
  Rhat[iter]<-max(gelman.diag(out)$psrf[,1])
  
} #end iteration loop



##############################################################################################
#### post-post-processing #####################################################################

##### check that all models converged
sum(sumsim[[5]]<=1.1) == niter ##TRUE if all converged

##### if some haven't converged, check which and how high the highest R-hat is
ch<-which(sumsim[[5]] >1.1)
sumsim[[5]][ch]


####################################################################################################
### continue after convergence check ###############################################################


##### get Bpvalue, mean and sd over all iterations

Bp<-matrix(NA, 2,2)
rownames(Bp)<-c('N', 'Obs')
colnames(Bp)<-c('Mean', 'SD')

Bp[1,1]<-mean(sumsim[[4]][1,,1])
Bp[1,2]<-sd(sumsim[[4]][1,,1])
Bp[2,1]<-mean(sumsim[[4]][2,,1])
Bp[2,2]<-sd(sumsim[[4]][2,,1])

write.table(round(Bp, dig=3), paste('Bp_Pois_', n.spec,'.csv', sep=''), sep=';' )



###### get regular summary stats for community level parameters

npar<-dim(sumsim[[3]])[1]

fin<-matrix(NA, nrow=npar, ncol=7)
rownames(fin)<-dimnames(sumsim[[3]])[[1]]
colnames(fin)<-c('AvgMean', 'RMSE', 'Bias','AvgMode', 'BiasMode', 'CIcov', 'TrueValue')


###mean across simulations
fin[,1]<-round(apply(sumsim[[3]][,,1],1, mean),dig=3)


###rmse

resid<-(sumsim[[3]][,,1] - matrix(inp.com, nrow=length(inp.com), ncol=niter))^2
fin[,2]<-round(sqrt(apply(resid,1,sum)/niter), dig=3)

###bias
fin[,3]<-round(apply(sumsim[[3]][,,4],1, mean), dig=3)

### avg mode
fin[,4]<-round(apply(sumsim[[3]][,,2],1, mean),dig=3)

###bias mode
fin[,5]<-round(apply(sumsim[[3]][,,5],1, mean), dig=3)

####CI coverage

fin[,6]<-apply(sumsim[[3]][,,6],1,sum)

###input values
fin[,7]<-inp.com

###write output table
write.csv(fin, paste('CommunityResults_', n.spec, 'spec.csv', sep=''))


###################################################################################################
############ species specific parameters ##########################################################

### total abundance

Nvec<-as.vector(sumsim[[1]][,,1])				# get abundance from all iterations
biasvec<-as.vector(sumsim[[1]][,,4])				# get relative bias
biasvec.a<-as.vector(sumsim[[1]][,,5])				# get absolute bias
Nmodevec<-as.vector(sumsim[[1]][,,2])				# get mode from all iterations
Ntrue<-NtrueX<-as.vector(sumsim[[1]][,,3])			# get true abundance
covvec<-as.vector(sumsim[[1]][,,8])				# get CI coverage from all iterations

###create abundance categories
catb<-c(0,10, 100, 1000, max(Ntrue))

Ncut<-cut(Ntrue,catb,include.lowest=T )

cats<-levels(Ncut)

### get bias and other stats per category
Nout<-matrix(NA, nrow=length(cats), ncol=6)
rownames(Nout)<-c(cats)
colnames(Nout)<-c('MeanEst','MeanTrue','MeanAbsBias','MeanRelBias', 'Coverage','#inCategory')

for (i in 1:length(cats)){
  Nout[i,1]<-mean(Nvec[Ncut==cats[i]])
  Nout[i,2]<-mean(Ntrue[Ncut==cats[i]])
  Nout[i,3]<-mean(biasvec[Ncut==cats[i]])
  Nout[i,4]<-mean(biasvec.a[Ncut==cats[i]])
  Nout[i,5]<-mean(covvec[Ncut==cats[i]])
}

#get number of cases in each abundance category
Nout[,6]<-table(Ncut)

##write output table
write.csv(round(Nout, dig=3),paste('Nresults_', n.spec, '.csv', sep='' ))


#### lambda intercept

indc<-1:n.spec

Nvec<-as.vector(sumsim[[2]][indc,,1])
biasvec<-as.vector(sumsim[[2]][indc,,4])
biasvec.a<-as.vector(sumsim[[2]][indc,,5])
Nmodevec<-as.vector(sumsim[[2]][indc,,2])
Ntrue<-as.vector(sumsim[[2]][indc,,3])
covvec<-as.vector(sumsim[[2]][indc,,8])



###make table with stats by abundance category

Ncut<-cut(NtrueX,catb,include.lowest=T )

cats<-levels(Ncut)

Nout<-matrix(NA, nrow=length(cats), ncol=6)
rownames(Nout)<-c(cats)
colnames(Nout)<-c('MeanEst','MeanTrue','MeanAbsBias','MeanRelBias', 'Coverage','#InCat')

for (i in 1:length(cats)){
  Nout[i,1]<-mean(Nvec[Ncut==cats[i]])
  Nout[i,2]<-mean(Ntrue[Ncut==cats[i]])
  Nout[i,3]<-mean(biasvec[Ncut==cats[i]])
  Nout[i,4]<-mean(biasvec.a[Ncut==cats[i]])
  Nout[i,5]<-mean(covvec[Ncut==cats[i]])
}
Nout[,6]<-table(Ncut)
write.csv(round(Nout, dig=3),paste('AlphaByN_', n.spec, '.csv', sep='' ))



#### second, sigma intercept

###index for where within output the parameter is
indc<-(n.spec+1):(2*n.spec)

Nvec<-as.vector(sumsim[[2]][indc,,1])
biasvec<-as.vector(sumsim[[2]][indc,,4])
biasvec.a<-as.vector(sumsim[[2]][indc,,5])
Nmodevec<-as.vector(sumsim[[2]][indc,,2])
Ntrue<-as.vector(sumsim[[2]][indc,,3])
covvec<-as.vector(sumsim[[2]][indc,,8])

Ncut<-cut(NtrueX,catb,include.lowest=T )

cats<-levels(Ncut)

### get bias and other stats per category
Nout<-matrix(NA, nrow=length(cats), ncol=6)
rownames(Nout)<-c(cats)
colnames(Nout)<-c('MeanEst','MeanTrue','MeanAbsBias','MeanRelBias', 'Coverage','#InCat')

for (i in 1:length(cats)){
  Nout[i,1]<-mean(Nvec[Ncut==cats[i]])
  Nout[i,2]<-mean(Ntrue[Ncut==cats[i]])
  Nout[i,3]<-mean(biasvec[Ncut==cats[i]])
  Nout[i,4]<-mean(biasvec.a[Ncut==cats[i]])
  Nout[i,5]<-mean(covvec[Ncut==cats[i]])
}
Nout[,6]<-table(Ncut)
write.csv(round(Nout, dig=3),paste('Sigmaresults_', n.spec, '.csv', sep='' ))



#### lam coefficient

indc<-(2*n.spec+1):(3*n.spec)

Nvec<-as.vector(sumsim[[2]][indc,,1])
biasvec<-as.vector(sumsim[[2]][indc,,4])
biasvec.a<-as.vector(sumsim[[2]][indc,,5])
Nmodevec<-as.vector(sumsim[[2]][indc,,2])
Ntrue<-as.vector(sumsim[[2]][indc,,3])
covvec<-as.vector(sumsim[[2]][indc,,8])


Ncut<-cut(NtrueX,catb,include.lowest=T )

cats<-levels(Ncut)

### get bias per category
Nout<-matrix(NA, nrow=length(cats), ncol=6)
rownames(Nout)<-c(cats)
colnames(Nout)<-c('MeanEst','MeanTrue','MeanAbsBias','MeanRelBias', 'Coverage', '#InCat')

for (i in 1:length(cats)){
  Nout[i,1]<-mean(Nvec[Ncut==cats[i]])
  Nout[i,2]<-mean(Ntrue[Ncut==cats[i]])
  Nout[i,3]<-mean(biasvec[Ncut==cats[i]])
  Nout[i,4]<-mean(biasvec.a[Ncut==cats[i]])
  Nout[i,5]<-mean(covvec[Ncut==cats[i]])
}
Nout[,6]<-table(Ncut)

write.csv(round(Nout, dig=3),paste('BetaByN_', n.spec, '.csv', sep='' ))



