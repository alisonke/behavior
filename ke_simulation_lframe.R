###Frishkoff framework simulation

#Generate data

#################################abundance process############################

#set number of sites, visits, behaviors
site<-10 #number of sites
visit<-3 #number of visits
behav<-4 #number of behaviors

#simulate forest covers
forest<-c(rep(0,site)) #empty forest vector
for (j in 1:site){
forest[j]<-runif(10,0,1)} #forest cover for each site

#make up a beta for relationship between forest and abundance
beta<-4

#simulate actual abundance data for single species
lambda<-c(rep(0,site)) #empty lambda vector
N<-c(rep(0,site)) #empty N vector
for (j in 1:site){
  lambda[j]<-exp(beta*forest[j]) #lambda for each site
    N[j]<-rpois(1,lambda[j]) }#N for each site

#simulate actual animals doing behaviors. each matrix is a different visit.
##create underlying behavior probabilities 
#Is this per species? how will we get the forest cover effect on behavior? We want the underlying probabilities to change because our question is whether behavhiors switch between habitats. Alternative: make separate models for just agriculture and just forest. Then compare those underlying behaviors for each species.

behav.prob<-c(.2,.1,.6,.1)
D<-array(c(rep(0,site*behav*visit)),dim=c(site,behav,visit)) #behavior matrix

for (j in 1:site){
 for (k in 1:visit){
    D[j,1:behav,k] <- rmultinom(1,N[j],behav.prob)}}


###########################detection process##############################

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


