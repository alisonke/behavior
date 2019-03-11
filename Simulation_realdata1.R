###try some real data
##for this exercise, use first only behavior observed. use LU as forst variable

Raw<-read.csv("~/Google Drive/FuturAves/Occupancy_Model/data/JZBirdSurveys_2017_rep.csv")

#remove those out of range
Raw<-Raw[Raw$Distance!="+",]

#find most observed species - RNWR, rufous-naped wren, indo
Raw<-Raw[Raw$Species=="RNWR",]

#repeat rows by number of observations
library(dplyr)
n.times <- Raw$Number 
repped<-Raw[rep(seq_len(nrow(Raw)), n.times),]

wren<-repped[,c(1,2,3,4,5,6,11,12,13,14,15),]

wren$B<-c(rep(NA,nrow(wren)))
wren$B[wren$Behavior %in% c("AF","BG","BP","E","FO","FV","OF","H","GR","PS","FE","SP")]

wren$B[wren$Behavior %in% c("AF","BG","BP","E","FO","FV","OF","H","GR","PS","FE","SP","SE","ST")]<-"Eat"

wren$B[wren$Behavior %in% c("S","C","A","D")]<-"Voc"

wren$B[wren$Behavior %in% c("CO","FY","M","NB","NS","PB")]<-"Mate"

wren$B[wren$Behavior %in% c("SA","MA","MO")]<-"Agg"




##try again with dove
Raw<-Raw[Raw$Species=="INDO",]

#repeat rows by number of observations
library(dplyr)
n.times <- Raw$Number 
repped<-Raw[rep(seq_len(nrow(Raw)), n.times),]

dove<-repped[,c(1,2,3,4,5,6,11,12,13,14,15),]

dove$B<-c(rep(NA,nrow(dove)))
dove$B[dove$Behavior %in% c("AF","BG","BP","E","FO","FV","OF","H","GR","PS","FE","SP")]

dove$B[dove$Behavior %in% c("AF","BG","BP","E","FO","FV","OF","H","GR","PS","FE","SP","SE","ST")]<-"Eat"

dove$B[dove$Behavior %in% c("S","C","A","D")]<-"Voc"

dove$B[dove$Behavior %in% c("CO","FY","M","NB","NS","PB")]<-"Mate"

dove$B[dove$Behavior %in% c("SA","MA","MO")]<-"Agg"

dove$B[dove$Behavior %in% c("F")]<-"Fly"

dove$B[dove$Behavior %in% c("PR","P","R")]<-"Pass"

###make Y matrix. need to have number for each site, behavior, and visit

#sites with 3 visits only #NARA, BARR, PANI, FRIO, MONN, PELP, PELC, VIEN, PURA, HOND, LOMA, MOOR

dove3<-dove[dove$Site %in% c("NARA","BARR","PANI","FRIO","MONN","PELP","PELC","VIEN","PURA","LOMA","MOOR"),]

Y<-array(rep(NA,66*3*4),dim=c(66,4,3))
sites<-c("BARR","FRIO","LOMA","MONN","MOOR","NARA","PANI","PELC","PELP","PURA","VIEN")

for(i in 1:11){
  for(point in 1:6){
  for(j in 1:4){
    for(k in 1:3){
      specifysite<-dove[dove$Site==sites[i],]
      specifypoint<-specifysite[specifysite$Point==point,]
      specifyvisit<-specifypoint[specifypoint$Rep==k,]
      Y[(i-1)*6+point,,k]<-summary(factor(specifyvisit$B,levels=c("Agg","Eat","Pass","Voc")))}
    }
  }
  }

####get forest covers

forest<-read.csv("~/Google Drive/FuturAves/DannyData/use/2017_Updates/LandUseRings_2017.csv")
fc<-forest %>% filter(buffer==50) 
fc<-fc[-c(1:6,13:42,49:60,73:78,91:102,127:144),]

forest<-fc$prop

#####model

data <- list(nSites = 66, nVisits = 3, forest = forest, nBehav = 4, Y=Y)

model.string <-"
model {
#### Priors 

#forest effect on abundance
beta ~ dnorm(0,0.001)

#forest effect on detection
Betadet ~ dnorm(0,0.0001)

#forest effect on behavior
for (j in 1 : nBehav) { b.beta[j] ~ dnorm(0, 0.0001)} # vague priors for behaviors

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
inits <- function(){list(Betadet=rnorm(1),N=rep(max(Y)*4,66),beta= rnorm(1))}
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

