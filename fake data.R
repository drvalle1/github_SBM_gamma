rm(list=ls(all=TRUE))
set.seed(3)

setwd('U:\\GIT_models\\github_SBM_gamma')
ngroup.loc=5
ngroup.spp=3

#get parameters
tmp=runif(ngroup.loc)
theta.true=theta=tmp/sum(tmp)
tmp=runif(ngroup.spp)
phi.true=phi=tmp/sum(tmp)

set.seed(4)
psi=matrix(c(0.05,0.5,0.95,
             0.5,0.05,0.95,
             0.05,0.95,0.5,
             0.5,0.95,0.05,
             0.1,0.5,0.05),ngroup.loc,ngroup.spp,byrow=T)
psi.true=psi

#get latent variables
nloc=1000
tmp=rmultinom(1,size=nloc,prob=theta)
tmp1=rep(1:ngroup.loc,times=tmp)
z.true=z=tmp1 #if not scrambled
# z=sample(tmp1,nind); 

nspp=50
tmp=rmultinom(1,size=nspp,prob=phi)
tmp1=rep(1:ngroup.spp,times=tmp)
w.true=w=tmp1 #if not scrambled
# w=sample(tmp1,nquest)

#generate data
y=matrix(NA,nloc,nspp)
for (i in 1:nloc){
  for (j in 1:nspp){
    y[i,j]=rbinom(1,size=1,prob=psi[z[i],w[j]])    
  }
}
image(y)

write.csv(y,'fake data.csv',row.names=F)