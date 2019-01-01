setwd('U:\\independent studies\\tsbp\\stochastic block model\\fake data')
nome=paste0('fake ',c('data ','psi ','z ','w '),jj,'sim ',ngroup[k],'ng','.csv')
psi.true=read.csv(nome[2],as.is=T)
z.true=read.csv(nome[3],as.is=T)
w.true=read.csv(nome[4],as.is=T)

plot(res$llk,type='l')

zestim=res$z[nrow(res$z),]
k=data.frame(ztrue=z.true,zestim=zestim)
k1=table(k)
ind.loc=numeric()
for (i in 1:nrow(k1)){
  tmp=which(k1[i,]==max(k1[i,]))
  ind.loc=c(ind.loc,colnames(k1)[tmp])
}
ind.loc=as.numeric(ind.loc)

theta.estim=res$theta[nrow(res$theta),]
plot(theta.estim,type='h')
theta.estim[ind.loc];theta.true

rango=range(c(theta.true,theta.estim[ind.loc]))
plot(theta.true,theta.estim[ind.loc],xlim=rango,ylim=rango)
lines(rango,rango)
#---------------------------------------
westim=res$w[nrow(res$w),]
k=data.frame(wtrue=w.true,westim=westim)
k1=table(k)
ind.spp=numeric()
for (i in 1:nrow(k1)){
  tmp=which(k1[i,]==max(k1[i,]))
  ind.spp=c(ind.spp,colnames(k1)[tmp])
}
ind.spp=as.numeric(ind.spp)

phi.estim=res$phi[nrow(res$phi),]
phi.estim[ind.spp];phi.true
plot(phi.estim,type='h')

rango=range(c(phi.estim[ind.spp],phi.true))
plot(phi.true,phi.estim[ind.spp],xlim=rango,ylim=rango)
lines(rango,rango)
#----------------------------------------
tmp=res$psi[nrow(res$psi),]
psi0=matrix(tmp,50,50)
psi=psi0[,ind.spp]
psi1=psi[ind.loc,]
rango=c(0,1)
plot(psi1,data.matrix(psi.true),xlim=rango,ylim=rango)
lines(rango,rango)
