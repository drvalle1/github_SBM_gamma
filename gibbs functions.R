#' Sample psi parameters
#' 
#' Sample psi parameters from their full conditional distributions
#' 
#' @param dat this matrix has L rows (e.g., locations) and S columns (e.g., species)
#'            and contains the presence-absence data 
#' @param ngroup.loc this is the maximum number of groups for locations
#' @param ngroup.spp this is the maximum number of groups for species
#' 
#' @param return this function returns a matrix of psi parameters
#' @export

sample.psi=function(z,w,dat,ngroup.loc,ngroup.spp){
  #summarize the data by calculating how many observations were assigned to each group
  #for which dat[i,j]=1 and for which dat[i,j]=0
  tmp=getql(z=z-1,w=w-1,dat=dat,ngrloc=ngroup.loc,ngrspp = ngroup.spp)
  
  #generate psi from beta distribution
  tmp1=rbeta(ngroup.loc*ngroup.spp,tmp$nql1+1,tmp$nql0+1)
  matrix(tmp1,ngroup.loc,ngroup.spp)
}

#' Sample theta and vk parameters
#' 
#' Sample vk parameters from their full conditional distributions and convert 
#' these vk parameters into theta parameters
#'
#' @param ngroup.loc this is the maximum number of groups for locations
#' @param gamma.v this is the truncated stick-breaking prior parameter for the 
#'                number of location groups. This value should be between 0 and 1, and
#'                small values enforce more parsimonius results (i.e., fewer groups)
#' @param return this function returns a list containing 2 vectors of parameters (theta and vk)
#' @export

sample.theta=function(ngroup.loc,gamma.v,burnin,gibbs.step,theta,psi,z){
  #re-order thetas. Based on that, re-order z's
  if(gibbs.step<burnin & gibbs.step%%50==0){
    ind=order(theta,decreasing=T)
    theta=theta[ind]
    psi=psi[ind,]
    
    #get z.new
    z.new=z; z.new[]=NA
    for (i in 1:ngroup.loc){
      cond=z==ind[i]
      z.new[cond]=i
    }
    z=z.new
  }
  
  #get the number of locations in each group
  nk=rep(0,ngroup.loc)
  tmp=table(z)
  nk[as.numeric(names(tmp))]=tmp
  
  #sample vk from a beta distribution
  ind=ngroup.loc:1
  invcumsum=cumsum(nk[ind])[ind]
  vk=rbeta(ngroup.loc,nk+1,invcumsum-nk+gamma.v)
  vk[vk>0.9999999]=0.9999999 #for numerical issues
  vk[ngroup.loc]=1
  
  #convert from vk to theta
  theta=convertSBtoNormal(vk)
  
  #to avoid numerical issues
  cond=vk>0.99999
  vk[cond]=0.99999
  
  #output both vk and theta
  list(vk=vk,theta=theta,z=z,psi=psi)
}

#' Sample phi and uk parameters
#' 
#' Sample uk parameters from their full conditional distributions and 
#' then calculate the implied phi parameters
#'
#' @param ngroup.spp this is the maximum number of groups for species
#' @param param list containing the most up-to-date parameter values
#' @param gamma.u this is the truncated stick-breaking prior parameter for the 
#'                number of species groups. This value should be between 0 and 1, and
#'                small values enforce more parsimonius results (i.e., fewer groups)
#' @param return this function returns a list with 2 vectors of parameters (phi and uk)
#' @export
#' 

sample.phi=function(ngroup.spp,gamma.u,burnin,gibbs.step,phi,psi,w){
  #re-order phi. Based on that, re-order w's
  if(gibbs.step<burnin & gibbs.step%%50==0){
    ind=order(phi,decreasing=T)
    phi=phi[ind]
    psi=psi[,ind]
    
    #get z.new
    w.new=w; w.new[]=NA
    for (i in 1:ngroup.spp){
      cond=w==ind[i]
      w.new[cond]=i
    }
    w=w.new
  }
  
  #calculate the number of species in each group
  mk=rep(0,ngroup.spp)
  tmp=table(w)
  mk[as.numeric(names(tmp))]=tmp
  
  #sample uk from a beta distribution
  ind=ngroup.spp:1
  invcumsum=cumsum(mk[ind])[ind]
  uk=rbeta(ngroup.spp,mk+1,invcumsum-mk+gamma.u)
  uk[uk>0.9999999]=0.9999999 #for numerical issues
  uk[ngroup.spp]=1

  #convert uk to phi
  phi=convertSBtoNormal(uk)
  # sum(phi)
  list(uk=uk,phi=phi,w=w,psi=psi)
}
#--------------------------
sample.gamma.u=function(uk,gamma.possib,ngroup.spp){
  ngamma=length(gamma.possib)
  soma=sum(log(1-uk[-ngroup.spp]))
  k=(ngroup.spp-1)*(lgamma(1+gamma.possib)-lgamma(gamma.possib))
  res=k+(gamma.possib-1)*soma
  # sum(dbeta(v[-ngroup.spp],1,gamma.possib[5],log=T))
  res=res-max(res)
  res1=exp(res)
  res2=res1/sum(res1)
  tmp=rmultinom(1,size=1,prob=res2)
  ind=which(tmp==1)
  gamma.possib[ind]
}
#--------------------------
sample.gamma.v=function(vk,gamma.possib,ngroup.loc){
  ngamma=length(gamma.possib)
  soma=sum(log(1-vk[-ngroup.loc]))
  k=(ngroup.loc-1)*(lgamma(1+gamma.possib)-lgamma(gamma.possib))
  res=k+(gamma.possib-1)*soma
  # sum(dbeta(v[-ngroups],1,gamma.possib[5],log=T))
  res=res-max(res)
  res1=exp(res)
  res2=res1/sum(res1)
  tmp=rmultinom(1,size=1,prob=res2)
  ind=which(tmp==1)
  gamma.possib[ind]
}
#--------------------------
sample.z=function(ltheta,dat,dat1m,lpsi,l1mpsi,ngroup.loc,ngroup.spp,nloc,nspp,w,z){
  #probability for groups that exist
  lprob.exist=matrix(NA,nloc,ngroup.loc)
  for (i in 1:ngroup.loc){
    lpsi1=matrix(lpsi[i,w],nloc,nspp,byrow=T)
    l1mpsi1=matrix(l1mpsi[i,w],nloc,nspp,byrow=T)
    lprob.exist[,i]=ltheta[i]+rowSums(dat*lpsi1+dat1m*l1mpsi1)
  }
  
  #probability for groups that do not exist yet
  tmp=rep(0,nloc)
  k=lgamma(2+nspp)
  for (i in 1:ngroup.spp){
    cond=w==i
    if (sum(cond)==0) nik=rep(0,nloc)
    if (sum(cond)==1) nik=dat[,cond]
    if (sum(cond)> 1) nik=rowSums(dat[,cond])
    tmp=tmp+lgamma(nik+1)+lgamma(nspp-nik+1)-k
  }
  lprob.nexist=matrix(ltheta,nloc,ngroup.loc,byrow=T)+matrix(tmp,nloc,ngroup.loc)
  
  #sample z
  tab=rep(0,ngroup.loc)
  tmp=table(z)
  tab[as.numeric(names(tmp))]=tmp
  
  for (i in 1:nloc){
    tab[z[i]]=tab[z[i]]-1
    lprob=rep(NA,ngroup.loc)
    cond=tab==0
    lprob[ cond]=lprob.nexist[i,cond]
    lprob[!cond]=lprob.exist[i,!cond]
    
    #get normalized probs
    tmp1=lprob-max(lprob) #for numerical stability
    tmp2=exp(tmp1) #exponentiate log probability
    prob=tmp2/sum(tmp2) #normalize to sum to 1
    
    #draw from multinomial distrib
    ind=rmultinom(1,size=1,prob=prob)
    ind1=which(ind==1)
    z[i]=ind1
    tab[ind1]=tab[ind1]+1
  }
  z
}
#--------------------------
sample.w=function(lphi,dat,dat1m,lpsi,l1mpsi,ngroup.spp,ngroup.loc,nloc,nspp,w,z){
  #probability for groups that exist
  lprob.exist=matrix(NA,nspp,ngroup.spp)
  for (i in 1:ngroup.spp){
    lpsi1=matrix(lpsi[z,i],nloc,nspp)
    l1mpsi1=matrix(l1mpsi[z,i],nloc,nspp)
    lprob.exist[,i]=lphi[i]+colSums(dat*lpsi1+dat1m*l1mpsi1)
  }
  
  #probability for groups that do not exist yet
  tmp=rep(0,nspp)
  k=lgamma(2+nloc)
  for (i in 1:ngroup.loc){
    cond=z==i
    if (sum(cond)==0) njk=rep(0,nspp)
    if (sum(cond)==1) njk=dat[cond,]
    if (sum(cond)> 1) njk=colSums(dat[cond,])
    tmp=tmp+lgamma(njk+1)+lgamma(nloc-njk+1)-k
  }
  lprob.nexist=matrix(lphi,nspp,ngroup.spp,byrow=T)+matrix(tmp,nspp,ngroup.spp)
  
  #sample w
  tab=rep(0,ngroup.spp)
  tmp=table(w)
  tab[as.numeric(names(tmp))]=tmp
  
  for (i in 1:nspp){
    tab[w[i]]=tab[w[i]]-1
    lprob=rep(NA,ngroup.spp)
    cond=tab==0
    lprob[ cond]=lprob.nexist[i,cond]
    lprob[!cond]=lprob.exist[i,!cond]
    
    #get normalized probs
    tmp1=lprob-max(lprob) #for numerical stability
    tmp2=exp(tmp1) #exponentiate log probability
    prob=tmp2/sum(tmp2) #normalize to sum to 1
    
    #draw from multinomial distrib
    ind=rmultinom(1,size=1,prob=prob)
    ind1=which(ind==1)
    w[i]=ind1
    tab[ind1]=tab[ind1]+1
  }
  w
}