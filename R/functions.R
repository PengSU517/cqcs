cqcor<-function(y,d,x,tau){
  n<-length(y)
  y<-as.numeric(y)
  x<-as.numeric(x)
  Te<-as.data.frame(cbind(y,d,x,P = NA))[order(y),]
  lambda<-((n-c(1:n))/(n-c(1:n)+1))^Te$d
  S<-cumprod(lambda)
  for(i in 1:n)
  { if((Te$d[i]==0)&((1-S[i])<tau)) Te$P[i]<-(tau-(1-S[i]))/S[i]
  else Te$P[i]<-(Te$y[i]<quantile(Te$y,tau)) }
  cqcov<-sum((tau-Te$P)*Te$x)/n
  covx<-sum((Te$x-mean(Te$x))^2)/n
  cqcov/sqrt((tau-tau^2)*covx)}

cqc<-function(y,d,x,tau){
  n<-length(y)
  cpcvec<-matrix(NA,nrow=1,ncol=dim(x)[2])
  for(k in 1:dim(x)[2]) cpcvec[k]<-cqcor(y,d,x[,k],tau)
  order(abs(cpcvec),decreasing = T)}

cqr<-function(y,d,x,tau){
  n<-length(y)
  p<-dim(x)[2]
  y<-as.numeric(y)
  linf<-10*max(y)
  x<-as.matrix(x)
  Te<-cbind(y,d,w = 1,x)[order(y),]

  for(i in 1:n){
    if(Te[i,2]==0){
      distance<-abs(x-x[rep(i,n),])/5#?˴?????
      weight<-15/16*(1-apply(distance^2,1,sum))^2#*(apply(distance,1,max)<1)
      oneminuslambda<-((sum(weight)-cumsum(weight))/
                         (sum(weight)-cumsum(c(0,weight[-n]))))^Te[1:n,2]
      S<-cumprod(oneminuslambda)
      Te[i,3]<-max(0,(tau-1+S[i])/S[i])

      Te<-rbind(Te,c(linf,Te[i,2],min(1,(1-tau)/S[i]),Te[i,4:(p+3)]))
    }}
  yy<-as.numeric(Te[,1])
  xx<-as.matrix(Te[,4:(p+3)])
  ww<-as.numeric(Te[,3])
  nn<-length(yy)
  fit1<-quantreg::rq(yy~xx,weights = ww,tau = tau)
  betahat<-fit1$coefficients
  betahat
}


cqrloss<-function(y,d,x,tau){
  n<-length(y)
  p<-dim(x)[2]
  y<-as.numeric(y)
  linf<-10*max(y)
  x<-as.matrix(x)
  Te<-cbind(y,d,w = 1,x)[order(y),]

  for(i in 1:n){
    if(Te[i,2]==0){
      distance<-abs(x-x[rep(i,n),])/5#?˴?????
      weight<-15/16*(1-apply(distance^2,1,sum))^2#*(apply(distance,1,max)<1)
      oneminuslambda<-((sum(weight)-cumsum(weight))/
                         (sum(weight)-cumsum(c(0,weight[-n]))))^Te[1:n,2]
      S<-cumprod(oneminuslambda)
      Te[i,3]<-max(0,(tau-1+S[i])/S[i])

      Te<-rbind(Te,c(linf,Te[i,2],min(1,(1-tau)/S[i]),Te[i,4:(p+3)]))
    }}
  yy<-as.numeric(Te[,1])
  xx<-as.matrix(Te[,4:(p+3)])
  ww<-as.numeric(Te[,3])
  nn<-length(yy)
  fit1<-quantreg::rq(yy~xx,weights = ww,tau = tau)
  betahat<-fit1$coefficients
  loss<-mean((yy-cbind(1,xx)%*%betahat)*(tau-((yy-cbind(1,xx)%*%betahat)<0)))
  loss
}


lossarray<-function(y,d,x,tau,cqcs){
  n<-length(y)
  xxx<-x[,cqcs]
  d<-dim(xxx)[2]
  selectedlabel<-vector(mode="numeric",length=0)
  selectedloss<-vector(mode="numeric",length=0)
  for(l in 1:d){
    possiblabel<-matrix(Inf,nrow=d,ncol=1)
    for(k in 1:d){
      if(!k%in%selectedlabel){
        possiblabel[k]<-cqrloss(y,d,as.matrix(xxx[,union(k,selectedlabel)]),tau)
      }
    }
    selectedlabel[l]<-which.min(possiblabel)
    selectedloss[l]<-min(possiblabel)
  }
  rbind(label=cqcs[selectedlabel],loss=selectedloss)
}


simulation<-function(m,n,p,tau,mu,mean,sigma,beta,disfunc,quantfunc,heter){

  AICtp<- vector(mode="numeric",length = m)
  BICtp<- vector(mode="numeric",length = m)
  EBIC1tp<- vector(mode="numeric",length = m)
  EBIC2tp<- vector(mode="numeric",length = m)


  AICfp<- vector(mode="numeric",length = m)
  BICfp<- vector(mode="numeric",length = m)
  EBIC1fp<- vector(mode="numeric",length = m)
  EBIC2fp<- vector(mode="numeric",length = m)


  beta<-c(3,3,0,0,3,rep(0,p-5))
  mean <- rep(0,p)
  sigma<-diag(rep(1,p))

  for(k in 1:m)
  {
    set.seed(k)
    x<-MASS::mvrnorm(n, mean, sigma)
    eps<-(disfunc(n,0,1)-quantfunc(tau))*(1+heter*exp(x[,1]))+quantfunc(tau)
    t<-x%*%beta+eps
    c<-rnorm(n,mu,5)
    y<-pmin(t,c)
    d<-as.numeric(t<c)

    ordercqc<-cqc(y,d,x,tau)
    cqcs<-ordercqc[1:(n/log(n))]
    fit<-lossarray(y,d,x,tau,cqcs)

    AIC<-fit[2,]+seq(1:(n/log(n)))*2/n
    BIC<-fit[2,]+seq(1:(n/log(n)))*log(n)/n
    EBIC1<-fit[2,]+seq(1:(n/log(n)))*log(n/log(n))*log(n)/n
    EBIC2<-fit[2,]+seq(1:(n/log(n)))*log(log(n/log(n)))*log(n)/n

    AICs<-fit[1,][1:which.min(AIC)]
    BICs<-fit[1,][1:which.min(BIC)]
    EBIC1s<-fit[1,][1:which.min(EBIC1)]
    EBIC2s<-fit[1,][1:which.min(EBIC2)]

    AICtp[k]<-1%in%AICs+2%in%AICs+5%in%AICs
    AICfp[k]<-length(AICs)-AICtp[k]
    BICtp[k]<-1%in%BICs+2%in%BICs+5%in%BICs
    BICfp[k]<-length(BICs)-BICtp[k]
    EBIC1tp[k]<-1%in%EBIC1s+2%in%EBIC1s+5%in%EBIC1s
    EBIC1fp[k]<-length(EBIC1s)-EBIC1tp[k]
    EBIC2tp[k]<-1%in%EBIC2s+2%in%EBIC2s+5%in%EBIC2s
    EBIC2fp[k]<-length(EBIC2s)-EBIC2tp[k]
  }
  c(mean(AICtp),mean(BICtp),mean(EBIC1tp),mean(EBIC2tp),
    mean(AICfp),mean(BICfp),mean(EBIC1fp),mean(EBIC2fp))
}
