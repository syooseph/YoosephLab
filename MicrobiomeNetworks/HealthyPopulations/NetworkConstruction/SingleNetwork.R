library(mvtnorm)
library(Matrix)
library(matrixcalc)
library(corpcor)
library(glasso)
library(cvTools)
library(huge)
library(mclust)
library(bayesm)

#X is Nxd sample taxa matrix

SingleNetwork<-function(X,F,nrho=50,beta=0.05,bsrep=50,bsthres=0.8){
  if(missing(X)){
    if(missing(F)){
      stop("Error: Need to provide input either through parameter X or parameter F")
    }else{
      X=as.matrix(read.csv(F,header=TRUE,sep=","))
      #      print(c("Dim X",dim(X)))
    }
  }
  if(nrho<=0 || (as.integer(nrho)!=nrho)){
    stop("Error: nrho should be positive integer - length of rho list")
  }
  if(bsrep<=0 || (as.integer(bsrep)!=bsrep)){
    stop("Error: bsrep should be positive integer for number of bootstrap replicates")
  }
  if(bsthres<0 || bsthres>1){
    stop("Error: bsthres should be between 0 and 1")
  }
  if(beta<0 || beta>1){
    stop("Error: beta should be between 0 and 1")
  }
  
  X<-as.matrix(X)
  taxa=colnames(X)
  
  result=huge(x=X,method="glasso",nlambda=nrho,lambda.min.ratio = 0.01,cov.output=TRUE)
  chosen=huge.select(est=result,criterion="stars",stars.thresh=beta)
  boot_res=bootstrap_glasso(X=t(X),Omega=chosen$opt.icov,rho=chosen$opt.lambda,rep=bsrep,threshold=bsthres)
  pboot=mypcor(boot_res)
  
  C=chosen$opt.cov
  O=chosen$opt.icov
  
  rownames(C)=taxa
  colnames(C)=taxa
  rownames(O)=taxa
  colnames(O)=taxa
  rownames(boot_res)=taxa
  colnames(boot_res)=taxa
  rownames(pboot)=taxa
  colnames(pboot)=taxa
  
  return(list("X"=X,"Sigma"=C,"Omega"=O,"BootOmega"=boot_res,"PCorBootOmega"=pboot,"rho"=chosen$opt.lambda,
              "nrho"=nrho,"lrho"=result$lambda,"beta"=beta,"bsrep"=bsrep,"bsthres"=bsthres))
}

bootstrap_glasso<-function(X,Omega,rho,threshold,rep){
  print(c("Bootstrap"))
  N=ncol(X)
  d=nrow(X)
  cnt=matrix(0,ncol=d,nrow=d)
  for(k in 1:rep){
    rd=sample.int(n=N,size=N,replace=TRUE)
    Xs=X[,rd]
    res.glasso = huge(x=t(Xs),lambda=rho,method="glasso",cov.output=TRUE,verbose=FALSE)
    cnt = cnt + res.glasso$path[[1]]
  }
  cnt = cnt/rep
  
  icov=matrix(0,ncol=d,nrow=d)
  diag(icov)=diag(Omega)
  for(i in 1:d){
    for(j in 1:d){
      if(cnt[i,j]>=threshold){
        icov[i,j]=Omega[i,j]
      }
    }
  }
  if(!isSymmetric(icov)){
    icov = (icov + t(icov))/2
  }
  
  return(icov)
}

mypcor<-function(Omega){
  n=ncol(Omega)
  s=sum(Omega==t(Omega))
  if(s!=(n^2)){
    Omega=(Omega+t(Omega))/2
  }
  pcor=matrix(0,ncol=n,nrow=n)
  for(i in 1:n){
    for(j in i:n){
      pcor[i,j] = -(Omega[i,j]/sqrt(Omega[i,i]*Omega[j,j]))
      pcor[j,i] = pcor[i,j]
    }
  }
  return(pcor)
}

myfscore<-function(T,P,thresh=1e-3){
  if(NROW(T)!=NCOL(T) || NROW(P)!=NCOL(P) || NROW(T)!=NROW(P)){
    stop ("Error in myfscore function")
  }
  Tc=1*(abs(T) >= thresh)
  Pc=1*(abs(P) >= thresh)
  n=NROW(T)
  diag(Tc)=0
  diag(Pc)=0
  TP=0 #True positives
  FP=0 #False positives
  Pos=0 #Total number of real positives (TP + FN)
  
  TP = sum((Tc+Pc)==2)/2
  FP = sum((Pc-Tc)==1)/2
  Pos = sum(Tc)/2
  
  #Recall or Sensitivity
  recall = TP/Pos
  
  #Precision or Positive predictive value
  precision = TP/(TP+FP)
  
  Fscore = 2*recall*precision/(recall + precision)
  
  return(list("recall"=recall,"precision"=precision,"fscore"=Fscore))
}
