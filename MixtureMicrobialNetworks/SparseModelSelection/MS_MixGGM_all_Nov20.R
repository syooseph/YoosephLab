library(mvtnorm)
library(Matrix)
library(matrixcalc)
library(corpcor)
library(glasso)
library(huge)
library(mclust)
#library(cvTools)
#library(bayesm)

# Using variational approach to identify number of components for a mixture of gaussian graphical models
# Can run in sparse or non-sparse mode
# Sparse mode: pcor, GLASSO (EBIC, CV), TIGER (EBIC, CV)
# Method: 0 (non-sparse), 1 (pcor), 2 (glasso EBIC), 3 (glasso CV), 4 (tiger ebic), 5 (tiger cv)

# M is number of components
# N is number of data points
# d is dimension of each data point; precision matrix is d x d

# X_dxN is data matrix
# Qm linked with mu and Tu
# QT linked with vT and VT
# pbar_MxN, pin_MxN, mu = list(), Tu =list(), vT=list(), VT=list
# Esij = pin, Emui = mu, EmimiT = list(), ETi = list(), ElnTi = list()

#comments are below
#for scale matrix V as diagonal matrix, larger values on diagonal tend to make Precision Matrix (IW in rwishart) to have smaller entries (i.e. making it sparse); and vice versa
#for v degrees of freedom, larger values of v (for fixed scale matrix V) results in smaller values in Precision Matrix (IW in Wishart)
#having very small values in Prec Mat seems to bias towards much smaller number of clusters?
#for beta to be small value (say 1e-6) means broad prior over mean mu

#In the sparse network version, we apply L1-penalty to invVT (that is, VT^{-1})
#We compute vT and VT using GLASSO (with STARS) and numerical optimization

best_model<-function(X,F,M,beta=1e-6,v,V,niter=30,seed=1001,start=10,method=0,alpha=0.05){
  if(missing(X)){
    if(missing(F)){
      stop("Error: Need to provide input either through parameter X or parameter F")
    }else{
      X=t(as.matrix(read.csv(F,header=TRUE,sep=",")))
#      print(c("Dim X",dim(X)))
    }
  }
  if(start<=0 || (as.integer(start)!=start)){
    stop("Error: number of starting points should be positive integer")
  }
  if(ncol(V)!=nrow(X) || nrow(V)!=nrow(X)){
    stop("Error: V should be square matrix with dimension d, where d is nrow(X)")
  }
  if(v<nrow(X)){
    stop("Error: v should be at least d, where d is nrow(X)")
  }
  if (!(method %in% 0:5)){
    stop("Error: incorrect method chosen; 0 (non-sparse), 1 (pcor), 2 (glasso EBIC), 3 (glasso CV), 4 (tiger ebic), 5 (tiger cv)")
  }
  
  X<-as.matrix(X)
  d=nrow(X)
  N=ncol(X)

  set.seed(seed)    
  
  print(c("Method", method))
  start.time=as.numeric(format(Sys.time(),"%s"))
  
  #run k-means
  print(c("Running K-means"))
  member=kmeans(x=t(X),centers=M,iter.max=100)$cluster
  result=model_select(X=X,M=M,v=v,V=V,beta=beta,niter=niter,member=member,method=method,alpha=alpha)
  print(c("K-means likelihood",result$likelihood))
  print(c("K-means Pi",result$Pi))
  print(c("Iterations",result$iterations))

  #run random starting points using seed  
  for(i in 1:start){
    print(c("Random start point iteration",i))
    member=sample.int(M,size=N,replace=TRUE)
    while(length(table(member))!=M){
      member=sample.int(M,size=N,replace=TRUE)
    }
    result_next=model_select(X=X,M=M,v=v,V=V,beta=beta,niter=niter,member=member,method=method,alpha=alpha)
    print(c("Iteration Likelihood",result_next$likelihood))
    print(c("Iteration Pi",result_next$Pi))
    print(c("Iterations",result_next$iterations))
    
    if(result_next$likelihood>result$likelihood){
      result=result_next
      print(c("Updated"))
    }
  }

  end.time=as.numeric(format(Sys.time(),"%s"))
  
  result$X=X
  result$M=M
  result$beta=beta
  result$v=v
  result$V=V
  result$niter=niter
  result$seed=seed
  result$start=start
  result$method=method
  result$alpha=alpha
  result$time=round(end.time - start.time)

  print(c("Total time", result$time))
  
  return(result)
}

model_select<-function(X,M,beta,v,V,niter,member,method,alpha){
  X<-as.matrix(X)
  d=nrow(X)
  N=ncol(X)
  
  if(M <= 0){
    stop("Error: provide positive integer for M")
  }
  
  if(N < 4*M){
    stop("Error: too few data points")
  }
  
  if(length(member)!=N){
    print(c("Error", N, length(member)))
    stop("Error: N != length(member)")
  }

  listv=list()
  listV=list()
  mu=list()
  Tu=list()
  vT=list()
  VT=list()
  invVT=list()
  ETi=list()
  EmimiT=list()
  ElnTi=list()
  
  const1=0.001
  for(m in 1:M){
    listv[[m]] = v
    listV[[m]] = V
    mu[[m]] = rep(0,d)
    EmimiT[[m]]=matrix(rep(0,d*d),nrow=d)
    VT[[m]] = const1*diag(d)
    invVT[[m]] = (1/const1)*diag(d)
  }

  #initializing pbar matrix
  pin=matrix(rep(1e-5,M*N),nrow=M)
  for(l in 1:M){
    ind=which(member==l)
    pin[l,ind]=1
  }
  
  for(i in 1:N){
    csum=sum(pin[,i])
    pin[,i]=pin[,i]/csum
  }
  
#  print(c("cls",apply(pin,2,function(x) which.max(x))))

  Pi=compute_Pi(M=M,N=N,pin=pin)
#  print(c("Pi",Pi))
  
  for(m in 1:M){
    ind = which(member==m)
    if(length(ind)>=4){
      Xs = X[,ind]
      mu[[m]] = rowMeans(Xs)
      EmimiT[[m]]=mu[[m]]%*%t(mu[[m]])
    }
  }
  
  temp=compute_vT_VT_invVT(v=v,pin=pin,M=M,V=V,X=X,N=N,mu=mu,EmimiT=EmimiT,method=method,alpha=alpha)
  vT=temp$vT
  VT=temp$VT
  invVT=temp$invVT
  
#  vT=compute_vT(v=v,pin=pin,M=M)
#  VT=compute_VT(V=V,X=X,M=M,N=N,pin=pin,mu=mu,EmimiT=EmimiT)
#  for(m in 1:M){
#    invVT[[m]]=solve(VT[[m]])
#  }
  
  ETi=compute_ETi(vT=vT,invVT=invVT,M=M)
  Tu=compute_Tu(beta=beta,ETi=ETi,pin=pin,M=M,d=d)
  mu=compute_mu(Tu=Tu,ETi=ETi,M=M,N=N,X=X,pin=pin)
  EmimiT=compute_EmimiT(Tu=Tu,mu=mu,M=M)
  ElnTi=compute_ElnTi(VT=VT,vT=vT,d=d,M=M)
    
  pbar=compute_pbar(Pi=Pi,ElnTi=ElnTi,X=X,ETi=ETi,mu=mu,EmimiT=EmimiT,M=M,N=N)
  pin=compute_pin(pbar)
  Pi=compute_Pi(M=M,N=N,pin=pin)
 # print(c("cls",apply(pin,2,function(x) which.max(x))))
 # print(c("Pi",Pi))
  print ("All Initialized")
  
  temp=compute_vT_VT_invVT(v=v,pin=pin,M=M,V=V,X=X,N=N,mu=mu,EmimiT=EmimiT,method=method,alpha=alpha)
  vT=temp$vT
  VT=temp$VT
  invVT=temp$invVT
  ETi=compute_ETi(vT=vT,invVT=invVT,M=M)
  Tu=compute_Tu(beta=beta,ETi=ETi,pin=pin,M=M,d=d)
  mu=compute_mu(Tu=Tu,ETi=ETi,M=M,N=N,X=X,pin=pin)
  EmimiT=compute_EmimiT(Tu=Tu,mu=mu,M=M)
  ElnTi=compute_ElnTi(VT=VT,vT=vT,d=d,M=M)
  pbar=compute_pbar(Pi=Pi,ElnTi=ElnTi,X=X,ETi=ETi,mu=mu,EmimiT=EmimiT,M=M,N=N)
  pin=compute_pin(pbar)
  Pi=compute_Pi(M=M,N=N,pin=pin)
#  print(c("Pi",Pi))
  
  likelihood=compute_likelihood(X=X,M=M,N=N,d=d,pin=pin,Pi=Pi,mu=mu,Tu=Tu,vT=vT,VT=VT,invVT=invVT,ElnTi=ElnTi,EmimiT=EmimiT,ETi=ETi,beta=beta,v=v,V=V,lv=listv,lV=listV)
  
  curr_pbar=pbar
  curr_pin=pin
  curr_Pi=Pi
  curr_mu=mu
  curr_Tu=Tu
  curr_vT=vT
  curr_VT=VT
  curr_invVT=invVT
  curr_ElnTi=ElnTi
  curr_ETi=ETi
  curr_EmimiT=EmimiT
  curr_likelihood=likelihood
  
  loop=TRUE
  iter=0
  while(loop && iter<niter){
    iter = iter + 1
    
    #E-step
    temp=compute_vT_VT_invVT(v=v,pin=pin,M=M,V=V,X=X,N=N,mu=mu,EmimiT=EmimiT,method=method,alpha=alpha)
    vT=temp$vT
    VT=temp$VT
    invVT=temp$invVT
    ETi=compute_ETi(vT=vT,invVT=invVT,M=M)
    Tu=compute_Tu(beta=beta,ETi=ETi,pin=pin,M=M,d=d)
    mu=compute_mu(Tu=Tu,ETi=ETi,M=M,N=N,X=X,pin=pin)
    EmimiT=compute_EmimiT(Tu=Tu,mu=mu,M=M)
    ElnTi=compute_ElnTi(VT=VT,vT=vT,d=d,M=M)
    pbar=compute_pbar(Pi=Pi,ElnTi=ElnTi,X=X,ETi=ETi,mu=mu,EmimiT=EmimiT,M=M,N=N)
    pin=compute_pin(pbar)
    
    #M step
    Pi=compute_Pi(M=M,N=N,pin=pin)
    
#    print(c("Pi",Pi))
    
    likelihood=compute_likelihood(X=X,M=M,N=N,d=d,pin=pin,Pi=Pi,mu=mu,Tu=Tu,vT=vT,VT=VT,invVT=invVT,ElnTi=ElnTi,EmimiT=EmimiT,ETi=ETi,beta=beta,v=v,V=V,lv=listv,lV=listV)
    if(likelihood > curr_likelihood){
      curr_pbar=pbar
      curr_pin=pin
      curr_Pi=Pi
      curr_mu=mu
      curr_Tu=Tu
      curr_vT=vT
      curr_VT=VT
      curr_invVT=invVT
      curr_ElnTi=ElnTi
      curr_ETi=ETi
      curr_EmimiT=EmimiT
      curr_likelihood=likelihood
    }else{
      loop=FALSE
    }
  }
  
  #Results in curr_ variables
  result=list("likelihood"=curr_likelihood,"Pi"=curr_Pi,"pbar"=curr_pbar,"pin"=curr_pin,"mu"=curr_mu,"Tu"=curr_Tu,"vT"=curr_vT,"VT"=curr_VT,"invVT"=curr_invVT,"ElnTi"=curr_ElnTi,"ETi"=curr_ETi,"EmimiT"=curr_EmimiT)
  return(result)
}

compute_likelihood<-function(X,M,N,d,pin,Pi,mu,Tu,vT,VT,invVT,ElnTi,EmimiT,ETi,beta,v,V,lv,lV){
  val = (compute_LPD(M=M,N=N,d=d,pin=pin,ETi=ETi,ElnTi=ElnTi,X=X,mu=mu,EmimiT=EmimiT) 
         + compute_LPs(M=M,N=N,pin=pin,Pi=Pi)
         + compute_LPmu(M=M,d=d,beta=beta,EmimiT=EmimiT)
         + compute_LPT(M=M,v=lv[[1]],V=lV[[1]],d=d,ElnTi=ElnTi,ETi=ETi) 
         - compute_LQs(M=M,pin=pin,N=N) 
         - compute_LQmu(M=M,Tu=Tu,d=d) 
         - compute_LPQT(M=M,vT=vT,VT=VT,ETi=ETi,ElnTi=ElnTi,d=d))

#  print(c("Likelihood",val))
  return(val)
}

compute_LPD<-function(M,N,d,pin,ETi,ElnTi,X,mu,EmimiT){
  value=0
  for(i in 1:M){
    for(n in 1:N){
      value = value + pin[i,n]*((1/2)*ElnTi[[i]] 
                              - (d/2)*log(2*pi)
                              - (1/2)*matrix.trace(ETi[[i]]%*%(X[,n]%*%t(X[,n])
                                                               - X[,n]%*%t(mu[[i]])
                                                               - mu[[i]]%*%t(X[,n])
                                                               + EmimiT[[i]])))
    }
  }
#  print(c("LPD",value))
  return(value)
}

compute_LPs<-function(N,M,pin,Pi){
  value=0
  for(i in 1:M){
    for(n in 1:N){
      if(pin[i,n]>=1e-30 && Pi[i]>=1e-30){
        value = value + pin[i,n]*log(Pi[i])
      }
    }
  }
#  print(c("LPs",value))
  return(value)
}

compute_LPmu<-function(M,d,beta,EmimiT){
  value = (M*d/2)*(log(beta/(2*pi)))
  for(i in 1:M){
    value = value - (beta/2)*matrix.trace(EmimiT[[i]])
  }
#  print(c("LPmu",value))
  return(value)
}

compute_LPT<-function(M,v,V,d,ElnTi,ETi){
  val=0
  for(s in 1:d){
    val = val + lgamma((v+1-s)/2)
  }
  value = M*(-(v*d/2)*log(2) - (d*(d-1)/4)*log(pi) - val + (v/2)*mydet(V))
  Ti=matrix(0,ncol=d,nrow=d)
  LTi=0
  for(i in 1:M){
    Ti = Ti + ETi[[i]]
    LTi = LTi + ElnTi[[i]]
  }
  value = value + ((v-d-1)/2)*LTi - (1/2)*matrix.trace(V%*%Ti)
#  print(c("LPT",value))
  return(value)
}

compute_LQs<-function(M,pin,N){
  val = 0
  for(i in 1:M){
    for(n in 1:N){
      if(pin[i,n]>=1e-30){
        val = val + (pin[i,n]*log(pin[i,n]))
      }
    }
  }
#  print(c("LQs",val))
  return(val)
}

compute_LQmu<-function(M,Tu,d){
  val=0
  for(i in 1:M){
    tmp = (-(d/2)*(1 + log(2*pi)) + (1/2)*mydet(Tu[[i]]))
    val = val + tmp
  }
#  print(c("LQmu",val))
  return(val)
}

compute_LPQT<-function(M,vT,VT,ETi,ElnTi,d){
  val=0
  for(i in 1:M){
    value=0
    for(s in 1:d){
      value = value + lgamma((vT[[i]]+1-s)/2)
    }
#    tmp = -((vT[[i]]*d)/2)*log(2) - (d*(d-1)/4)*log(pi) - value + (vT[[i]]/2)*mydet(VT[[i]]) + ((vT[[i]]-d-1)/2)*ElnTi[[i]] - (1/2)*matrix.trace(VT[[i]]%*%ETi[[i]])
    tmp = -((vT[[i]]*d)/2)*log(2) - (d*(d-1)/4)*log(pi) - value + (vT[[i]]/2)*mydet(VT[[i]]) + ((vT[[i]]-d-1)/2)*ElnTi[[i]] - (1/2)*vT[[i]]
    val = val + tmp
  }
#  print(c("LPQT",val))
  return(val)
}

compute_Pi<-function(pin,N,M){
  Pi=matrix(rep(0,M),ncol=1)
  for(i in 1:M){
    Pi[i] = (1/N)*(sum(pin[i,]))
  }
#  print(c("PI",Pi))
  return(Pi)
}

compute_Tu<-function(beta,ETi,pin,M,d){
  Tu=list()
  for(i in 1:M){
    Tu[[i]] = make.positive.definite(diag(x=beta,nrow=d) + ETi[[i]]*(sum(pin[i,])),tol=1e-10)
#    print(c("Tu, i",i,Tu[[i]]))
  }
  return(Tu)
}

compute_ETi<-function(vT,invVT,M){
  ETi=list()
  for(i in 1:M){
    ETi[[i]] = vT[[i]]*(invVT[[i]])
#    print(c("ETi i",i,ETi[[i]]))
  }
  return (ETi)
}

compute_mu<-function(Tu,ETi,M,N,X,pin){
  mu=list()
  for(i in 1:M){
# initialize vector x
#    x=X[,1]*0
    x=vector("numeric",length=nrow(X))
    for(n in 1:N){
      x = x + X[,n]*pin[i,n]
    }
    
    if(kappa(Tu[[i]])<1000){
      mu[[i]] = (solve(Tu[[i]]))%*%(ETi[[i]]%*% x)
    }else{
#handle case when condition number kappa of Tu[[i]] is too large
      mu[[i]]=x/sum(pin[i,])
    }  
    
#    print(c("ETi",i,ETi[[i]]))
#    print(c("Tu[[i]]",i,Tu[[i]]))
#    print(c("solve(Tu[[i]])",i,solve(Tu[[i]])))
#    print(c("mu x", x))
#    print(c("condition number solve(Tu[[i]]) i", i, kappa(solve(Tu[[i]]))))
#    print(c("mu i",i))
#    print(mu[[i]])
  }
  return(mu)
}

compute_vT<-function(v,pin,M){
  vT=list()
  for(i in 1:M){
    vT[[i]] = v + sum(pin[i,]) 
#    print(c("vT i",i,vT[[i]]))
  }
  return(vT)
}

compute_pbar<-function(Pi,ElnTi,X,ETi,mu,EmimiT,M,N){
  pbar=matrix(rep(1,M*N),nrow=M)
  for(i in 1:M){
    x=ElnTi[[i]]/2 + log(Pi[i] + 1e-20)
    for(n in 1:N){
      y = matrix.trace(ETi[[i]] %*% (X[,n]%*%t(X[,n]) - mu[[i]]%*%t(X[,n]) - X[,n]%*%t(mu[[i]]) + EmimiT[[i]]))
      val = x - (1/2)*y 
      if(is.infinite(val)){
        stop("Infinite value")
        pbar[i,n] = 1e-30
      }else{
        pbar[i,n] = val
      }
    }
  }
  return(pbar)  
}

compute_VT<-function(V,X,M,N,pin,mu,EmimiT){
  VT=list()
  for(i in 1:M){
    z=0
    for(n in 1:N){
      z = z + X[,n]%*%(t(X[,n])*pin[i,n]) - X[,n]%*%(pin[i,n]*t(mu[[i]])) - mu[[i]]%*%(t(X[,n])*pin[i,n])   
    }
    VT[[i]] = make.positive.definite((V + z + EmimiT[[i]]*(sum(pin[i,]))),tol=1e-10)
  }
  return(VT)
}

compute_vT_VT_invVT<-function(v,pin,M,V,X,N,mu,EmimiT,method,alpha){
#  print(c("Compute_vT_VT_invVT_rho"))
#  print(c("cls",apply(pin,2,function(x) which.max(x))))
  
  vT=list()
  VT=list()
  invVT=list()  
  
  if (method==0){
    # No sparse
    vT=compute_vT(v=v,pin=pin,M=M)
    VT=compute_VT(V=V,X=X,M=M,N=N,pin=pin,mu=mu,EmimiT=EmimiT)
    for(i in 1:M){
      invVT[[i]]=solve(VT[[i]])
    }
  }else if (method==1){
    # pcor
    vT=compute_vT(v=v,pin=pin,M=M)
    VT=compute_VT(V=V,X=X,M=M,N=N,pin=pin,mu=mu,EmimiT=EmimiT)
    for(i in 1:M){
      omega=pcor_thresh(Sigma=VT[[i]],n=N,alpha=alpha)
      invVT[[i]]=omega$thresh.omega
    }
  }else {
    # method: 2 (glasso EBIC), 3 (glasso CV), 4 (tiger ebic), 5 (tiger cv)
    for(i in 1:M){
      temp=solve_eta_C_rho(method=method,d=nrow(X),N=N,v=v,V=V,s=pin[i,],X=X,mu=mu[[i]],EmimiT=EmimiT[[i]])
      vT[[i]]=temp$eta
      VT[[i]]=temp$C
      invVT[[i]]=temp$invC
    }
  }
    
  result=list("vT"=vT,"VT"=VT,"invVT"=invVT)
  return(result)
}

compute_EmimiT<-function(Tu, mu, M){
  EmimiT=list()
  for(i in 1:M){
    mumut=mu[[i]]%*%t(mu[[i]])
    if(kappa(Tu[[i]])<1000){
      EmimiT[[i]] = make.positive.definite(solve(Tu[[i]]) + mumut,tol=1e-5)
    }else{
      EmimiT[[i]] = mumut
    }  

#    print(c("EmimiT i",i,EmimiT[[i]]))
  }
  return(EmimiT)
}

compute_ElnTi<-function(VT,vT,d,M){
  ElnTi=list()
  for(i in 1:M){
    dg=0
    for(s in 1:d){
      dg = dg + digamma((vT[[i]] + 1 - s)/2)
    }
    ElnTi[[i]] = dg + d*log(2) - mydet(VT[[i]])
#    print(c("ElnTi i", i, ElnTi[[i]]))
  }
  return(ElnTi)
}

mydet<-function(Mat){
  return(as.numeric(determinant(Mat,logarithm=TRUE)$modulus))
}

compute_pin<-function(P){
  M=nrow(P)
  N=ncol(P)
  
  pin=matrix(rep(0,M*N),nrow=M)
  for(n in 1:N){
    x=P[,n]
    for(i in 1:M){
      pin[i,n] = 1/(sum(exp(x-P[i,n])))
    }
  }
  
  for(i in 1:M){
#    print(c("pin[i,] i",i,pin[i,]))
  }
  
  return(pin)
}

pcor_thresh<-function(Sigma,n,alpha=0.05){
  d=dim(Sigma)[1]
  s=d-1

#  print(c("condition number Sigma",kappa(Sigma)))
#  print(c("Printing sigma"))
#  print(Sigma)
  
#  while(kappa(Sigma)>10000){
#    Sigma=Sigma + diag(10000,nrow=d)
#    print(c("Updated condition number Sigma",kappa(Sigma)))
#  }
  
  #MLE estimate for Omega = inverse(Sigma)
  Omega=solve(Sigma)
  pc=mypcor(Omega) 
  
#  print(c("condition number sigma and omega",kappa(Sigma),kappa(Omega)))
  
  thresh_parcors=ci_par_cor(par_cors=pc, alpha=alpha, n=n, s=s) 
  
  thresh_omega= (thresh_parcors$mat) * Omega
 
  #ensuring thresh_omega is well conditioned
  while(kappa(thresh_omega)>kappa(Omega)){
    thresh_omega=thresh_omega + diag(1,nrow=d)
  } 
  
#  print(c("condition number thresh_omega",kappa(thresh_omega)))
  
  return(list("thresh.omega"=thresh_omega))
}

#From Pcor method Williams and Rast 2018
ci_par_cor <- function(par_cors,alpha,n,s){
  # n: sample size 
  # s: d - 1 (controlled for)
  # alpha: confidence level 
  # par_cors: partial correlations
  mat <- matrix(0,nrow =  s + 1, ncol = s + 1)
  CI_ls <- list()
  par_cor <- par_cors[upper.tri(par_cors)]
  cov <- list()
  for(i in 1:length(par_cor)){
    # critical value
    z_crit <- qnorm(1 - alpha/2)
    # standard error
    se <- sqrt(1/(max(n - s - 3, 1)))
#    se <- sqrt(1/(n - s - 3))
    # z transformation
    z <- log((1 + par_cor[i])/(1 - par_cor[i]))/2
    # z lower bound
    Z_L <- z - z_crit * se
    # Z upper bound 
    
    Z_U <- z + z_crit * se
    rho_L <- (exp(2*Z_L) - 1)/(exp(2*Z_L) + 1)
    rho_U <- (exp(2*Z_U) - 1)/(exp(2*Z_U) + 1)
    CI <- c(rho_L, rho_U)
    CI_ls[[i]] <- CI
    cov[[i]] <- ifelse(CI[1] < 0 & CI[2] > 0, 0, 1)
  }
  
  diag(mat) <- 1
  mat[upper.tri(mat)] <- unlist(cov)
  mat <- as.matrix(forceSymmetric(mat) ) 
  return(list("mat"=mat))
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

solve_eta_C_rho<-function(method,d,N,v,V,s,X,mu,EmimiT){
  val=v+sum(s)
  eta=val
  res=compute_prec_rho(method=method,d=d,N=N,v=v,V=V,s=s,eta=eta,X=X,mu=mu,EmimiT=EmimiT)
  tomega=res$icov
  Mat=res$Mat
  if(!isSymmetric(tomega)){
    tomega = (tomega + t(tomega))/2
  }
  invC=make.positive.definite(tomega,tol=1e-6)
  C=solve(invC)
  rho = 0.5*(v+sum(s))*(res$rho)
  score=penalized_elbo(invC=invC,Mat=Mat,eta=eta,d=d,val=val,rho=rho)
  
  if((0.5*(d-matrix.trace(Mat%*%invC)) + 0.25*sum(trigamma(seq(from=val,to=val+1-d,by=-1)/2))*(val-d))<0){
    new_eta=val
    #    print(invC)
    #    print("Mat")
    #    print(Mat)
    #    print("Trace MinvC")
    #    print(matrix.trace(M%*%invC))
    #    print("TraceSinvC")
    #    print(matrix.trace(S%*%invC))
    #    stop("stopping")
  }else{
    new_eta=mysolve_eta(d=d,c=val,t=matrix.trace(Mat%*%invC))
  }
  
  res=compute_prec_rho(method=method,d=d,N=N,v=v,V=V,s=s,eta=new_eta,X=X,mu=mu,EmimiT=EmimiT)
  Mat = res$Mat
  tomega=res$icov
  if(!isSymmetric(tomega)){
    tomega = (tomega + t(tomega))/2
  }
  new_invC=make.positive.definite(tomega,tol=1e-6)
  new_C=solve(new_invC)
  new_rho = 0.5*(v+sum(s))*(res$rho)
  new_score=penalized_elbo(invC=new_invC,Mat=Mat,eta=new_eta,d=d,val=val,rho=new_rho)
  
  cnt=0
  while(new_score>score){
    cnt=cnt+1
    #    print(c("cnt ",cnt))
    
    C=new_C
    invC=new_invC
    eta=new_eta
    rho=new_rho
    score=new_score
    
    if((0.5*(d-matrix.trace(Mat%*%invC)) + 0.25*sum(trigamma(seq(from=val,to=val+1-d,by=-1)/2))*(val-d))<0){
      new_eta=eta
    }else{
      new_eta=mysolve_eta(d=d,c=val,t=matrix.trace(Mat%*%invC))
    }
    
    res=compute_prec_rho(method=method,d=d,N=N,v=v,V=V,s=s,eta=new_eta,X=X,mu=mu,EmimiT=EmimiT)
    Mat = res$Mat
    tomega=res$icov
    if(!isSymmetric(tomega)){
      tomega = (tomega + t(tomega))/2
    }
    new_invC=make.positive.definite(tomega,tol=1e-6)
    new_C=solve(new_invC)
    new_rho = 0.5*(v+sum(s))*(res$rho)
    new_score=penalized_elbo(invC=new_invC,Mat=Mat,eta=new_eta,d=d,val=val,rho=new_rho)
  }
  
  #  print("DONE solve_eta")
  result=list("eta"=eta,"C"=C,"invC"=invC,"rho"=rho)
  return(result)
}

penalized_elbo<-function(invC,Mat,eta,d,val,rho){
  score = 0.5*val*mydet(invC) - 0.5*eta*matrix.trace(Mat%*%invC) - rho*sum(abs(invC))  + 0.5*(val-eta)*sum(digamma(seq(from=eta,to=eta+1-d,by=-1)/2)) + 0.5*eta*d + sum(lgamma(seq(from=eta,to=eta+1-d,by=-1)/2))
  return(score)
}

compute_prec_rho<-function(method,d,N,v,V,s,eta,X,mu,EmimiT){
  val=v+sum(s)
  nrho=5 # can be parameterized
  z=0
  for(i in 1:N){
    z = z + X[,i]%*%(t(X[,i])*s[i]) - X[,i]%*%(s[i]*t(mu)) - mu%*%(t(X[,i])*s[i])   
  }
  #S is empirical covariance matrix
  S = make.positive.definite((V + z + EmimiT*(sum(s))),tol=1e-10)
  S = (eta/val)*S
  
  Sm=S
  diag(Sm)=0
  rho.max=max(max(abs(Sm)) + 1e-6, pi*sqrt(log(d)/N))
  rho.min.ratio=0.01 # can be turned into parameter
  rho.min=rho.min.ratio * rho.max
  lrho=10^seq(log10(rho.min), log10(rho.max), length = nrho)
  lrho=sort(lrho,decreasing=TRUE)
  
  if(method==2 || method==4){
    # ebic
    if(method==2){
      # glasso
      res = huge(x=S,lambda=lrho,method="glasso",cov.output=FALSE,verbose=FALSE)
    }else{
      # tiger
      res = huge(x=S,lambda=lrho,method="tiger",cov.output=FALSE,verbose=FALSE)
    }
    ebic = vector(mode="numeric",length=nrho)
    thr=1e-6 #threshold for zero; can be turned into parameter
    g=0.5 # gamma value in EBIC
    for(i in 1:nrho){
      k=0
      k=sum(abs(res$icov[[i]])>thr)
      vll= (N/2)*(-d*log(2*pi) + compute_loglike_glasso_tiger(Omega=res$icov[[i]], Sigma=S, rho=lrho[i]))
      vebic= k*log(N) - 2*vll + 4*g*k*log(d)
      ebic[i] = vebic
    }
    ind=which.min(ebic)
    return(list("icov" = res$icov[[ind]], "Mat" = S, "rho"= lrho[ind]))
  }else{
    # cv
    # three fold stratified cross validation
    fld=3
    #stratified
    in_comp=sample(which(s>=0.4))
    out_comp=sample(which(s<0.4))
#    print(c("in",length(in_comp),"out",length(out_comp)))
    fold=vector(mode="numeric",length=N)
    for(ij in 1:length(in_comp)){
      mij = (ij %% fld) + 1
      fold[in_comp[ij]] = mij
    }
    for(ij in 1:length(out_comp)){
      mij = (ij %% fld) + 1
      fold[out_comp[ij]] = mij
    }
   
#    print(fold)
    
    llcv=vector(mode="numeric",length=nrho) # log likelihood for cv
    for(kind in 1:fld){
      #Strain
      z=0
      for(i in which(fold!=kind)){
        z = z + X[,i]%*%(t(X[,i])*s[i]) - X[,i]%*%(s[i]*t(mu)) - mu%*%(t(X[,i])*s[i])   
      }
      #Strain is empirical covariance matrix
      Strain = make.positive.definite((V + z + EmimiT*(sum(s))),tol=1e-10)
      Strain = (eta/val)*Strain
      
      #Stest
      z=0
      for(i in which(fold==kind)){
        z = z + X[,i]%*%(t(X[,i])*s[i]) - X[,i]%*%(s[i]*t(mu)) - mu%*%(t(X[,i])*s[i])   
      }
      #Stest is empirical covariance matrix
      Stest = make.positive.definite((V + z + EmimiT*(sum(s))),tol=1e-10)
      Stest = (eta/val)*Stest
      
      if(method==3){
        # glasso
        res = huge(x=Strain,lambda=lrho,method="glasso",cov.output=FALSE,verbose=FALSE)
      }else{
        # tiger
        res = huge(x=Strain,lambda=lrho,method="tiger",cov.output=FALSE,verbose=FALSE)
      }
      
      for(i in 1:nrho){
        llcv[i] = llcv[i] + compute_loglike_glasso_tiger(Omega=res$icov[[i]], Sigma=Stest, rho=lrho[i])
      }
    }
    ind=which.max(llcv)
    if(method==3){
      # glasso
      res = huge(x=S,lambda=lrho,method="glasso",cov.output=FALSE,verbose=FALSE)
    }else{
      # tiger
      res = huge(x=S,lambda=lrho,method="tiger",cov.output=FALSE,verbose=FALSE)
    }
    return(list("icov" = res$icov[[ind]], "Mat" = S, "rho" = lrho[ind]))
  }
}


compute_loglike_glasso_tiger<-function(Omega,Sigma,rho){
  #  val = mydet(Omega) - matrix.trace(Sigma%*%Omega) - rho*(sum(abs(Omega)))
  val = mydet(Omega) - matrix.trace(Sigma%*%Omega)
  return(val)
}

# Newton-Raphson method from numerical recipes in C after modifications to fit the function we have
# solving for eta in Wishart distribution
mysolve_eta<-function(d,c,t,xacc=1e-3){
  #c is (sum_si + v)
  #t is Tr(MC^{-1})
  
  #  print(c("d=",d,"c=",c,"t=",t))  
  
  if(c<d){
    stop("Error: c cannot be less than d")
  }
  
  xl=d
  xh=c
  fl= 0.5*(d-t) + 0.25*sum(trigamma(seq(from=xl,to=xl+1-d,by=-1)/2))*(c-xl)
  fh= 0.5*(d-t) + 0.25*sum(trigamma(seq(from=xh,to=xh+1-d,by=-1)/2))*(c-xh)
  
  cnt=1
  while(fh>=0 && cnt<=12){
    cnt=cnt+1
    xh=2*xh
    fh= 0.5*(d-t) + 0.25*sum(trigamma(seq(from=xh,to=xh+1-d,by=-1)/2))*(c-xh)
  }
  
  if(cnt>12){
    print(c("c= ",c," xh= ",xh))
    print(c("d= ",d, " t" = t))
    stop(c("Error: Check c as xh is pretty high"))
  }
  
  if(fl<0){
    print(c("xl = ", xl, "fl = ",fl, ", xh = ",xh,", fh = ",fh))
    warning(c("Error in root finding"," d= ",d," c= ",c, " t= ",t))
    return(xl)
  }else{
    #    print(c("xl = ", xl, "fl = ",fl, ", xh = ",xh,", fh = ",fh))
  }
  
  
  rts=0.5*(xl+xh)
  dxold=abs(xh-xl)
  dx=dxold
  f= 0.5*(d-t) + 0.25*sum(trigamma(seq(from=rts,to=rts+1-d,by=-1)/2))*(c-rts)
  df= -0.5*sum(trigamma(seq(from=rts,to=rts+1-d,by=-1)/2)) + 0.25*sum(psigamma(seq(from=rts,to=rts+1-d,by=-1)/2, deriv=2))*(c-rts)
  
  
  MAXIT=100
  j=0
  while(j<=MAXIT){
    #    print(rts)
    j=j+1
    if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0) || (abs(2*f)>abs(dxold*df))){
      #      print("HERE1")
      dxold=dx
      dx=0.5*(xh-xl)
      rts=xl+dx
      if(xl==rts){
        return (rts)
      }
    }else{
      #      print("HERE2")
      dxold=dx
      dx=f/df
      temp=rts
      rts=rts - dx
      if(temp==rts){
        return (rts)
      }
    }
    if(abs(dx) < xacc) {
      return (rts)
    }
    #    print(c("rts=",rts, "dx = ",dx))
    f= 0.5*(d-t) + 0.25*sum(trigamma(seq(from=rts,to=rts+1-d,by=-1)/2))*(c-rts)
    df= -0.5*sum(trigamma(seq(from=rts,to=rts+1-d,by=-1)/2)) + 0.25*sum(psigamma(seq(from=rts,to=rts+1-d,by=-1)/2, deriv=2))*(c-rts)
    
    #set the new left and right boundaries
    if(f<0){
      xh=rts
    }else{
      xl=rts
    }
  }
  return (xl)
}


