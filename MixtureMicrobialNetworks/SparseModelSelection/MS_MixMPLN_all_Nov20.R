library(mvtnorm)
library(Matrix)
library(matrixcalc)
library(corpcor)
library(glasso)
#library(cvTools)
library(huge)
library(mclust)
#library(bayesm)


#######################################################################
# Variational approximation algorithm for model selection for MixMPLN

# d: number of taxa (precision matrices are d x d)
# n: number of samples (each sample is a d-dimensional count vector)
# K: number of MixMPLN components
# X_dxn is the data matrix (the taxa-sample matrix)
# Pi is K dimensional vector of mixing coefficients

# EOmega=list() is list of expected values of the K MPLN precision matrices, each of size dxd
# EMu=list() is list of expected values of the K MPLN mean vectors, each of size d
# ELambda=list() is list of expected values of the K MPLN lambda matrices, each of size dxn
# The expected values of the corresponding RVs wrt their respective factor functions

# beta is constant in Gaussian prior for each mean vector Mu; Mu~Norm(0,beta*I)
# v, V are part of Wishart prior for each precision matrix Omega; Omega~Wishart(v,V)

#Es_li = pli Expected value of s_li
#variables i=1:n, j=1:d, l=1:K

#parameters of the factor functions
#q(Si): pli
#q(Mu_l): EMu, Tu
#q(Omega_l): wO, WOmega
#q(Lambda_l:i): a_lji, b_lji where a is mean, b is variance, c=1/b; ELambda will store the a_lji values and list B will
# the l dxn variance values; each matrix in B is a diagonal matrix

#ELambda: list with K elements, each element is a dxn matrix
#(TOO LARGE TO CARRY AROUND) ELlLlT: list with K elements, where each element is a list with n elements, each of which is a dxd matrix
#EexpLambda: list with K elements, each element is a dxn matrix
#B: list with K elements, each element is a dxn matrix (These are matrices of variances as opposed to matrices of precision values)
#UA: same as ELambda
#UB: same as B
#EMu: list with K elements, each element is a dx1 matrix
#EOmega: list with K elements, each element is a dxd precision matrix
#ElnOmega: list with K elements, each element is a number
#EMlMlT: list with K elements, each element is a dxd matrix
#wO: list with K elements, each element is Wishart degree freedom scale value
#WOmega: list with K elements, each element is dxd matrix Wishart scale matrix
#Tu: list with K elements each element is a dxd precision matrix
#pli: Kxn matrix
#Pi: Kx1 matrix

# Can run in sparse or non-sparse mode
# Sparse mode: pcor, GLASSO (EBIC, CV), TIGER (EBIC, CV)
# Method: 0 (non-sparse), 1 (pcor), 2 (glasso EBIC), 3 (glasso CV), 4 (tiger ebic), 5 (tiger cv)


best_model<-function(X,F,K,beta=1e-6,v,V,niter=30,seed=1001,start=10,method=0,alpha=0.05){
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
  n=ncol(X)
  
  set.seed(seed)
  
  print(c("Method", method))
  start.time=as.numeric(format(Sys.time(),"%s"))
  
  #run k-means
  print(c("Running K-means"))
  member=kmeans(x=t(X),centers=K,iter.max=100)$cluster
  print(c("Table member",table(member)))
  result=model_select(X=X,K=K,v=v,V=V,beta=beta,niter=niter,member=member,method=method,alpha=alpha)
  print(c("K-means likelihood",result$likelihood))
  print(c("K-means Pi",result$Pi))
  print(c("Iterations",result$iterations))
  gc()
  
  #run random starting points using seed  
  for(i in 1:start){
    writeLines("\n")
    print(c("Random start point iteration",i))
    member=sample.int(K,size=n,replace=TRUE)
    while(length(table(member))!=K){
      member=sample.int(K,size=n,replace=TRUE)
    }
    print(c("Table member",table(member)))
    result_next=model_select(X=X,K=K,v=v,V=V,beta=beta,niter=niter,member=member,method=method,alpha=alpha)
    print(c("Iteration Likelihood",result_next$likelihood))
    print(c("Iteration Pi",result_next$Pi))
    print(c("Iterations",result_next$iterations))
    
    if(result_next$likelihood>result$likelihood){
      result=result_next
      print(c("Updated"))
    }
    gc()
  }
  
  end.time=as.numeric(format(Sys.time(),"%s"))
  
  result$X=X
  result$K=K
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
  
  gc()
  return(result)
}

model_select<-function(X,K,beta,v,V,niter,member,method,alpha){
  X<-as.matrix(X)
  d=nrow(X)
  n=ncol(X)
  
  if(K <= 0){
    stop("Error: provide positive integer for K")
  }
  
  if(n < 4*K){
    stop("Error: too few data points")
  }
  
  if(length(member)!=n){
    print(c("Error", n, length(member)))
    stop("Error: n != length(member)")
  }
  
  EMu=list()
  Tu=list()
  wO=list()
  WOmega=list()
  invWOmega=list()
  EOmega=list()
  EMlMlT=list()
  ElnOmega=list()
  ELambda=list()
  EexpLambda=list()
  B=list()

    
  pli=matrix(rep(1e-5,K*n),nrow=K)
  for(l in 1:K){
    ind=which(member==l)
    pli[l,ind]=1
  }

  for(i in 1:n){
    csum=sum(pli[,i])
    pli[,i]=pli[,i]/csum
  }
  
  Pi=compute_Pi(pli=pli,n=n,K=K)
  
#  print(c("Pi",Pi))
  
  for(l in 1:K){
    ind=which(member==l)
#   print(c("ind l",l,ind))
    if(length(ind)==1){
      tmu=log(X[,ind]+1)
    }else{
      tmu=rowMeans(log(X[,ind]+1))
    }

    B[[l]]=matrix(rep(1e-5,d*n),ncol=n)
    ELambda[[l]]=matrix(rep(1,d*n),ncol=n)*tmu
    EexpLambda[[l]] = ELambda[[l]] + 0.5*B[[l]]
    EMu[[l]]=tmu
    EMlMlT[[l]]=tmu%*%t(tmu)
  }

  temp=compute_wO_WOmega_invWOmega(v=v,V=V,K=K,n=n,pli=pli,EMu=EMu,EMlMlT=EMlMlT,ELambda=ELambda,B=B,d=d,method=method,alpha=alpha)
  wO=temp$wO
  WOmega=temp$WOmega
  invWOmega=temp$invWOmega
  EOmega=compute_EOmega(n=n,wO=wO,invWOmega=invWOmega,K=K)
  ElnOmega=compute_ElnOmega(WOmega=WOmega,wO=wO,d=d,K=K)
  
  pbar=compute_pbar(Pi=Pi,ElnOmega=ElnOmega,X=X,EOmega=EOmega,EMu=EMu,EMlMlT=EMlMlT,ELambda=ELambda,B=B,EexpLambda=EexpLambda,K=K,n=n,d=d)
  pli=compute_pli(pbar)
  Pi=compute_Pi(pli=pli,n=n,K=K)
  print ("All Initialized")
  
#  likelihood=compute_likelihood(X=X,K=K,n=n,d=d,pli=pli,Pi=Pi,EMu=EMu,EOmega=EOmega,wO=wO,WOmega=WOmega,
#                                ElnOmega=ElnOmega,EMlMlT=EMlMlT,Tu=Tu,ELambda=ELambda,EexpLambda=EexpLambda,
#                                B=B,beta=beta,v=v,V=V)
#  print(c("Pi",Pi))
  
#  print_all(K=K,n=n,d=d,pli=pli,Pi=Pi,EMu=EMu,EOmega=EOmega,wO=wO,WOmega=WOmega,invWOmega=invWOmega,
#            ElnOmega=ElnOmega,EMlMlT=EMlMlT,Tu=Tu,ELambda=ELambda,EexpLambda=EexpLambda,
#            B=B,likelihood=0)

  Tu=compute_Tu(beta=beta,EOmega=EOmega,pli=pli,K=K,d=d)  
  EMu=compute_EMu(Tu=Tu,EOmega=EOmega,K=K,n=n,ELambda=ELambda,pli=pli)
  EMlMlT=compute_EMlMlT(Tu=Tu, EMu=EMu, K=K)
  Lambdavals=compute_LambdaValues(ELambda=ELambda,EexpLambda=EexpLambda,B=B,X=X,EOmega=EOmega,EMu=EMu,pli=pli,d=d,K=K,n=n)
  ELambda=Lambdavals$ELambda
  EexpLambda=Lambdavals$EexpLambda
  B=Lambdavals$B
  temp=compute_wO_WOmega_invWOmega(v=v,V=V,K=K,n=n,pli=pli,EMu=EMu,EMlMlT=EMlMlT,ELambda=ELambda,B=B,d=d,method=method,alpha=alpha)
  wO=temp$wO
  WOmega=temp$WOmega
  invWOmega=temp$invWOmega
  EOmega=compute_EOmega(n=n,wO=wO,invWOmega=invWOmega,K=K)
  ElnOmega=compute_ElnOmega(WOmega=WOmega,wO=wO,d=d,K=K)
  pbar=compute_pbar(Pi=Pi,ElnOmega=ElnOmega,X=X,EOmega=EOmega,EMu=EMu,EMlMlT=EMlMlT,ELambda=ELambda,B=B,EexpLambda=EexpLambda,K=K,n=n,d=d)
  pli=compute_pli(pbar)
  Pi=compute_Pi(pli=pli,n=n,K=K)
  
  likelihood=compute_likelihood(X=X,K=K,n=n,d=d,pli=pli,Pi=Pi,EMu=EMu,EOmega=EOmega,wO=wO,WOmega=WOmega,
                                ElnOmega=ElnOmega,EMlMlT=EMlMlT,Tu=Tu,ELambda=ELambda,EexpLambda=EexpLambda,
                                B=B,beta=beta,v=v,V=V)

# print(c("Next"))
# print_all(K=K,n=n,d=d,pli=pli,Pi=Pi,EMu=EMu,EOmega=EOmega,wO=wO,WOmega=WOmega,invWOmega=invWOmega,
#           ElnOmega=ElnOmega,EMlMlT=EMlMlT,Tu=Tu,ELambda=ELambda,EexpLambda=EexpLambda,
#            B=B,likelihood=likelihood)
  
  curr_pli=pli
  curr_Pi=Pi
  curr_ELambda=ELambda
  curr_EexpLambda=EexpLambda
  curr_B=B
  curr_EOmega=EOmega
  curr_ElnOmega=ElnOmega
  curr_EMu=EMu
  curr_EMlMlT=EMlMlT
  curr_Tu=Tu
  curr_wO=wO
  curr_WOmega=WOmega
  curr_invWOmega=invWOmega
  curr_likelihood=likelihood
  
  loop=TRUE
  iter=0
  while(loop && iter<niter){
    iter = iter + 1
    
    #E-step
    Tu=compute_Tu(beta=beta,EOmega=EOmega,pli=pli,K=K,d=d)
    EMu=compute_EMu(Tu=Tu,EOmega=EOmega,K=K,n=n,ELambda=ELambda,pli=pli)
    EMlMlT=compute_EMlMlT(Tu=Tu, EMu=EMu, K=K)
    Lambdavals=compute_LambdaValues(ELambda=ELambda,EexpLambda=EexpLambda,B=B,X=X,EOmega=EOmega,EMu=EMu,pli=pli,d=d,K=K,n=n)
    ELambda=Lambdavals$ELambda
    EexpLambda=Lambdavals$EexpLambda
    B=Lambdavals$B
    temp=compute_wO_WOmega_invWOmega(v=v,V=V,K=K,n=n,pli=pli,EMu=EMu,EMlMlT=EMlMlT,ELambda=ELambda,B=B,d=d,method=method,alpha=alpha)
    wO=temp$wO
    WOmega=temp$WOmega
    invWOmega=temp$invWOmega
    EOmega=compute_EOmega(n=n,wO=wO,invWOmega=invWOmega,K=K)
    ElnOmega=compute_ElnOmega(WOmega=WOmega,wO=wO,d=d,K=K)
    pbar=compute_pbar(Pi=Pi,ElnOmega=ElnOmega,X=X,EOmega=EOmega,EMu=EMu,EMlMlT=EMlMlT,ELambda=ELambda,B=B,EexpLambda=EexpLambda,K=K,n=n,d=d)
    pli=compute_pli(pbar)

    #M-step
    Pi=compute_Pi(pli=pli,n=n,K=K)
    
    likelihood=compute_likelihood(X=X,K=K,n=n,d=d,pli=pli,Pi=Pi,EMu=EMu,EOmega=EOmega,wO=wO,WOmega=WOmega,
                                  ElnOmega=ElnOmega,EMlMlT=EMlMlT,Tu=Tu,ELambda=ELambda,EexpLambda=EexpLambda,
                                  B=B,beta=beta,v=v,V=V)
  
#    print(c("Loop"))
#    print_all(K=K,n=n,d=d,pli=pli,Pi=Pi,EMu=EMu,EOmega=EOmega,wO=wO,WOmega=WOmega,invWOmega=invWOmega,
#              ElnOmega=ElnOmega,EMlMlT=EMlMlT,Tu=Tu,ELambda=ELambda,EexpLambda=EexpLambda,
#              B=B,likelihood=likelihood)
       
    if(likelihood > curr_likelihood){
      curr_pli=pli
      curr_Pi=Pi
      curr_ELambda=ELambda
      curr_EexpLambda=EexpLambda
      curr_B=B
      curr_EOmega=EOmega
      curr_ElnOmega=ElnOmega
      curr_EMu=EMu
      curr_EMlMlT=EMlMlT
      curr_Tu=Tu
      curr_wO=wO
      curr_WOmega=WOmega
      curr_invWOmega=invWOmega
      curr_likelihood=likelihood
    }else{
      loop=FALSE
    }
  }
  
  #Results in curr_ variables
  result=list("likelihood"=curr_likelihood,"iterations"=iter,"ELambda"=curr_ELambda,"EexpLambda"=curr_EexpLambda,
              "B"=curr_B,"wO"=curr_wO,"WOmega"=curr_WOmega,"invWOmega"=curr_invWOmega,"EOmega"=curr_EOmega,"ElnOmega"=curr_ElnOmega,
              "Tu"=curr_Tu,"EMu"=curr_EMu,"EMlMlT"=curr_EMlMlT,"pli"=curr_pli,"Pi"=curr_Pi)
  return(result)
}

print_all<-function(K,n,d,pli,Pi,EMu,EOmega,wO,WOmega,invWOmega,ElnOmega,EMlMlT,Tu,ELambda,EexpLambda,B,likelihood){
  print(c("All temp"))
  print(c("Likelihood",likelihood))
  print(c("Pi",Pi))
  for(l in 1:K){
    print(c("Component l",l))
    print(c("pli",l,pli[l,]))
    print(c("EOmega ElnOmega", ElnOmega[[l]]))
    print(EOmega[[l]])
    print("EMu")
    print(EMu[[l]])
    print("wO")
    print(wO[[l]])
    print("WOmega")
    print(WOmega[[l]])
    print("invWOmega")
    print(invWOmega[[l]])
    print(c("ELambda, B,EexpLambda"))
    for(i in 1:n){
      print(c("ELambda",l,i,ELambda[[l]][,i]))
      print(c("Lamsigma",l,i,B[[l]][,i]))
      print(c("EexpLamda",l,i,EexpLambda[[l]][,i]))
   }
  }
}

compute_likelihood<-function(X,K,n,d,pli,Pi,EMu,EOmega,wO,WOmega,ElnOmega,EMlMlT,Tu,ELambda,EexpLambda,B,beta,v,V){
  val = (compute_LPData(K=K,n=n,d=d,pli=pli,EOmega=EOmega,ElnOmega=ElnOmega,X=X,EMu=EMu,
                        EMlMlT=EMlMlT,ELambda=ELambda,EexpLambda=EexpLambda,B=B,Pi=Pi)
         + compute_LPOmega(K=K,v=v,V=V,d=d,ElnOmega=ElnOmega,EOmega=EOmega)
         + compute_LPMu(K=K,d=d,beta=beta,EMlMlT=EMlMlT)
         - compute_LQs(K=K,pli=pli,n=n)
         - compute_LQMu(K=K,Tu=Tu,d=d)
         - compute_LQOmega(K=K,wO=wO,WOmega=WOmega,EOmega=EOmega,ElnOmega=ElnOmega,d=d)
         - compute_LQLambda(K=K,n=n,d=d,B=B))

  return(val)
}


compute_LPData<-function(K,n,d,pli,EOmega,ElnOmega,X,EMu,EMlMlT,ELambda,EexpLambda,B,Pi){
  val=0
  for(l in 1:K){
    ELl=ELambda[[l]]
    EexpLambdal=EexpLambda[[l]]
    Bl=B[[l]]
    for(i in 1:n){
      if(pli[l,i]<1e-10){
        x=0
      }else{
        z=0
        ELlLlTi = make.positive.definite(diag(Bl[,i], nrow=d) + ELl[,i]%*%t(ELl[,i]),tol=1e-5)
        for (j in 1:d){
          z = z + (-exp(min(50,EexpLambdal[j,i])) - lfactorial(X[j,i]) +ELl[j,i]*X[j,i]) 
        }
        y = matrix.trace(EOmega[[l]] %*% (ELlLlTi - EMu[[l]]%*%t(ELl[,i]) - ELl[,i]%*%t(EMu[[l]]) + EMlMlT[[l]]))
        x = pli[l,i]*(0.5*ElnOmega[[l]] -(d/2)*log(2*pi) - (1/2)*y +z) 
#        print(c("LPData l,i",l,i,pli[l,i],x,z))
      }
      val = val + x
    }
  }
  val_LPs = compute_LPs(n=n,K=K,pli=pli,Pi=Pi)
  rval = val + val_LPs
#  print(c("LPData",rval))
  return(rval)  
}

compute_LQLambda<-function(K,n,d,B){
  val=0
  for(l in 1:K){
    Bl=B[[l]]
    for(i in 1:n){
      val = val + 0.5*(sum(log(Bl[,i]))) -(d/2)*(1+log(2*pi))
    }
  }
#  print(c("LQLambda", val))
  return(val)
}

compute_LQOmega<-function(K,wO,WOmega,EOmega,ElnOmega,d){
  val=0
  for(l in 1:K){
    value=0
    for(s in 1:d){
      value = value + lgamma((wO[[l]]+1-s)/2)
    }
    tmp = - ((wO[[l]]*d)/2)*log(2) - (d*(d-1)/4)*log(pi) - value + (wO[[l]]/2)*mydet(WOmega[[l]]) + ((wO[[l]]-d-1)/2)*ElnOmega[[l]] - (1/2)*matrix.trace(WOmega[[l]]%*%EOmega[[l]])
#    tmp = - ((wO[[l]]*d)/2)*log(2) - (d*(d-1)/4)*log(pi) - value + (wO[[l]]/2)*mydet(WOmega[[l]]) + ((wO[[l]]-d-1)/2)*ElnOmega[[l]] - (1/2)*wO[[l]]   
    val = val + tmp
  }
#  print(c("LQOmega", val))
  return(val)
}

compute_LPOmega<-function(K,v,V,d,ElnOmega,EOmega){
  val=0
  for(s in 1:d){
    val = val + lgamma((v+1-s)/2)
  }
  value = K*(-(v*d/2)*log(2) - (d*(d-1)/4)*log(pi) - val + (v/2)*mydet(V))
  Tl=matrix(0,ncol=d,nrow=d)
  LTl=0
  for(l in 1:K){
    Tl = Tl + EOmega[[l]]
    LTl = LTl + ElnOmega[[l]]
  }
  value = value + ((v-d-1)/2)*LTl - (1/2)*matrix.trace(V%*%Tl)
#  print(c("LPOmega",value))
  return(value)
}

compute_LQMu<-function(K,Tu,d){
  val=0
  for(l in 1:K){
    tmp = (-(d/2)*(1 + log(2*pi)) + (1/2)*mydet(Tu[[l]]))
    val = val + tmp
  }
#  print(c("LQmu",val))
  return(val)
}

compute_LPMu<-function(K,d,beta,EMlMlT){
  value = (K*d/2)*(log(beta/(2*pi)))
  for(l in 1:K){
    value = value - (beta/2)*matrix.trace(EMlMlT[[l]])
  }
#  print(c("LPMu",value))
  return(value)
}

compute_LQs<-function(K,pli,n){
  val = 0
  for(l in 1:K){
    for(i in 1:n){
      val = val + (pli[l,i]*log(pli[l,i]+(1e-20)))
    }
  }
#  print(c("LQs",val))
  return(val)
}

compute_LPs<-function(n,K,pli,Pi){
#  print(c("LPs Pi", Pi))
# print(c("LPs pli"))
#  print(t(pli))
  value=0
  for(l in 1:K){
    for(i in 1:n){
      if(pli[l,i]>=1e-30 && Pi[l]>=1e-30){
        value = value + pli[l,i]*log(Pi[l]+(1e-20))
      }
    }
  }
#  print(c("LPs",value))
  return(value)
}

compute_LambdaValues<-function(ELambda,EexpLambda,B,X,EOmega,EMu,pli,d,K,n){
  UA=list()
  UB=list()
  UEexpLambda=list()
  for(l in 1:K){
    Atmp=ELambda[[l]]
    Btmp=B[[l]]
    EexpLambdatmp=EexpLambda[[l]]
    EOmegal=EOmega[[l]]
    EMul=EMu[[l]]
    for(i in 1:n){
      #gradient search to find optimal values of Atmp[,i] and Btmp[,i]
      val = compute_valLambdafn(d=d,Xi=X[,i],Eei=EexpLambdatmp[,i],Ai=Atmp[,i],Bi=Btmp[,i],plili=pli[l,i],EOmegal=EOmegal,EMul=EMul)
      loop=TRUE
      iter=0
      while(loop){
        iter=iter+1
        
        curr_Ai=Atmp[,i]
        curr_Bi=Btmp[,i]

#        print(c("curr_Ai",i, curr_Ai))
#        print(c("curr_Bi",i, curr_Bi))
#        print(c("Eomegal l",l))
#        print(EOmegal)
                
        for(j in 1:d){
          ta = exp(-Btmp[j,i]/2)
          ba = EOmegal[j,j]*ta
          ca = (X[j,i] + sum(EOmegal[j,]*EMul) - sum(EOmegal[j,]*Atmp[,i]) + EOmegal[j,j]*Atmp[j,i])*ta
#          print(c("i d ba ca",i,d,ba,ca))
#          print(c("pli l,i",l,i,pli[l,i]))
#          print(c("EMul",EMu))
          curr_Ai[j] = mysolve_a(b=ba,c=ca)
#          print(c("currAi",curr_Ai[j]))

          if(Atmp[j,i] < -50){
            if(pli[l,i] >= 1e-10){
              curr_Bi[j] = 1/(EOmegal[j,j]*pli[l,i])
            }else{
              curr_Bi[j] = (10^10)/(EOmegal[j,j])
            }
#            print(c("currBi",curr_Bi[j]))
          }else{        
            tb = exp(-Atmp[j,i])
            if(pli[l,i] >= 1e-10){
              bb = tb/pli[l,i]
            }else{
              bb = tb*(10^10)
            }
#            print(c("EOmegajj", j, EOmegal[j,j], tb, bb))
            cb = -EOmegal[j,j]*tb
            curr_Bi[j] = mysolve_b(b=bb,c=cb)
#            print(c("i d bb cb",i,d,bb,cb))
#            print(c("currBi",curr_Bi[j]))
          }
        }
        
        curr_Eexpi=curr_Ai + 0.5*curr_Bi
#        print(c("curr_Eexpi",curr_Eexpi))
        curr_val = compute_valLambdafn(d=d,Xi=X[,i],Eei=curr_Eexpi,Ai=curr_Ai,Bi=curr_Bi,plili=pli[l,i],EOmegal=EOmegal,EMul=EMul)
        
        
        if(curr_val>val){
 #         print(c("prev_val curr_val", val, curr_val))
          Atmp[,i]=curr_Ai
          Btmp[,i]=curr_Bi
          EexpLambdatmp[,i]=curr_Eexpi
          val=curr_val
        }else{
          loop = FALSE
        }
#        print(c("Final curr_Ai", curr_Ai))
#        print(c("Final curr_Bi",curr_Bi))
      }
#     print(c("Lambda search",i,l,iter,pli[l,i]))
    }
    
    UA[[l]]=Atmp
    UB[[l]]=Btmp
    UEexpLambda[[l]]=EexpLambdatmp
  }
  
  return(list("ELambda"=UA,"B"=UB,"EexpLambda"=UEexpLambda))
}
 
compute_valLambdafn<-function(Eei,Ai,Bi,Xi,EOmegal,EMul,plili,d){
#  print(c("Eei",Eei))
#  print(c("Ai",Ai))
#  print(c("Bi",Bi))
#  print(c("Xi",Xi))
#  print(c("plili",plili))
    
  z=0
  val = 0.5*(sum(log(Bi)))
  if(plili < 1e-10){
#    print(c("val",val))
    return (val)
  }
#  ELlLlTi=make.positive.definite(diag(Bi, nrow=d) + Ai%*%t(Ai),tol=1e-5)
  ELlLlTi=diag(Bi, nrow=d) + Ai%*%t(Ai)
  for (j in 1:d){
    z = z + (-exp(Eei[j]) + Ai[j]*Xi[j]) 
  }
  y = -0.5*matrix.trace(EOmegal %*% (ELlLlTi - EMul%*%t(Ai) - Ai%*%t(EMul)))
  val = val + plili*(z + y)
#  print(c("val",val))
  return(val)
}

#BELOW function is not used
compute_Lambda_aux<-function(ELambda,B,K,d,n){
  EexpLambda=list()
  ELlLlT=list()
  for(l in 1:K){
    ELl = ELambda[[l]]
    Bl = B[[l]]
    EexpLambda[[l]] = Elambda[[1]]+0.5*B[[l]]
    tmp=list()
    for(i in 1:n){
      tmp[[i]] = make.positive.definite(diag(Bl[,i],nrow=d) + ELl[,i]%*%t(ELl[,i]),tol=1e-5)
    }
    ELlLlT[[l]]=tmp
  }
  return(list("EexpLambda"=EexpLambda,"ELlLlT"=ELlLlT))
}

compute_Pi<-function(pli,n,K){
  Pi=matrix(rep(0,K),ncol=1)
  for(l in 1:K){
    Pi[l] = (1/n)*(sum(pli[l,]))
  }
  #  print(c("PI",Pi))
  return(Pi)
}

compute_EOmega<-function(n,wO,invWOmega,K){
  EOmega=list()
  for(l in 1:K){
      EOmega[[l]] = wO[[l]]*invWOmega[[l]]
  }
  return (EOmega)
}

compute_Tu<-function(beta,EOmega,pli,K,d){
  Tu=list()
  for(l in 1:K){
    Tu[[l]] = make.positive.definite(diag(x=beta,nrow=d) + EOmega[[l]]*(sum(pli[l,])),tol=1e-10)
  }
  return(Tu)
}

compute_EMu<-function(Tu,EOmega,K,n,ELambda,pli){
  EMu=list()
  for(l in 1:K){
    X=ELambda[[l]]
    x=vector("numeric",length=nrow(X))
    for(i in 1:n){
      x = x + X[,i]*pli[l,i]
    }
    EMu[[l]] = (solve(Tu[[l]]))%*%(EOmega[[l]]%*% x)
  }
  return(EMu)
}

compute_wO<-function(v,pli,K){
  wO=list()
  for(l in 1:K){
    wO[[l]] = v + sum(pli[l,]) 
  }
  return(wO)
}

compute_WOmega<-function(V,K,n,pli,EMu,EMlMlT,ELambda,B,d){
  WOmega=list()
  for(l in 1:K){
    z=0
    ELl=ELambda[[l]]
    Bl=B[[l]]
    for(i in 1:n){
#      ELlLlTi = make.positive.definite(diag(Bl[,i], nrow=d) + ELl[,i]%*%t(ELl[,i]),tol=1e-5)
      ELlLlTi = diag(Bl[,i], nrow=d) + ELl[,i]%*%t(ELl[,i])
      z = z + ELlLlTi*pli[l,i] - ELl[,i]%*%(pli[l,i]*t(EMu[[l]])) - EMu[[l]]%*%(t(ELl[,i])*pli[l,i])   
    }
    WOmega[[l]] = make.positive.definite((V + z + EMlMlT[[l]]*(sum(pli[l,]))),tol=1e-10)
  }
  return(WOmega)
}


compute_wO_WOmega_invWOmega<-function(v,V,K,n,pli,EMu,EMlMlT,ELambda,B,d,method,alpha){
  wO=list()
  WOmega=list()
  invWOmega=list()
  
  if (method==0){
    # No sparse
    wO=compute_wO(v=v,pli=pli,K=K)
    WOmega=compute_WOmega(V=V,K=K,n=n,pli=pli,EMu=EMu,EMlMlT=EMlMlT,ELambda=ELambda,B=B,d=d)
    for(i in 1:K){
      invWOmega[[i]]=solve(make.positive.definite(WOmega[[i]],tol=1e-10))
    }
  }else if (method==1){
    # pcor
    wO=compute_wO(v=v,pli=pli,K=K)
    WOmega=compute_WOmega(V=V,K=K,n=n,pli=pli,EMu=EMu,EMlMlT=EMlMlT,ELambda=ELambda,B=B,d=d)
    for(i in 1:K){
      omega=pcor_thresh(Sigma=make.positive.definite(WOmega[[i]],tol=1e-10),n=n,alpha=alpha)
      invWOmega[[i]]=omega$thresh.omega
    }
  }else {
    # method: 2 (glasso EBIC), 3 (glasso CV), 4 (tiger ebic), 5 (tiger cv)
    for(i in 1:K){
      temp=solve_eta_C_rho(method=method,d=nrow(ELambda[[i]]),N=n,v=v,V=V,s=pli[i,],X=ELambda[[i]],mu=EMu[[i]],EmimiT=EMlMlT[[i]])
      wO[[i]]=temp$eta
      WOmega[[i]]=temp$C
      invWOmega[[i]]=temp$invC
    }
  }
  
  result=list("wO"=wO,"WOmega"=WOmega,"invWOmega"=invWOmega)
  return(result)
}

compute_EMlMlT<-function(Tu, EMu, K){
  EMlMlT=list()
  for(l in 1:K){
    mumut = EMu[[l]]%*%t(EMu[[l]])
    EMlMlT[[l]] = make.positive.definite(solve(Tu[[l]]) + mumut,tol=1e-5)
  }
  return(EMlMlT)
}

compute_ElnOmega<-function(WOmega,wO,d,K){
  ElnOmega=list()
  for(l in 1:K){
    dg=0
    for(s in 1:d){
      dg = dg + digamma((wO[[l]] + 1 - s)/2)
    }
    ElnOmega[[l]] = dg + d*log(2) - mydet(WOmega[[l]])
  }
  return(ElnOmega)
}

compute_pbar<-function(Pi,ElnOmega,X,EOmega,EMu,EMlMlT,ELambda,B,EexpLambda,K,n,d){
  pbar=matrix(rep(1,K*n),nrow=K)
  lpi=(-d/2)*log(2*pi)
  for(l in 1:K){
    x=ElnOmega[[l]]/2 + log(Pi[l]+(1e-20)) + lpi
    ELl=ELambda[[l]]
    EexpLambdal=EexpLambda[[l]]
    Bl=B[[l]]
    for(i in 1:n){
      z=0
#      ELlLlTi = make.positive.definite(diag(Bl[,i], nrow=d) + ELl[,i]%*%t(ELl[,i]),tol=1e-5)
      ELlLlTi = diag(Bl[,i], nrow=d) + ELl[,i]%*%t(ELl[,i])
      for (j in 1:d){
        z = z + (-exp(min(50,EexpLambdal[j,i])) - lfactorial(X[j,i]) +ELl[j,i]*X[j,i]) 
      }
      y = matrix.trace(EOmega[[l]] %*% (ELlLlTi - EMu[[l]]%*%t(ELl[,i]) - ELl[,i]%*%t(EMu[[l]]) + EMlMlT[[l]]))
      val = x + z - (1/2)*y 
      if(is.infinite(val)){
        print(c("Pbar Infinite",l,i,x,z,y))
        print(c(ElnOmega[[l]],Pi[l],lpi))
        pbar[l,i] = -100 # so that pli[l,i] will be close to 0
      }else{
        pbar[l,i] = val
      }
    }
  }
  return(pbar)  
}

mydet<-function(M){
  return(as.numeric(determinant(M,logarithm=TRUE)$modulus))
}

compute_pli<-function(P){
  K=nrow(P)
  n=ncol(P)
  
  pli=matrix(rep(0,K*n),nrow=K)
  for(i in 1:n){
    x=P[,i]
    for(l in 1:K){
      pli[l,i] = 1/(sum(exp(x-P[l,i])))
    }
  }
  return(pli)
}

pcor_thresh<-function(Sigma,n,alpha){
  d=dim(Sigma)[1]
  s=d-1
  
  #  print(c("condition number Sigma",kappa(Sigma)))
  #  print(c("Printing sigma"))
  #  print(Sigma)

  while(kappa(Sigma)>100000){
	Sigma=Sigma + diag(10*max(abs(Sigma)),nrow=d)
      print(c("Updated condition number Sigma",kappa(Sigma)))
  }
  
  #MLE estimate for Omega = inverse(Sigma)
  Omega=solve(Sigma)
  pc=mypcor(Omega) 

  #  print(c("condition number sigma and omega",kappa(Sigma),kappa(Omega)))
  
  thresh_parcors=ci_par_cor(par_cors=pc, alpha=alpha, n=n, s=s) 
  thresh_omega= (thresh_parcors$mat) * Omega

#  print(c("Sigma"))
#  print(Sigma)
#  print(c("Omega"))
#  print(Omega)
#  print(c("pc"))
#  print(pc)
#  print(c("thresh_parcors"))
# print(thresh_parcors$mat)
#  print(c("thresh_omega"))
#  print(thresh_omega)
#  print(c("kappa", kappa(thresh_omega),kappa(Omega)))
  
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
      Clean=res$icov[[ind]]
      Clean[!is.finite(Clean)]=1e-5
      res$icov[[ind]]=Clean
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

##################################################################
# -exp(x) - b*x + c
# Newton-Raphson method from numerical recipes in C after modifications to fit the function we have
mysolve_a<-function(b,c,xacc=1e-3){
#  print(c("Solve_a b c",b,c))
  if(b<0){
    stop("Error: b cannot be <0")
  }
  
  if(b<1e-10){
    if(c<=0 || abs(c)<1e-10){
      return (-50) # pegged to small negative value
    }else{
      return (log(c))
    }
  }
  
  rts=0
  #bounding the interval where the root is
  if(c==0){
    x1=-800
    x2=0
  }
  if(c>0){
    if(c>=1){
      x1=-log(1+c)
      x2=log(1+c)
    }else{
      x1=log(c)
      x2=log(1+c)
    }
  }else{
    if(c<0){
      x1=-((abs(c)+2)/b)
      x2=0
    }
  }
  
  fl=-exp(x1) - b*x1 + c
  fh=-exp(x2) - b*x2 + c
  
  
  #  if((fl>0 && fh>0)||(fl<0 && fh<0)){
  #    stop (c("Error in root finding"," b= ",b," c= ",c))
  #  }
  
  if(fl==0){
    return (x1)
  }
  if(fh==0){
    return (x2)
  }
  
  if(fl<0){
    xl=x1
    xh=x2
  }else{
    xh=x1
    xl=x2
  }
  rts=0.5*(x1+x2)
  dxold=abs(x2-x1)
  dx=dxold
  f=-exp(rts) - b*rts + c
  df=-exp(rts) - b
  
  MAXIT=100
  j=0
  while(j<=MAXIT){
    j=j+1
    if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0) || (abs(2*f)>abs(dxold*df))){
      dxold=dx
      dx=0.5*(xh-xl)
      rts=xl+dx
      if(xl==rts){
        return (rts)
      }
    }else{
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
    f=-exp(rts) - b*rts + c
    df=-exp(rts) - b
    if(f<0){
      xl=rts
    }else{
      xh=rts
    }
  }
  return (x1)
}

##################################################################
# -exp(x/2) + b*(1/x) + c
# Newton-Raphson method from numerical recipes in C after modifications to fit the function we have
mysolve_b<-function(b,c,xacc=1e-3){
#  print(c("Solve_b b c",b,c))
  if(b<0){
    stop("Error: b cannot be <0")
  }
  rts=0
  
  #bounding the interval where the root is
  # lower bound for root since we are solving for variance (non-negative quantity) and original function
  # is of form -c1*exp(x/2) -c2*x + (1/2)log(x) + constant
  
  x1=0 
  x2=min(2*(abs(b)+abs(c)), 50)
  
  fl=-exp(x1/2) + b*(1/x1) + c
  fh=-exp(x2/2) + b*(1/x2) + c
  
  
  #  if((fl>0 && fh>0)||(fl<0 && fh<0)){
  #    stop (c("Error in root finding"," b= ",b," c= ",c))
  #  }
  
  if(fl==0){
    return (x1)
  }
  if(fh==0){
    return (x2)
  }
  
  if(fl<0){
    xl=x1
    xh=x2
  }else{
    xh=x1
    xl=x2
  }
  rts=0.5*(x1+x2)
  dxold=abs(x2-x1)
  dx=dxold
  f=-exp(rts/2) + b*(1/rts) + c
  df=-(1/2)*exp(rts/2) - b*(1/(rts*rts))
  
  MAXIT=100
  j=0
  while(j<=MAXIT){
    j=j+1
    if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0) || (abs(2*f)>abs(dxold*df))){
      dxold=dx
      dx=0.5*(xh-xl)
      rts=xl+dx
      if(xl==rts){
        return (rts)
      }
    }else{
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
    f=-exp(rts/2) + b*(1/rts) + c
    df=-(1/2)*exp(rts/2) - b*(1/(rts*rts))
    if(f<0){
      xl=rts
    }else{
      xh=rts
    }
  }
  return (x1)
}

