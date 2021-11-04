

Generate.S1=function(p,n,c)
{
  S=diag(p)
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      S[i,j]=2^(-(j-i) )
      S[j,i]=S[i,j]
    }
  }
  ###################################
  
  
  prop_p=function(x){return((1+exp(-x[1]-0.5*x[2]-0.25*x[3]-0.125*x[4]))^(-1))}
  prop_delta=function(x){return((1+exp(c-x[1]-0.5*x[2]-0.25*x[3]-0.125*x[4]))^(-1))}
  mu=function(x){return(x[1]+0.5*x[2]+0.25*x[3]+0.125*x[4])}
  ###################################
  
  X=mvtnorm::rmvnorm(n,sigma = S)
  Y=apply(X[,1:4],1,mu)+rnorm(n)
  P_p=apply(X[,1:4],1,prop_p)
  P_delta=apply(X[,1:4],1,prop_delta)
  P_pi=P_delta/P_p
  Trt=rbinom(n,1,P_p)
  R=rbinom(n,1,P_pi)
  for (i in 1:n) {
    if(R[i]==0)
    {Trt[i]=0}
  }
  
  return(list(D=tibble(Y=Y,X=X,X.p=X[,1:4],X.mu=X[,1:4],Trt=Trt,R=R,inv_Trt=abs(Trt-1),U=Trt*R),
              prop_p=prop_p,prop_delta=prop_delta,mu.f=mu))
}

Generate.S2=function(p,n,c)
{
  S=diag(p)
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      S[i,j]=2^(-(j-i) )
      S[j,i]=S[i,j]
    }
  }
  ######################################
  X=mvtnorm::rmvnorm(n,sigma = S)
  #####################################
  tran=function(a){return(a+( max(a+1,0) )^2)}
  Z=matrix(0,n,4)
  for (i in 1:n ) {
    for(j in 1:4)
    {
      
      Z[i,j]=tran(X[i,j])
    }
  }
  
  Z=scale(Z)
  ######################################
  prop_p=function(x){return((1+exp(-x[1]-0.5*x[2]-0.25*x[3]-0.125*x[4]))^(-1))}
  prop_delta=function(x){return((1+exp(c-x[1]-0.5*x[2]-0.25*x[3]-0.125*x[4]))^(-1))}
  mu=function(x){return(x[1]+0.5*x[2]+0.25*x[3]+0.125*x[4])}
  ###################################

  Y=apply(Z[,1:4],1,mu)+rnorm(n)
  
  P_p=apply(X[,1:4],1,prop_p)
  P_delta=apply(X[,1:4],1,prop_delta)
  P_pi=P_delta/P_p
  Trt=rbinom(n,1,P_p)
  R=rbinom(n,1,P_pi)
  for (i in 1:n) {
    if(R[i]==0)
    {Trt[i]=0}
  }
  
  
  
  return(list(D=tibble(Y=Y,X=X,X.p=X,X.mu=Z,Trt=Trt,R=R,inv_Trt=abs(Trt-1),U=Trt*R),
              prop_p=prop_p,prop_delta=prop_delta,mu.f=mu))
}

Generate.S3=function(p,n,c)
{
 
  S=diag(p)
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      S[i,j]=2^(-(j-i) )
      S[j,i]=S[i,j]
    }
  }
  ######################################
  X=mvtnorm::rmvnorm(n,sigma = S)
  #####################################
  tran=function(a){return(a+( max(a+1,0) )^2)}
 Z=matrix(0,n,4)
   for (i in 1:n ) {
   for(j in 1:4)
   {
     
     Z[i,j]=tran(X[i,j])
   }
 }

 Z=scale(Z)
  #####################################
  #####################################
  prop_p=function(x){return(0.7)}
  prop_delta=function(x){return(0.7*(1+exp(c-x[1]-0.5*x[2]-0.25*x[3]-0.125*x[4]))^(-1))}
  mu=function(x){return(x[1]+0.5*x[2]+0.25*x[3]+0.125*x[4])}
  ##################################
  
  Y=apply(X[,1:4],1,mu)+rnorm(n)
  
  P_p=apply(Z[,1:4],1,prop_p)
  P_delta=apply(Z[,1:4],1,prop_delta)
  P_pi=P_delta/P_p
  Trt=rbinom(n,1,P_p)
  R=rbinom(n,1,P_pi)
  for (i in 1:n) {
    if(R[i]==0)
    {Trt[i]=0}
  }

  return(list(D=tibble(Y=Y,X=X,X.p=Z,X.mu=X,Trt=Trt,R=R,inv_Trt=abs(Trt-1),U=Trt*R),
              prop_p=prop_p,prop_delta=prop_delta,mu.f=mu))
  
}

##############################


Generate.S4=function(p,n,c)
{
  S=diag(p)
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      S[i,j]=2^(-(j-i) )
      S[j,i]=S[i,j]
    }
  }
  

  ######################################
  X=mvtnorm::rmvnorm(n,sigma = S)
  #####################################
  tran=function(a){return(a+( max(a+1,0) )^2)}
  Z=matrix(0,n,4)
  for (i in 1:n ) {
    for(j in 1:4)
    {
      
      Z[i,j]=tran(X[i,j])
    }
  }
  
  Z=scale(Z)
  ######################################
  prop_p=function(x){return((1+exp(-x[1]-0.5*x[2]-0.25*x[3]-0.125*x[4]))^(-1))}
  prop_delta=function(x){return((1+exp(c-x[1]-0.5*x[2]-0.25*x[3]-0.125*x[4]))^(-1))}
  mu=function(x){return(x[1]+0.5*x[2]+0.25*x[3]+0.125*x[4])}
  ###################################
  
  Y=apply(X[,1:4],1,mu)+rnorm(n)
  
  P_p=apply(Z[,1:4],1,prop_p)
  P_delta=apply(Z[,1:4],1,prop_delta)
  
  P_pi=P_delta/P_p
  Trt=rbinom(n,1,P_p)
  R=rbinom(n,1,P_pi)
  for (i in 1:n) {
    if(R[i]==0)
    {Trt[i]=0}
  }
  
  
  
  return(list(D=tibble(Y=Y,X=X,X.p=Z,X.mu=X,Trt=Trt,R=R,inv_Trt=abs(Trt-1),U=Trt*R),
              prop_p=prop_p,prop_delta=prop_delta,mu.f=mu))
}


  
  
  

