require(mvtnorm)
require(tibble)
require(dplyr)
require(RCAL)
require(glmnet)
library(xtable)
#######################
######################

Simulation=function()
{
  Out2=NULL
  Out3=NULL
  Out4=NULL
  Out5=NULL
  for (p in p.all) {
    #####TRUE ATT###
    A=Generate(n=1000000,p=p,c )
    mu.f=A$mu.f
    D=A$D
    D.temp=D%>%select(X.mu)%>%rowwise()%>%mutate(mu.true=mu.f(x=X.mu[,1:4]))
    ATT.true= mean(D.temp$mu.true[D$Trt==1])

  for (N in N.all) {
  
  #########
  Out=NULL
  Out.sd.est=NULL
  Cover=NULL
  pip=NULL
  pr=NULL
  ahat=NULL
  ahat0=NULL
  a.true=NULL
  a.op=NULL
  a.boot=NULL
  n.temp=NULL
   for (t in 1:T) {
    
  
    
    t1=Sys.time()
    ####generate data#######
    A=Generate(n=N,p=p ,c)
    D=A$D
    mean.X=apply(D$X,2 ,mean)
    D=D%>%mutate(X.c=X-mean.X)
    
    prop_p=A$prop_p
    prop_delta=A$prop_delta
    mu.f=A$mu.f
    
    
    D.in=D%>%dplyr::filter(R==1)
    pip=c(pip,mean(D$Trt))
    n.temp=c(n.temp, sum(D$R))
    ###########Ture##############
    D=D%>%rowwise()%>%mutate(p.true=prop_p(x=X.p[,1:4]))
    D=D%>%rowwise()%>%mutate(delta.true=prop_delta(x=X.p[,1:4]))
    D=D%>%rowwise()%>%mutate(mu.true=mu.f(x=X.mu[,1:4]))
    #######################delta######################

    fit.cv=glm.regu.cv(fold=5,y=D$inv_Trt,x=D$X,nrho = 10,ann.init = FALSE)
    
     D.temp= D%>%select(X)%>% rowwise()%>%
     mutate(delta=(1+exp(fit.cv$sel.bet[1,1]+sum(X*fit.cv$sel.bet[2:(p+1),1])))^(-1))
    
     D=D%>%tibble(delta=D.temp$delta)
    
    D.temp= D%>% select(X)%>%rowwise()%>%
      mutate(weight=exp(-sum(X*fit.cv$sel.bet[2:(p+1),1])))
    
    D=D%>%tibble(weight=D.temp$weight)

    #######################p#######################
    fit.p=glm.regu.cv(fold=5,y=D.in$inv_Trt,x=D.in$X,nrho = 10,ann.init = FALSE)
    
   
    D.temp= D%>%select(X)%>%  rowwise()%>% 
      mutate(p=(1+exp(fit.p$sel.bet[1,1]+sum(X*fit.p$sel.bet[2:(p+1),1])))^(-1))

      D=D%>%tibble(p=D.temp$p)

   
    D.temp= D%>%select(X)%>%  rowwise()%>%
      mutate(weight.beta=exp(-sum(X*fit.p$sel.bet[2:(p+1),1])))
    
    D=D%>%tibble(weight.beta=D.temp$weight.beta)
     
    #####################mu#######################
    D.0=D%>%dplyr::filter(Trt==0)
    mean.X=apply(D.0$X,2 ,mean)
    D.0=D.0%>%mutate(X.c=X-mean.X)
    alpha.0=mean(D.0$Y)
    
    
    
    fit.mu.cv=cv.glmnet(y=D.0$Y-alpha.0,x=D.0$X.c,
           intercept=FALSE,weights = D.0$weight,nfolds = 5,family="gaussian" )
    
    D=D%>%mutate(mu=predict(fit.mu.cv,newx=D$X)-predict(fit.mu.cv,newx=mean.X )[1]   +mean(D.0$Y))
    
     ##################mu.1#############
    
    fit.mu.cv.1=cv.glmnet(y=D.0$Y,x=D.0$X,
         intercept=TRUE,weights = D.0$weight,nfolds = 5,family="gaussian" )
    
    D=D%>%mutate(mu.1=predict(fit.mu.cv.1,newx = D$X))
    #####################mu.beta#######################
    
    D.01=D%>%dplyr::filter(Trt==0,R==1)
    mean.X=apply(D.01$X,2 ,mean)
    D.01=D.01%>%mutate(X.c=X-mean.X)
    fit.mu.cv.p=cv.glmnet(y=D.01$Y-mean(D.01$Y),x=D.01$X.c,weights = D.01$weight.beta ,
                          intercept=FALSE,nfolds = 5,family="gaussian" )
    D=D%>%mutate(mu.beta=predict(fit.mu.cv.p,newx=D$X)-predict(fit.mu.cv.p,newx=mean.X )[1]   +mean(D.01$Y))
    
      #####################Control Effect######################
    phi.eff=function(U,Y,mu,delta)
    {return(U*mu/pip[t]+(1-U)*delta*(Y-mu)/pip[t]/(1-delta) )}
    phi.nv=function(R,Trt,Y,mu,p)
    {return(R*Trt*mu/pip[t]+R*(1-Trt)*p*(Y-mu)/pip[t]/(1-p) )}
    phi.or= function(U,mu) {
      return(U*mu/pip[t])
    }

    ##############Eff##########################
    D=D%>%rowwise()%>%mutate(phi.eff=phi.eff(U=U,Y=Y,mu=mu,delta=delta))
    D=D%>%rowwise()%>%mutate(phi.eff.1=phi.eff(U=U,Y=Y,mu=mu.1,delta=delta))
    ###############NV###########################
    
    D=D%>%rowwise()%>%mutate(phi.nv=phi.nv(R=R,Trt=Trt,Y=Y,mu=mu.beta,p=p))
    D=D%>%rowwise()%>%mutate(phi.nv.0=phi.nv(R=R,Trt=Trt,Y=Y,mu=mu,p=p))
    
    ###############NV.old###########################
    
    #################OR##########################
    D=D%>%rowwise()%>%mutate(phi.or=phi.or(U=U,mu=mu))
    D=D%>%rowwise()%>%mutate(phi.or.1=phi.or(U=U,mu=mu.1))
    D=D%>%rowwise()%>%mutate(phi.or.beta=phi.or(U=U,mu=mu.beta))
    
    ##################a##################
    aN=function(R,Trt,Y,delta,p,mu)
    {return(R*(1-Trt)*(1-delta/p)*p^2*(Y-mu)^2/(1-delta)/(1-p)^2 /pip[t] )}
    aD=function(R,Trt,Y,delta,p,mu)
    {return(            
      ((1-R)*(delta/p)^2*(1-p)^2+R*(1-Trt)*(1-delta/p)^2)*p^2*(Y-mu)^2/(1-delta)^2/(1-p)^2 /pip[t]
    )}
    
    D=D%>%rowwise()%>%mutate(aN=aN(R=R,Trt=Trt,Y=Y,delta=delta,p=p,mu=mu.beta),
                             aD=aD(R=R,Trt=Trt,Y=Y,delta=delta,p=p,mu=mu.beta))
    
    a=mean(D$aN)/mean(D$aD)
    ahat=c(ahat,a)
   
     D=D%>%rowwise()%>%mutate(aN=aN(R=R,Trt=Trt,Y=Y,delta=delta,p=p,mu=mu),
                             aD=aD(R=R,Trt=Trt,Y=Y,delta=delta,p=p,mu=mu))
    
    a3=mean(D$aN)/mean(D$aD)
    ahat0=c(ahat0,a3)
    
    
    
    
    D=D%>%rowwise()%>%mutate(aN.true=aN(R=R,Trt=Trt,Y=Y,delta=delta.true,p=p.true,mu=mu.true),
                             aD.true=aD(R=R,Trt=Trt,Y=Y,delta=delta.true,p=p.true,mu=mu.true))
    
    a2=mean(D$aN.true)/mean(D$aD.true)
    a.true=c(a.true,a2)
    
    
    ###################
    
    f.mu=function(a)
    {return( mean(( (1-a)*(D$phi.nv.0-mean(D$phi.nv.0)  )
                    +a*(D$phi.eff-mean(D$phi.eff)    )     )^2 ) )}
    
    a4=optimize(f=f.mu,lower=0,upper = 3,maximum = FALSE)$minimum
    
    a.op=c(a.op,a4)
    
    ########################a.boot########################
    a.can=(-15):15/10
    a.T=1000
    NN=dim(D)[1]
    a.temp=NULL
    
    for (k in 1:a.T) {
      S=sample(1:NN,replace = TRUE)
      EFF=mean(D$phi.eff[S])
      NV.0=mean(D$phi.nv[S])
      a.temp= rbind(a.temp, a.can*EFF+(1-a.can)*NV.0)
    }
    
    a.boot.temp=a.can[which.min( apply(a.temp,2,sd))]
    a.boot=c(a.boot,  a.boot.temp)
    
    
    
    
    
    ################################################
    
    D=D%>%rowwise()%>%mutate(phi.a=a*phi.eff+(1-a)*phi.nv)
    D=D%>%rowwise()%>%mutate(phi.a0=a3*phi.eff+(1-a3)*phi.nv.0)
    D=D%>%rowwise()%>%mutate(phi.a.star=a2*phi.eff+(1-a2)*phi.nv)
    D=D%>%rowwise()%>%mutate(phi.a.op=a4*phi.eff+(1-a4)*phi.nv.0)
    D=D%>%rowwise()%>%mutate(phi.a.boot=a.boot.temp*phi.eff+(1-a.boot.temp)*phi.nv.0)
    ####################################
    Method=c("phi.or","phi.or.1","phi.or.beta","phi.eff","phi.eff.1",
             "phi.nv","phi.nv.0","phi.a","phi.a0","phi.a.star","phi.a.op","phi.a.boot"
    )
    
    
    ###################################
   ATT=D%>%select(Method)%>% group_by()%>%summarise_all(mean)
   colnames(ATT)=c("OR","OR.1","OR.beta","EFF","EFF.1","NV","NV.0","A","A.0","A.star","A.op","A.boot")
   
     ###################################

  #################sd####################
   sd.est=D%>%select(Method)%>% group_by()%>%summarise_all(sd)
   
   sd.est= ((sd.est)*sqrt(N-1)/sqrt(N) ) /sqrt(N)
   sd.est[1]= sqrt(mean((D$phi.or-ATT.true)^2)/(N))
   sd.est[2]= sqrt(mean((D$phi.or.1-ATT.true)^2)/(N))

   colnames(sd.est)=  colnames(ATT)
   #################################
   
   CI.U=ATT+qnorm(0.975)*sd.est
   CI.L=ATT-qnorm(0.975)*sd.est
   
   Cover=  bind_rows(Cover,(sign(( CI.U-ATT.true)* (CI.L-ATT.true))*(-1)+1)/2)
    

   
   
  ####################################
    Out=bind_rows(Out,ATT)
    Out.sd.est=bind_rows(Out.sd.est,sd.est)
    
    pr=c(pr,mean(D.in$Trt))
    
    ##################################
    
    
    
    
    
    t2=Sys.time()
    
 #   p.err=c(p.err,mean(((D$p.true-D$p)/D$p.true)^2)  )
    #delta.err=c(delta.err,mean(((D$delta.true-D$delta)/D$delta.true)^2)  )
    
    #eff.err=c(eff.err,mean(D$eff.test))
    #nv.err=c(nv.err,nv.test=mean(D$nv.test))
    
    print(paste( "The time period of", "N=",N  ,"p=", p   ,  t,"-th",  " simulation is ",t2-t1))
    
  }
  
  


  
  ATT=Out%>%group_by()%>%summarise_all(mean)
  SD.est=Out.sd.est%>%group_by()%>%summarise_all(mean)
  SD=Out%>%group_by()%>%summarise_all(sd)
  Bias.r=100*(ATT-ATT.true)/ATT.true
  Bias=(ATT-ATT.true)
  
  Cover.95=apply(Cover,2,mean)
  
  ATT=round(ATT,digits = 3)
  SD.est=round(SD.est,digits = 3)
  SD=round(SD,digits = 3)
  Bias=round(Bias,digits = 3)
  Bias.r=round(Bias.r,digits = 3)
  Cover.95=round(Cover.95,digits=3)
  
  
  Out2.temp=NULL
  Out2.temp=rbind(Out2.temp,c(Case=Case,c=c,N=N,n=round(mean(n.temp),digits = 3),
                              p=p,pr=round(mean(pr),digits = 3), Bias
  ))
  
  
  Out2.temp=rbind(Out2.temp,c(Case=Case,c=c,N=N,n=round(mean(n.temp),digits = 3),
                              p=p,pr=round(mean(pr),digits = 3), Bias.r) )
  
  Out2.temp=rbind(Out2.temp,c(Case=Case,c=c,N=N,n=round(mean(n.temp),digits = 3),
                              p=p,pr=round(mean(pr),digits = 3), SD.est)) 
  Out2.temp=rbind(Out2.temp,c(Case=Case,c=c,N=N,n=round(mean(n.temp),digits = 3),
                              p=p,pr=round(mean(pr),digits = 3), SD))
  Out2.temp=rbind(Out2.temp,c(Case=Case,c=c,N=N,n=round(mean(n.temp),digits = 3),
                              p=p,pr=round(mean(pr),digits = 3), Cover.95))
  
  rownames(Out2.temp)=c("Bias","Bias.r","SD.est","SD","Cover.95")
  
  
  
  Out3.temp= tibble(Case=Case,N=N,n=round(mean(n.temp),digits = 3),
                    p=p,pr=round(mean(pr),digits = 3),ATT.true=ATT.true,a=mean(ahat),
                    a.true=mean(a.true),
                    a0=mean(ahat0),a.op=mean(a.op),
                    a.boot=mean(a.boot))
  
   Out4.temp=Out2.temp[,c(1:6,9,11,13,14,15)]
  
   Out2=rbind(Out2,Out2.temp)
   Out3=rbind(Out3,Out3.temp)
   Out4=rbind(Out4,Out4.temp)
   write.table(Out2, file =  paste0(Case,ratio,name,".csv"), sep = ",", col.names=NA)
   
  }
  }
  
print(Out2)
  print(Out3)
  print(Out4)
  
  
  print(xtable(Out2.temp,digits=3))
  print(xtable(Out3,digits=3))
  print(xtable(Out4,digits=3))
  
  
  write.table(Out2, file =  paste0(Case,ratio,name,".csv"), sep = ",", col.names=NA)
  
# A=read_csv("S10.4.csv")
  
  
}
