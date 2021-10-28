
data.in<-foreach(j=1:length(links)) %do% {
  data.in=as.numeric(FourierAnom[[j]][,2]!=0);
}
data.out<-foreach(j=1:length(links)) %do% {
  data.out=as.numeric(FourierAnom[[j]][,3]!=0);
}

data<-foreach(j=1:length(links)) %do% {
  data=as.data.frame(matrix(as.numeric(FourierAnom[[j]][,2:3]!=0),ncol=2,dimnames=list(list(),list("Incoming","Outgoing"))));
}
names(data)=links

spaces.in<-foreach(j=1:length(links)) %do% {
  spaces<-c();
  newspace<-1;
  oldspace<-1;
  
  for (i in 1:length(data.in[[j]])) {
    
    if (data.in[[j]][i]==1) {
      
      newspace<-i;
      
    }
    
    if (newspace!=oldspace) {
      
      spaces<-c(spaces,newspace-oldspace);
      oldspace<-newspace;
      
    }
    
  }
  spaces<-c(spaces,length(data.in[[j]][sum(spaces):length(data.in[[j]])]));
  spaces.in<-spaces
  
}



spaces.out<-foreach(j=1:length(links)) %do% {
  spaces<-c();
  newspace<-1;
  oldspace<-1;
  
  for (i in 1:length(data.out[[j]])) {
    
    if (data.out[[j]][i]==1) {
      
      newspace<-i;
      
    }
    
    if (newspace!=oldspace) {
      
      spaces<-c(spaces,newspace-oldspace);
      oldspace<-newspace;
      
    }
    
  }
  spaces<-c(spaces,length(data.out[[j]][sum(spaces):length(data.out[[j]])]));
  spaces.out<-spaces
  
}




#################generate alternating renewal processes
m=2;


R.in<-foreach(j=1:length(links)) %do% {
  R.in=spaces.in[[j]][spaces.in[[j]]>m];
}

R.out<-foreach(j=1:length(links)) %do% {
  R.out=spaces.out[[j]][spaces.out[[j]]>m];
}


A.in<-foreach(j=1:length(links)) %do% {
  A=c();
  nzeros=0;
  counter=1;
  n=1;
  while (n <=length(spaces.in[[j]]) ) {
    while (spaces.in[[j]][n]<=m & n<=length(spaces.in[[j]])) {
      counter=counter+spaces.in[[j]][n];
      if (spaces.in[[j]][n]>1) {
        nzeros=nzeros+spaces.in[[j]][n]-1;
      }
      n=n+1;
    }
    A=c(A,counter);
    counter=1;
    nzeros=0;
    n=n+1;
  }
  A=A[2:length(A)];
  A.in=A;
}

A.out<-foreach(j=1:length(links)) %do% {
  A=c();
  nzeros=0;
  counter=1;
  n=1;
  while (n <=length(spaces.out[[j]]) ) {
    while (spaces.out[[j]][n]<=m & n<=length(spaces.out[[j]])) {
      counter=counter+spaces.out[[j]][n];
      if (spaces.out[[j]][n]>1) {
        nzeros=nzeros+spaces.out[[j]][n]-1;
      }
      n=n+1;
    }
    A=c(A,counter);
    counter=1;
    nzeros=0;
    n=n+1;
  }
  A=A[2:length(A)];
  A.out=A;
}


S.in=foreach(j=1:length(links)) %do% {
  S.in=R.in[[j]][-1]+A.in[[j]]-1;
}

S.out=foreach(j=1:length(links)) %do% {
  S.out=R.out[[j]][-1]+A.out[[j]]-1;
}

T.in=foreach(j=1:length(links)) %do% {
  S=c(S.in[[j]][1],rep(0,length(S.in[[j]])-1));
  for (i in 2:length(S.in[[j]]) ) {
    S[i]=S.in[[j]][i]+S[i-1];
  }
  T.in=S-(1:length(S))*mean(S.in[[j]]);
}


T.out=foreach(j=1:length(links)) %do% {
  S=c(S.out[[j]][1],rep(0,length(S.out[[j]])-1));
  for (i in 2:length(S.out[[j]]) ) {
    S[i]=S.out[[j]][i]+S[i-1];
  }
  T.out=S-(1:length(S))*mean(S.out[[j]]);
}







par(mfrow=c(2,2));
acf(log(S));
acf(S^(-1));
acf(S^(-2));
acf(S^(-3));


##################find strings during activity periods
M=20;
aic<-rep(0,M-1);
y_length=rep(0,M-1);
for (m in 2:M) {
  R=c();
  A=c();
  props=c();
  nzeros=0;
  m=2;
  R=spaces[spaces>m];
  
  counter=1;
  
  n=1;
  
  while (n <=length(spaces) ) {
    
    
    while (spaces[n]<=m & n<=length(spaces)) {
      counter=counter+spaces[n];
      if (spaces[n]>1) {
        nzeros=nzeros+spaces[n]-1;
      }
      n=n+1;
    }
    
    A=c(A,counter);
    props=c(props,1-nzeros/counter);
    counter=1;
    nzeros=0;
    n=n+1;
    
  }
  A=A[2:length(A)];
  props=props[2:length(props)]
  
  
  
  
  
  a=max(A);
  X=matrix(2,ncol=length(A),nrow=a);
  S=R[1];
  for (i in 1:(length(A)-1)) {
    X[1:A[i],i]=data[(S+1):(S+A[i])];
    S=S+R[i+1]+A[i]-1;
  }
  
  X[1:A[length(A)],length(A)]=data[(S+1):(S+A[length(A)])];
  
  
  #####################logit regression
  time<-c();
  length<-c();
  y<-c();
  for (i in 1:(length(A)-1)) {
    if (A[i]>1) {
      y<-c(y,X[2:A[i],i]);
      time<-c(time,2:A[i]);
      length<-c(length,rep(A[i],A[i]-1));
    }  
    
  }
  
  mat<-data.frame(
    y = y,
    t = time,
    A = length
  )
  
  fit<-glm(y ~ t + A,data=mat,family=binomial);
  
  #y_prob=as.vector(predict(fit,newdata=NULL,type="response"));
  #y_pred=as.numeric(y_prob>.5);
  #sum(abs(y_pred-y))/length(y);
  
  #plot(y_prob[1:(sum(A[1:20]))])
  y_length[m-1]=length(y);
  aic[m-1]=fit$aic;
}

par(mfrow=c(1,2));
plot(aic);
plot(y_length);

plot(aic/log(y_length));


#############likelihood calculation/estimation
#############for now we can ignore modeling the traffic during anomaly periods
j=1; #which link to use
start_time=Sys.time()
M=5;
my_par=matrix(0,nrow=M,ncol=4);
hieu_par=matrix(0,nrow=M,ncol=5);
my_llik=rep(0,M);
hieu_llik=rep(0,M);

for (m in 1:M) {
  #note m=1 corresponds to anomally period containing consecutive ones
  
  R=c();
  A=c();
  props=c();
  nzeros=0;
  R=spaces[spaces>m];
  
  counter=1;
  
  n=1;
  
  while (n <=length(spaces) ) {
    
    
    while (spaces[n]<=m & n<=length(spaces)) {
      counter=counter+spaces[n];
      if (spaces[n]>1) {
        nzeros=nzeros+spaces[n]-1;
      }
      n=n+1;
    }
    
    A=c(A,counter);
    props=c(props,1-nzeros/counter);
    counter=1;
    nzeros=0;
    n=n+1;
    
  }
  A=A[2:length(A)];
  props=props[2:length(props)]
  
  
  S=R[-1]+A-1;
  
  
  my_log_likelihood=function(x) {
    if(x[1]<=0 || x[2]<=1 || x[3]<=0 || x[3]>=min(R) || x[4]<=0) {
      return(10^(10));
    } else {
      
      result=-sum(log(x[1]*x[2]*(log(R/x[3]))^(x[2]-1)*exp(-x[1]*(log(R/x[3]))^x[2])/R));
      result=result-sum(log(x[4]/(A^(x[4]+1))));
      
      return(result);
    }
  }
  temp_par=matrix(0,nrow=100,ncol=4);
  temp=rep(0,100);
  for(k in 1:100) {
    u=runif(4,0,5);
    
    bs=my_log_likelihood(u);
    if (bs==Inf || is.nan(bs)==TRUE) {
      temp[k]=-100000;
    } else {
      fit_my=optim(u,my_log_likelihood);
      temp[k]=-fit_my$value;
      temp_par[k,]=fit_my$par
    }
  }
  my_llik[m]=max(temp);
  my_par[m,]=temp_par[which.max(temp),];
  
  hieu_log_likelihood=function(x) {
    
    if(x[1]<0 || x[1]>1 || x[2]<=0 || x[3]<=0 || x[4]<=1 || x[5]<=0) {
      return(10^(20));
    } else {
      
      result=-sum(log(x[1]*(x[2]/x[3])*(S/x[3])^(x[2]-1)*exp(-(S/x[3])^x[2])+(1-x[1])*2*gamma((x[4]+1)/2)/(gamma(x[4]/2)*sqrt(x[4]*pi*x[5]))*(1+S^2/(x[4]*x[5]))^(-(x[4]+1)/2)));
      
      return(result);
      
    }
  }
  
  temp_par=matrix(0,nrow=100,ncol=5);
  temp=rep(0,100);
  for(k in 1:100) {
    u1=runif(2,0,1);
    u=runif(3,0,10);
    u=c(u1[1],u[1],u[2],u1[2],u[3]);
    bs=hieu_log_likelihood(u);
    if (bs!=Inf || is.nan(bs)==FALSE) {
      fit_hieu=optim(u,hieu_log_likelihood);
      temp[k]=-fit_hieu$value;
      temp_par[k,]=fit_hieu$par
    } else {
      temp[k]=-100000;
    }
  }
  hieu_llik[m]=max(temp);
  hieu_par[m,]=temp_par[which.max(temp),];
  
}
plot(my_llik,type="b",col="blue",pch='+',ylim=c(-3000,0),xlab="M",ylab="Log Likelihood");
lines(hieu_llik,type="b",col="black");
stop_time=Sys.time()
stop_time-start_time
###############################################


hieu_model=function(y,x) {
  
  result=x[1]*(x[2]/x[3])*(y/x[3])^(x[2]-1)*exp(-(y/x[3])^x[2])+(1-x[1])*2*gamma((x[4]+1)/2)/(gamma(x[4]/2)*sqrt(x[4]*pi*x[5]))*(1+y^2/(x[4]*x[5]))^(-(x[4]+1)/2)
  return(result);
}

hieu_model_d=function(y) {
  return(hieu_model(y,fit_hieu$par));
}


hist(S,breaks=seq(1,max(S)+10,10),freq=FALSE);
t=1:4000;
lines(t,hieu_model(t,fit_hieu_sim),col='blue',lwd=2);

######Simulation Study to test likelihood
##my model####
N=100;

n=length(A);
est_my=my_par[3,]

Lam=rep(0,N);
Sig=rep(0,N);
Nu=rep(0,N);
Beta=rep(0,N);
b0=rep(0,N); b1=b0; b2=b0;

for ( j in 1:N ) {
  
  R_sim=rPPS(n,lam=est_my[1],sc=est_my[3],v=est_my[2]);
  A_sim=round(rpareto(n,shape=est_my[4]),0);
  
  a=max(A_sim);
  X_sim=matrix(2,ncol=length(A_sim),nrow=a);
  S=R_sim[1];
  for (i in 1:(length(A)-1)) {
    X_sim[1:A_sim[i],i]=rbinom(A_sim[i],1,inv.logit(est_my[5]+est_my[6]*(1:A_sim[i])+est_my[7]*A_sim[i]));
    S=S+R_sim[i+1]+A_sim[i]-1;
  }
  
  X_sim[1:A_sim[length(A_sim)],length(A_sim)]=rbinom(A_sim[length(A_sim)],1,inv.logit(est_my[5]+est_my[6]*(1:A_sim[i])+est_my[7]*A_sim[i]));
  
  
  logistic_log_lik_sim<-function(b) {
    result=0;
    for (i in 1:length(A_sim)) {
      if (A_sim[i]>1) {
        result=result-sum(dbinom(X_sim[2:A_sim[i],i],1,inv.logit(b[1]+b[2]*(2:A_sim[i])+b[3]*A_sim[i]),log=TRUE));
      }
    }
    return(result);
  }
  
  my_log_likelihood_sim=function(x) {
    if(x[1]<=0 || x[2]<=1 || x[3]<=0 || x[3]>=min(R_sim) || x[4]<=0) {
      return(10^(10));
    } else {
      
      result=-sum(log(x[1]*x[2]*(log(R_sim/x[3]))^(x[2]-1)*exp(-x[1]*(log(R_sim/x[3]))^x[2])/R_sim));
      result=result-sum(log(x[4]/(A_sim^(x[4]+1))));
      result=result+logistic_log_lik_sim(x[5:7]);
      return(result);
    }
  }
  
  temp=rep(0,50);
  temp_par=matrix(0,nrow=50,ncol=7);
  for (k in 1:50) {
    
    u=c(runif(4,0,2),0,0,0);
    
    bs=my_log_likelihood_sim(u);
    if (bs==Inf || is.nan(bs)==TRUE) {
      temp[k]=-100000;
    } else {
      fit_my=optim(u,my_log_likelihood_sim);
      temp_par[k,]=fit_my$par;
      temp[k]=-fit_my$value;
    }
  }
  pars=temp_par[which.max(temp),]
  
  Lam[j]=pars[1];
  Sig[j]=pars[3];
  Nu[j]=pars[2];
  Beta[j]=pars[4];
  b0[j]=pars[5];
  b1[j]=pars[6];
  b2[j]=pars[7];
  
}

mean(Lam);    sd(Lam);
mean(Sig);    sd(Sig);
mean(Nu);     sd(Nu);
mean(Beta);   sd(Beta);
mean(b0);     sd(b0);
mean(b1);     sd(b1);
mean(b2);     sd(b2);


####Hieu Model
N=20;
est_hieu=hieu_par[1,];
est_hieu[4]=3;
est_hieu[5]=1;
P=rep(0,N); K=P; Lambda=P; New=P; Sigma=P;

for (j in 1:N) {
  S_sim=rep(0,n);
  y=rbinom(n,1,est_hieu[1]);
  r=sum(y);
  
  if(r==0) {
    S_sim=abs(rt(n,est_hieu[4])*sqrt(est_hieu[5]));
  } else {
    S_sim[1:r]=rweibull(r,shape=est_hieu[2],scale=est_hieu[3]);
    S_sim[(r+1):n]=abs(rt(n-r,est_hieu[4])*sqrt(est_hieu[5]));
  }
  
  
  hieu_log_likelihood_sim=function(x) {
    
    if(x[1]<0 || x[1]>1 || x[2]<=0 || x[3]<=0 || x[4]<=1 || x[5]<=0 ) {
      return(10^(20));
    } else {
      
      result=-sum(log(x[1]*(x[2]/x[3])*(S_sim/x[3])^(x[2]-1)*exp(-(S_sim/x[3])^x[2])+(1-x[1])*2*exp(lgamma((x[4]+1)/2)-lgamma(x[4]/2))/sqrt(pi*x[5]*x[4])*(1+S_sim^2/(x[4]*x[5]))^(-(x[4]+1)/2)));
      
      return(result);
      
    }
  }
  
  
  
  temp=rep(0,50);
  temp_par=matrix(0,nrow=50,ncol=5);
  for (k in 1:50) {
    
    u1=runif(1,0,1);
    u=runif(4,0,10);
    u=c(u1,u);
    bs=hieu_log_likelihood_sim(u);
    if (bs==Inf || is.nan(bs)==TRUE) {
      temp[k]=-100000;
    } else {
      fit_hieu=optim(u,hieu_log_likelihood_sim);
      temp_par[k,]=fit_hieu$par;
      temp[k]=-fit_hieu$value;
    }
  }
  pars=temp_par[which.max(temp),]
  
  
  P[j]=pars[1];
  K[j]=pars[2];
  Lambda[j]=pars[3];
  New[j]=pars[4];
  Sigma[j]=pars[5];
  j
}

mean(P);	sd(P);
mean(K);	sd(K);
mean(Lambda);	sd(Lambda);
mean(New);	sd(New);
mean(Sigma);	sd(Sigma);


sigs=rep(0,100);
nus=rep(0,100);
xmax=rep(0,100);
for (i in 1:100) {
  x=abs(rt(500,10)*sqrt(5));
  xmax[i]=max(x);
  g=function(x,nu,sig) {
    if (nu<=0 || sig <=0 || nu>=12) {
      return(10000000000);
    } else {
      result=-sum(log(2*exp(lgamma((nu+1)/2)-lgamma(nu/2))/sqrt(pi*nu*sig)*(1+x^2/(nu*sig))^(-(nu+1)/2)));
    }
    return(result);
  }
  
  pars=optim(c(10,5),g)$par
  nus[i]=pars[1];
  sigs[i]=pars[2]
}





##############################################################
#my overall likelihood########################################

start_time=Sys.time()

M=2;
my_par.in=foreach(j=1:length(links)) %do% {
  my_par.in=matrix(0,nrow=M,ncol=7);
}
my_par.out=foreach(j=1:length(links)) %do% {
  my_par.out=matrix(0,nrow=M,ncol=7);
}

my_llik.in=foreach(j=1:length(links)) %do% {
  my_llik.in=rep(0,M);
}

my_llik.out=foreach(j=1:length(links)) %do% {
  my_llik.out=rep(0,M);
}



for (m in 1:M) {
  
  if (my_llik.in[[1]][M]!=0) {
    break;
  }
  
  
  
  
  R.in<-foreach(j=1:length(links)) %do% {
    R.in=spaces.in[[j]][spaces.in[[j]]>m];
  }
  
  R.out<-foreach(j=1:length(links)) %do% {
    R.out=spaces.out[[j]][spaces.out[[j]]>m];
  }
  
  
  A.in<-foreach(j=1:length(links)) %do% {
    A=c();
    nzeros=0;
    counter=1;
    n=1;
    while (n <=length(spaces.in[[j]]) ) {
      while (spaces.in[[j]][n]<=m & n<=length(spaces.in[[j]])) {
        counter=counter+spaces.in[[j]][n];
        if (spaces.in[[j]][n]>1) {
          nzeros=nzeros+spaces.in[[j]][n]-1;
        }
        n=n+1;
      }
      A=c(A,counter);
      counter=1;
      nzeros=0;
      n=n+1;
    }
    A=A[2:length(A)];
    A.in=A;
  }
  
  A.out<-foreach(j=1:length(links)) %do% {
    A=c();
    nzeros=0;
    counter=1;
    n=1;
    while (n <=length(spaces.out[[j]]) ) {
      while (spaces.out[[j]][n]<=m & n<=length(spaces.out[[j]])) {
        counter=counter+spaces.out[[j]][n];
        if (spaces.out[[j]][n]>1) {
          nzeros=nzeros+spaces.out[[j]][n]-1;
        }
        n=n+1;
      }
      A=c(A,counter);
      counter=1;
      nzeros=0;
      n=n+1;
    }
    A=A[2:length(A)];
    A.out=A;
  }
  
  
  foreach(j=1:length(links)) %do% {
    ####extract zeros and 1's
    A=A.in[[j]];
    R=R.in[[j]]-m;
    a=max(A);
    X=matrix(2,ncol=length(A),nrow=a);
    S=R[1];
    for (i in 1:(length(A)-1)) {
      X[1:A[i],i]=data.in[[j]][(S+1):(S+A[i])];
      S=S+R[i+1]+A[i]-1;
    }
    
    if(length(R)>length(A)) {
      X[1:A[length(A)],length(A)]=data.in[[j]][(S+1):(S+A[length(A)])];
    } else {
      X[1:(length(data.in[[j]])-S),length(A)]=data.in[[j]][(S+1):length(data.in[[j]])];
    }
    
    
    
    
    #####calculate log likelihood
    inv.logit<-function(x) {
      result=1/(1+exp(-x));
      return(result);
    }
    
    logistic_log_lik<-function(b) {
      result=0;
      for (i in 1:length(A)) {
        if (A[i]>2) {
          result=result-sum(dbinom(X[2:(A[i]-1),i],1,inv.logit(b[1]+b[2]*(2:(A[i]-1))+b[3]*A[i]),log=TRUE));
        }
      }
      return(result);
    }
    
    
    
    my_log_likelihood=function(x) {
      if(x[1]<=0 || x[2]<1 || x[3]<=0 || x[3]>=min(R) || x[4]<=0) {
        return(10^(10));
      } else {
        
        result=-sum(log(x[1]*x[2]*(log(R/x[3]))^(x[2]-1)*exp(-x[1]*(log(R/x[3]))^x[2])/R));
        result=result-sum(log(x[4]/(A^(x[4]+1))));
        if (m>2) {
          result=result+logistic_log_lik(x[5:7]);
        }
        return(result);
      }
    }
    temp_par=matrix(0,nrow=200,ncol=7);
    temp=rep(0,200);
    for(k in 1:200) {
      u=c(runif(4,0,5),0,0,0);
      
      bs=my_log_likelihood(u);
      if (bs==Inf || is.nan(bs)==TRUE) {
        temp[k]=-100000;
      } else {
        fit_my=optim(u,my_log_likelihood);
        temp[k]=-fit_my$value;
        temp_par[k,]=fit_my$par
      }
    }
    my_llik.in[[j]][m]=max(temp);
    my_par.in[[j]][m,]=temp_par[which.max(temp),];
    
    
    
    
    ################################3outgoing
    
    
    A=A.out[[j]];
    R=R.out[[j]]-m;
    a=max(A);
    X=matrix(2,ncol=length(A),nrow=a);
    S=R[1];
    for (i in 1:(length(A)-1)) {
      X[1:A[i],i]=data.out[[j]][(S+1):(S+A[i])];
      S=S+R[i+1]+A[i]-1;
    }
    
    if(length(R)>length(A)) {
      X[1:A[length(A)],length(A)]=data.out[[j]][(S+1):(S+A[length(A)])];
    } else {
      X[1:(length(data.out[[j]])-S),length(A)]=data.out[[j]][(S+1):length(data.out[[j]])];
    }
    
    
    #####calculate log likelihood
    inv.logit<-function(x) {
      result=1/(1+exp(-x));
      return(result);
    }
    
    logistic_log_lik<-function(b) {
      result=0;
      for (i in 1:length(A)) {
        if (A[i]>2) {
          result=result-sum(dbinom(X[2:(A[i]-1),i],1,inv.logit(b[1]+b[2]*(2:(A[i]-1))+b[3]*A[i]),log=TRUE));
        }
      }
      return(result);
    }
    
    
    ###x[1]: lambda ###x[2]: nu ###x[3]: sigma ###x[4]: beta
    my_log_likelihood=function(x) {
      if(x[1]<=0 || x[2]<1 || x[3]<=0 || x[3]>=min(R) || x[4]<=0) {
        return(10^(10));
      } else {
        
        result=-sum(log(x[1]*x[2]*(log(R/x[3]))^(x[2]-1)*exp(-x[1]*(log(R/x[3]))^x[2])/R));
        result=result-sum(log(x[4]/(A^(x[4]+1))));
        if (m>2) {
          result=result+logistic_log_lik(x[5:7]);
        }
        return(result);
      }
    }
    
    temp_par=matrix(0,nrow=200,ncol=7);
    temp=rep(0,200);
    for(k in 1:200) {
      u=c(runif(4,0,5),0,0,0);
      
      bs=my_log_likelihood(u);
      if (bs==Inf || is.nan(bs)==TRUE) {
        temp[k]=-100000;
      } else {
        fit_my=optim(u,my_log_likelihood);
        temp[k]=-fit_my$value;
        temp_par[k,]=fit_my$par
      }
    }
    my_llik.out[[j]][m]=max(temp);
    my_par.out[[j]][m,]=temp_par[which.max(temp),];
    
  }
  
}
stop_time=Sys.time();
stop_time-start_time;


for (j in 1:length(links)) {
  st<-paste(c('loglik',links[j]),collapse='_');
  
  jpeg(paste(c(st,'.jpg'),collapse=''));
  plot(my_llik.in[[j]],type='b',pch='*',xlab="M",ylab="Log-Likelihood",main=paste(c("Log-Likelihood of",links[j],"Data"),collapse=' '),
       ylim=c(min(c(my_llik.in[[j]],my_llik.out[[j]]))-200,max(c(my_llik.in[[j]],my_llik.out[[j]]))+200));
  lines(my_llik.out[[j]],col='blue',type='b',pch='+',xlab="M",ylab="Log-Likelihood",main=paste(c("Log-Likelihood of Outgoing",links[j],"Data"),collapse=' '));
  legend("bottomleft",legend=c("Incoming","Outgoing"),col=c("black","blue"),pch=c('*','+'),lty=c(5,5))
  
  dev.off();
}

......................................................................................

######convolution bullshit
logf=function(pi,p,nu,sig,lam,S) {
  #pi=.4; p=.2; nu=2;lam=2;
  N=max(S);
  x=seq(1.5,N+.5,1);
  F_pps=pPPS(x,lam,sig,nu);
  f_par=mgeom(x-.5,pi,p);
  f_pps=rep(0,N);
  f_pps[1]=F_pps[1]; f_pps[2:N]=F_pps[2:N]-F_pps[1:(N-1)];

  C=convolve(f_pps,rev(f_par),type="o");
  result=sum(log(abs(C[S])));
  return(result);
}

l=length(links);
conv.lik.in=rep(0,l);
conv.lik.out=rep(0,l);
for (i in 1:l) {
  if (length(A.in[[i]])==length(R.in[[i]])) {
    S=A.in[[i]]+R.in[[i]]-1;
  } else {
    S=A.in[[i]]+R.in[[i]][-1]-1;
  }
  p=my_par.in[[i]][m,1:5];
  conv.lik.in[i]=logf(p[4],p[5],p[3],p[2],p[1],S);
  if (length(A.out[[i]])==length(R.out[[i]])) {
    S=A.out[[i]]+R.out[[i]]-1;
  } else {
    S=A.out[[i]]+R.out[[i]][-1]-1;
  }
  p=my_par.out[[i]][m,1:4];
  conv.lik.out[i]=logf(p[4],p[5],p[3],p[2],p[1],S);
}



hieu.lik.in=rep(0,l);
hieu.lik.out=rep(0,l);

for (i in 1:l) {
  
  S=A.in[[i]]+R.in[[i]][-1]-1;
  
  hieu_log_likelihood=function(x) {
    
    if(x[1]<0 || x[1]>1 || x[2]<=0 || x[3]<=0 || x[4]<=1 || x[5]<=0) {
      return(10^(20));
    } else {
      
      result=-sum(log(x[1]*(x[2]/x[3])*(S/x[3])^(x[2]-1)*exp(-(S/x[3])^x[2])+(1-x[1])*2*gamma((x[4]+1)/2)/(gamma(x[4]/2)*sqrt(x[4]*pi*x[5]))*(1+S^2/(x[4]*x[5]))^(-(x[4]+1)/2)));
      
      return(result);
      
    }
  }
  
  temp_par=matrix(0,nrow=100,ncol=5);
  temp=rep(0,100);
  for(k in 1:100) {
    u1=runif(2,0,1);
    u=runif(3,0,10);
    u=c(u1[1],u[1],u[2],u1[2],u[3]);
    bs=hieu_log_likelihood(u);
    if (bs!=Inf || is.nan(bs)==FALSE) {
      fit_hieu=optim(u,hieu_log_likelihood);
      temp[k]=-fit_hieu$value;
      temp_par[k,]=fit_hieu$par
    } else {
      temp[k]=-100000;
    }
  }
  hieu.lik.in[i]=max(temp);
  
  
  
  S=A.out[[i]]+R.out[[i]][-1]-1;
  
  hieu_log_likelihood=function(x) {
    
    if(x[1]<0 || x[1]>1 || x[2]<=0 || x[3]<=0 || x[4]<=1 || x[5]<=0) {
      return(10^(20));
    } else {
      
      result=-sum(log(x[1]*(x[2]/x[3])*(S/x[3])^(x[2]-1)*exp(-(S/x[3])^x[2])+(1-x[1])*2*gamma((x[4]+1)/2)/(gamma(x[4]/2)*sqrt(x[4]*pi*x[5]))*(1+S^2/(x[4]*x[5]))^(-(x[4]+1)/2)));
      
      return(result);
      
    }
  }
  
  temp_par=matrix(0,nrow=100,ncol=5);
  temp=rep(0,100);
  for(k in 1:100) {
    u1=runif(2,0,1);
    u=runif(3,0,10);
    u=c(u1[1],u[1],u[2],u1[2],u[3]);
    bs=hieu_log_likelihood(u);
    if (bs!=Inf || is.nan(bs)==FALSE) {
      fit_hieu=optim(u,hieu_log_likelihood);
      temp[k]=-fit_hieu$value;
      temp_par[k,]=fit_hieu$par
    } else {
      temp[k]=-100000;
    }
  }
  hieu.lik.out[i]=max(temp);
  
  
}


par(mfrow=c(1,2));


A=matrix(c(links,round(conv.lik.in,0),round(hieu.lik.in,0)),nrow=l,ncol=3,byrow=FALSE);
barplot(A,beside=TRUE);


A=matrix(c(conv.lik.in,hieu.lik.in),nrow=2,ncol=l,byrow=TRUE);
barplot(A,beside=TRUE,col=c('black','blue'));

A=matrix(c(conv.lik.out,hieu.lik.out),nrow=2,ncol=l,byrow=TRUE);
barplot(A,beside=TRUE,col=c('black','blue'));


plot(conv.lik.in-hieu.lik.in)


Z=data.frame(Links=links,
             "Our Likelihood"=A[1,],
             "Hieu Likelihood"=A[2,])
xtable(Z,row.names=FALSE)


