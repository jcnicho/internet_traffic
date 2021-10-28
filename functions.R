
##The purpose of this section is to transform the data into something that we can deal with
##more easily. First, we transform the data into binary sequences for each direction in each
##link. Then, we calculate the "spaces", meaning the number of data points between one's, noting
##that consecutive one's yields the space of 1. Using the spaces data, we can easily calculate 
##the active and inactive periods for arbitrary choice of m.
data.in<-foreach(j=1:length(links)) %do% {
  data.in=as.numeric(FourierAnom[[j]][,2]!=0); #transform to binary
}
data.out<-foreach(j=1:length(links)) %do% {
  data.out=as.numeric(FourierAnom[[j]][,3]!=0);
}



spaces.in<-foreach(j=1:length(links)) %do% { #transform to spaces
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


inv.logit<-function(x) { #inverse logit function for logistic likelihood
  result=1/(1+exp(-x));
  return(result);
}

logistic_log_lik<-function(b) { #calculates -log likelihood of logistic component
  result=0;
  for (i in 1:(length(A)-1)) {
    if (A[i]>2) {
      result=result-sum(dbinom(X[2:(A[i]-1),i],1,inv.logit(b[1]+b[2]*(2:(A[i]-1))+b[3]*A[i]),log=TRUE));
    }
  }
  return(result);
}


start_time=Sys.time()

M=30;
my_par.in=foreach(j=1:length(links)) %do% { #initialize vectors of parameters to be estimated
  my_par.in=matrix(0,nrow=M,ncol=8);
}
my_par.out=foreach(j=1:length(links)) %do% {
  my_par.out=matrix(0,nrow=M,ncol=8);
}

my_llik.in=foreach(j=1:length(links)) %do% {
  my_llik.in=rep(0,M);
}

my_llik.out=foreach(j=1:length(links)) %do% {
  my_llik.out=rep(0,M);
}



##Here we will maximize the log likelihood function and estimate parameters for each direction
##of each link, for M=1:30. Since the likelihood function can be broken into several different
##pieces such that different parameters do not depend on one another, it is sufficient to maximize
##the likelihood for the inactive periods, the active periods, and the strings of ones and zeros
##within the active periods. Thus, three different likelihoods are calculated
for (m in 1:M) {
  R.in<-foreach(j=1:length(links)) %do% { #inactive periods for given m
    R.in=spaces.in[[j]][spaces.in[[j]]>m];
  }
  
  R.out<-foreach(j=1:length(links)) %do% {
    R.out=spaces.out[[j]][spaces.out[[j]]>m];
  }
  
  
  A.in<-foreach(j=1:length(links)) %do% { #active periods for given m
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
  
  
  for(j in 1:length(links)) {
  
  
    ####incoming links
    ####extract zeros and 1's
    A=A.in[[j]];
    R=R.in[[j]];
    a=max(A);
    X=matrix(2,ncol=length(A),nrow=a); #a matrix of zeros and ones for logistic regression
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
    R=R-m;
    
    
    log.1=0; 
    if (m>=2) {
    logit.fit=optim(c(0,0,0),logistic_log_lik); #calculate logistic likelihood
    log.1=-logit.fit$value;
    my_par.in[[j]][m,6:8]=logit.fit$par;
    }
    
    #to guarantee stability of the pps likelihood, we establish lower bounds for parameters 
    log.pps=function(x) {
      if (x[2]<10e-5 || x[1]<=10e-7 || x[3] <=0) {
        return(10e10);
      }
      result=-sum(log(mPPS(R,x[1],x[2],x[3])));
      return(result);
    }
    pps.fit=optim(c(1,1,1),log.pps); #calculate likelihood of inactive periods
    log.2=-pps.fit$value;
    my_par.in[[j]][m,1:3]=pps.fit$par
    
    
    mgeom=function(x,pi,p) { #geometric mixture for active periods
      result=rep(0,length(x));
      for (i in 1:length(x)) {
        result[i]=pi*as.numeric(x[i]==1)+(1-pi)*p*(1-p)^(x[i]-1);
      }
      return(result);
    }
    
    log.mg=function(x) { #calculate likelihood of active periods
      result=0;
      if (x[1] < 0 || x[1]>1 || x[2]<0 || x[2]>1) {
        return(10^10);
      }
      for (i in 1:length(A)) {
        result=result-log(mgeom(A[i],x[1],x[2]));
      }
      return(result);
    }
    
    mg.fit=optim(c(.5,.5),log.mg); #maximize likelihood of active periods
    log.3=-mg.fit$value;
    my_par.in[[j]][m,4:5]=mg.fit$par;
    
    my_llik.in[[j]][m]=log.1+log.2+log.3;
    
    ####incoming links
    ####extract zeros and 1's
    A=A.out[[j]];
    R=R.out[[j]];
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
    R=R-m;
    
    
    log.1=0; 
    if (m>=2) {
      logit.fit=optim(c(0,0,0),logistic_log_lik);
      log.1=-logit.fit$value;
      my_par.in[[j]][m,6:8]=logit.fit$par;
    }
    
    
    log.pps=function(x) {
      if (x[2]<10e-5 || x[1]<=10e-7 || x[3] <=0) {
        return(10e10);
      }
      result=-sum(log(mPPS(R,x[1],x[2],x[3])));
      return(result);
    }
    pps.fit=optim(c(1,1,1),log.pps);
    log.2=-pps.fit$value;
    my_par.out[[j]][m,1:3]=pps.fit$par
    
    
    log.mg=function(x) {
      result=0;
      if (x[1] < 0 || x[1]>1 || x[2]<0 || x[2]>1) {
        return(10^10);
      }
      for (i in 1:length(A)) {
        result=result-log(mgeom(A[i],x[1],x[2]));
      }
      return(result);
    }
    
    mg.fit=optim(c(.5,.5),log.mg);
    log.3=-mg.fit$value;
    my_par.out[[j]][m,4:5]=mg.fit$par;
    
    my_llik.out[[j]][m]=log.1+log.2+log.3;
    
    
}
   
}

stop_time=Sys.time();
stop_time-start_time;

