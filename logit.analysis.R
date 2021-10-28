###logistic regression analysis
logistic_log_lik_red=function (x) {
  return(logistic_log_lik(c(x[1],0,x[2])));
}
M=30; 
for (m in 1:M) {
  
}
m=5
R.in<-foreach(j=1:length(links)) %do% {
  R.in=spaces.in[[j]][spaces.in[[j]]>m];
}

R.out<-foreach(j=1:length(links)) %do% {
  R.out=spaces.out[[j]][spaces.out[[j]]>m];
}


A.in<-foreach(j=1:length(links)) %do% {
  j=1
  A=c();
  nzeros=0;
  props=c();
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
    props=c(props,1-nzeros/counter);
    counter=1;
    nzeros=0;
    n=n+1;
  }
  A=A[2:length(A)];
  A.in=A; props=props[-1];
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

A=A.in[[l]];
R=R.in[[l]];
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
summary(fit);


llik=rep(0,M);
llik_red=llik;

for (m in 1:M) {
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
  A=A.in[[l]];
  R=R.in[[l]];
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
  llik[m]=-optim(c(0,0,0),logistic_log_lik)$value;
  llik_red[m]=-optim(c(0,0),logistic_log_lik_red)$value;
}

par(mfrow=c(1,2),mar=c(3,2,1,1));

R=sort(R.in[[1]]-m);
hist(R,breaks=seq(0,max(R)+50,50),freq=F,main="")

fit=PPS.fit(R)$estimate;
lines(R,dPPS(R,fit$lambda,fit$sigma,fit$nu),lwd=2);


zeta.f=function(x) {
  result=-sum(dzeta(R,x,log=TRUE));
  return(result);
}

fit=optim(1,zeta.f,method="Brent",lower=0,upper=10)$par

lines(R,dzeta(R,fit),lwd=2,lty=2);

A=A.in[[1]];
hist(A[A<50],breaks=seq(0,min(max(A),50),1),freq=FALSE,main="");
fit=optim(c(.5,.5),log.mg)$par;
t=seq(1,max(A)+1,1);
y=mgeom(t,fit[1],fit[2]);
lines((t-.5),y,lwd=2);

zeta.f=function(x) {
  result=-sum(dzeta(A,x,log=TRUE));
  return(result);
}

fit=optim(1,zeta.f,method="Brent",lower=0,upper=10)$par
lines((t-.5),dzeta(t,fit),lwd=2,lty=2);



