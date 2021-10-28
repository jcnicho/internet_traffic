##this makes the stupid fucking plot package in R from the PPS
par(mfrow=c(2,2),mar=c(2,3,1,1))
x=sort(R-m);
est=fit$estimate;
t=seq(0,3000,1);
y=dPPS(t,est$lambda,est$sigma,est$nu);
par(mfrow=c(2,2));
hist(x,freq=FALSE,ylab="",main="",xlab="")

lines(t,y,lty=2,lwd=2,col='red');

Fn=ecdf(x)
plot(R-m,Fn(R-m),ylab="",xlab="");
y=pPPS(x,est$lambda,est$sigma,est$nu);
lines(x,y,lty=2,lwd=2,col='red');

plot(log(x),log(rank(1-Fn(x))),ylab="",xlab="");
lines(log(x),log(rank(1-y)),lty=2,lwd=2,col='red');

plot(log(log(x/est$sigma)),log(-log(1-Fn(x))),ylab="",xlab="");
lines(log(log(x/est$sigma)),log(-log(1-y)),lty=2,lwd=2,col='red');

#####Active Periods##########################################

A=A.in[[1]]
par(mfrow=c(2,2),mar=c(2,3,1,1))
x=sort(A);

log.mg=function(x) {
  result=0;
  for (i in 1:length(A)) {
    result=result-log(mgeom(A[i],x[1],x[2]));
  }
  return(result);
}

fit=optim(c(.5,.5),log.mg)$par;
t=seq(1,max(A)+1,1);
y=mgeom(t,fit[1],fit[2]);

par(mfrow=c(2,2));
hist(x,freq=FALSE,breaks=seq(0,max(A),1),ylab="",main="",xlab="")

lines((t-.5),y,lty=2,lwd=2,col='red');

Fn=ecdf(x)
plot(A,Fn(A),ylab="",xlab="");
y=pmgeom(x,fit[1],fit[2]);
lines(x,y,lty=2,lwd=2,col='red');

plot(log(x),log(rank(1-Fn(x))),ylab="",xlab="");
lines(log(x),log(rank(1-y)),lty=2,lwd=2,col='red');

plot(log(log(x/est$sigma)),log(-log(1-Fn(x))),ylab="",xlab="");
lines(log(log(x/est$sigma)),log(-log(1-y)),lty=2,lwd=2,col='red');




