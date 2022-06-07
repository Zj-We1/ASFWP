library(R.matlab,warn.conflicts = FALSE)
library(pracma)
library(igraph)
w.norm=function(y,w,beta){
  return(sum(y^2*w^(beta)))
}
ASFWP=function(X,s=-1,k,lambda,eta=1.04,beta=1,tmax=200,tol=10^(-3)){
  n=dim(X)[1]
  p=dim(X)[2]
  phi=matrix(0,nrow=n,ncol=k)
  m=matrix(0,nrow=n,ncol=k)
  m1=matrix(0,nrow=n,ncol=k)
  label=numeric(n)
  dd=numeric(k)
  #initialization
  sam=sample(n,k)
  theta=X[sam,]
  for(i in 1:p){
    theta[,i]=theta[,i]+rnorm(k,0,0.1)
  }
  W=rep(1/p,p)
  W1=rep(1/p,p)
  D=rep(0,p)
  D1=rep(0,n)
  M=rep(0,p)
  h=rep(1/n,n)
  val1=theta
  con=rep(0,tmax)
  for(t in 1:tmax){
    #update membership
    for(i in 1:n){
      for(l in 1:k){
        m[i,l]=h[i]*w.norm((X[i,]-theta[l,]),W,beta)
      }
    }
    for(i in 1:n){
      for(l in 1:k){
        m1[i,l]=w.norm((X[i,]-theta[l,]),W,beta) 
      }
    }
    for(i in 1:n){
      for(l in 1:k){
        phi[i,l]=((m[i,l]+0.0001)^(s-1))*(sum(m[i,]^s)+0.0001)^((1/s)-1)
      }
    }
    #update centroids
    for(l in 1:k){
      for(d in 1:p){
        theta[l,d]=(sum(phi[,l]*X[,d]))/sum(phi[,l])
      }
    }
    #update s
    if(t%%2==0){
      s=eta*s
    }
    #update feature weights
    for(l in 1:p){
      s1=0
      for(i in 1:n){
        for(j in 1:k){
          s1=s1+phi[i,j]*h[i]*(X[i,l]-theta[j,l])^2
        }
      }
      tmp=-s1/lambda
      D[l]=exp(tmp)
    } 
    W=D/sum(D)
    cat('\n')
    #update sample weights
    for(i in 1:n)
    {
      s2=0
      for(j in 1:k)
      {
        s2=s2+phi[i,j]*m1[i,j]
      }
      D1[i]=1/s2
    }
    h=D1/sum(D1)
    cat(t)
    cat('\n')
    t=t+1
    val2=theta
    a=norm(val1-val2)
    cat(a)
    cat('\n')
    con[t]=a
    a[is.na(a)]<-0
    cat('\n')
    if(a<tol){
      break
    }
    else{
      val1=val2
    }
  }
  
  for(i in 1:n){
    for(j in 1:k){
      dd[j]=sum((X[i,]-theta[j,])^2*W)
    }
    label[i]=which.min(dd)
  }
  list2=list(theta,label,W,h,con)
  names(list2)=c('theta','label','f_weight',"s_weight",'con')
  return(list2)
}
data_generate=function(n,M,prob1,sigma,sigma2){
  p=dim(M)[2]
  X=matrix(0,n,p)
  k=dim(M)[1]
  label=numeric(n)
  for(i in 1:n){
    s=sample(1:k,size=1,prob=prob1)
    for(l in 1:p){
      if(M[s,l]==0){
        X[i,l]=rnorm(1,M[s,l],sigma2)
      }else{
        X[i,l]=rnorm(1,M[s,l],sigma)
      }
      
    }
    label[i]=s
  }
  ls=list(X,label)
  names(ls)=c('data','label')
  return(ls)
}



k=5
p=20
M=rand(k,p)
prop=0.95
s=sample(p,floor(p*prop))
s=6:20
M[,s]=0
X=data_generate(1000,M,rep(1/k,k),0.02,1)
label=X$label
X=X$data
n=dim(X)[1]
p=dim(X)[2]
a=200
for(i in 1:a)
{
  X[i,]=X[i,]+rnorm(p,mean=0,sd=0.5)
}
for(i in 1:p){
  X[,i]=(X[,i]-mean(X[,i]))/sd(X[,i])
}
l0=ASFWP(X,k=5,s=-1,eta=1,lambda=1.03,tmax=100,tol=0.001)
compare(label,l0$label,'nmi')
compare(label,l0$label,'adjusted.rand')
a1=rep(sum(l0$s_weight[1:200])/200,200)
a2=rep(sum(l0$s_weight[201:1000])/800,800)
plot(1:n,l0$s_weight,type="l",xlab='sample',ylab='s_weight')
lines(1:200,a1,col=2,lwd=5)
lines(201:1000,a2,col=4,lwd=5)
plot(1:p,l0$f_weight,type="b",xlab="feature",ylab="f_weight")

