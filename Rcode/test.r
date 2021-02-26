library(LindleyR)
#############  generate the Bin-Lindley random numbers
Gene<-function(m,pi,beta,n){
  m<-m
  n<-n
  #beta<-2
  #pi<-0.3
  theta<-1/beta
  alpha<- 1/(beta*pi)-1/beta
  ########################  Generate the Bin-Lindley data
  Z<-numeric(n)
  X<-numeric(n)
  Y<-numeric(n)
  M<-numeric(n)
  W<-matrix(0,nrow=n,ncol=2)
  for(i in 1:n){
    Z[i]<- rslindley(1, theta, alpha, mixture = TRUE)
    X[i]<- rbinom(1, m, exp(-Z[i]))
    Y[i]<- lamb(m, pi,beta,X[i])#*beta/(0.5025*beta+0.6643)
    M[i]<-m
    W[i,]<-c(X[i],M[i])
  }
  return(W)
}

####      Compute the conditional expectation of Lindley R.V given Bin-Lindley R.V.
###################################################################
lamb=function(m, pi, beta,x){
  k<-m-x
  num<-numeric(1)
  deno<-numeric(1)
  lamb<-numeric(1)
  for (i in 0:k){
    num <-num+ choose(k,i)*(-1)**i*(2*beta-pi*beta+(x+i)*pi*beta**2)/((1+(x+i)*beta)**3)
    deno<-deno+choose(k,i)*(-1)**i*(1+(x+i)*pi*beta)/((1+(x+i)*beta)**2)
  }
  
  lamb=num/(deno)
  return(lamb)
}

##############  conditional expectation of Z given X
zed=function(m, pi, beta,x){
  k<-m-x
  num<-numeric(1)
  deno<-numeric(1)
  zed<-numeric(1)
  for (i in 0:k){
    num <-num+ choose(k,i)*(-1)**i*pi/(1+(x+i)*beta)
    deno<-deno+choose(k,i)*(-1)**i*(1+(x+i)*pi*beta)/((1+(x+i)*beta)**2)
  }
  
  zed=num/(deno)
  return(zed)
}

################  COmpute the MLE for Bin-Lindley Model  #################
###########################################################################
Est5 = function(Y){
  ## This function is to calculate the MLEs in the Lindley model
  ## This EM algorithm is my method
  n <- length(Y[,1])
  p=1
  X<-numeric(n)
  M<-numeric(n)
  Z<-numeric(n)
  W<-numeric(n)
  X<-Y[,1]
  M<-Y[,2]
 
  pi <- numeric(p)
  pi0<- numeric(p)
  beta <- numeric(p)
  beta0<- numeric(p)
  pi0<-0.5
  #beta0<-1.0
  beta0<-mean(M)/mean(X)-1
  #beta0
  pi.new <- numeric(p)
  beta.new<-numeric(p)
  beta1.new<-numeric(p)
  pi <- numeric(p)
  beta<-numeric(p)
  
  for(i in 1:n){
    Z[i]=   zed(M[i],pi0,beta0,X[i])
    W[i]=lamb(M[i],pi0,beta0,X[i])
      pi.new<-  pi.new+Z[i]/n
    beta1.new<-beta1.new+W[i]
  }
    beta.new<-beta1.new/(2*n-n*pi.new)
 #pi.new
 #beta.new
 #mean(Z)
 #mean(W)
  repe <- 1 
  while ((abs(pi.new-pi)>1e-5)|(abs(beta.new-beta)>1e-5)){
    beta <- beta.new; pi <- pi.new; repe <- repe+1
    pi.new <- numeric(p)
    beta1.new<-numeric(p)
    beta.new<-numeric(p)
    for(i in 1:n){
      W[i]=lamb(M[i],pi,beta,X[i])
      Z[i]=zed(M[i],pi,beta,X[i])
      pi.new<-pi.new+Z[i]/n
      beta1.new<-beta1.new+W[i]
    }
    beta.new<-beta1.new/(2*n-n*pi.new)
    #beta.new
    #beta
    #pi.new
    #pi
    #mean(Z)
    #mean(W)
    #m*(1+pi.new*beta.new)/(1+beta.new)^2
    #pi.new*beta.new+2*(1-pi.new)*beta.new
  }
  
  para<-c(pi.new,beta.new)
  return(para)
}

MZIB.CI <- function(Y,alpha,G){
  ## i-th row is for the i-th parameter
  ## The 1st is the lower bound; 2nd is the upper bound, 3rd is the width;
  ##  4th is the mean, and 5th is the median
  n<- dim(Y)[1]; p <- dim(Y)[2]
  #n<- dim(ZMM1)[1]; p <- dim(ZMM1)[2]
  M <- matrix(0,nrow=n,ncol=1)
  M<-Y[,2]
  # m<-ZMM1[,(s+1):p]
  #m<- Y[, (s+1):p]
  e <-  matrix(0,nrow=G,ncol=2)
  ee <- matrix(0,nrow=G,ncol=2)
  e2 <- matrix(0,nrow=2,ncol=6)
  a <- Est5(Y)
  pi.est <- a[1]
  beta.est <- a[2] 
  
  for (i in 1:G){
    YY <- matrix(0,nrow=n,ncol=2)
    for (k in 1:n){
      YY[k,] <- Gene(M[k],pi.est,beta.est,1)
    }
    e[i,] <- Est5(YY)      
  }
  q1 = round(G*alpha/2); q2 = G-q1+1
  for (j in 1:2){
    ee <- sort(e[,j])
    
    aa=0
    for (k in 1:G){
      if (ee[k]==0) aa=aa+1
    }

    e2[j,1] <- ee[q1];  e2[j,2] <- ee[q2]
    e2[j,3] <- ee[q2]-ee[q1]; e2[j,4] <- mean(ee)
    e2[j,5] <- median(ee)
    e2[j,6] <- aa/G
  }

  print(e2)
  return(e2)
}

m=10
n <-100
G <-200
t <- 500
alpha <- 0.05
s <- length(pi)
count <- numeric(2)
wid <- numeric(2)
pi<-0.2
beta<-0.5
para<-c(pi,beta)
est1 <- matrix(0, ncol = 2, nrow= t)
time0 = proc.time()
for (i in 1:t){
  Y <- Gene(m, pi, beta, n)   
  est1[i,] <- Est5(Y)
  b <- MZIB.CI(Y,alpha,G) 
  for (j in 1:2){      
    if  ((b[j,1]<=para[j])&(b[j,2]>=para[j])) {count[j] <- count[j]+1}
    wid[j] <- wid[j]+b[j,3]
  }
  
}
C<-c(mean(est1[,1]),mean(est1[,2]),count[1],count[2],wid[1],wid[2])
proc.time()-time0
result<-A
result<-rbind(A,B)
proc.time()-time0
apply(est1,2,mean)