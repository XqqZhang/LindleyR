library(LindleyR)
library(nleqslv)

Dla<-function(p){
dla=-n/(p[2]+p[1])
for(i in 1:n){
    s=0
    t=0

    for(j in 0:x[i]){
        a=choose(x[i],j)*((-1)^j)*(p[2]+j+m-x[i]+p[1])/(p[2]+j+m-x[i])^2
        s=s+a
        j=j+1
    }
    for(j in 0:x[i]){
        b=choose(x[i],j)*((-1)^j)/(p[2]+j+m-x[i])^2 
        t=t+b
        j=j+1
    }
            
    c=t/s
    dla=dla+c
    i=i+1
} 
return(dla)
}
         
Dlt<-function(p){
   dlt=2*n/p[2]-n/(p[2]+p[1])
for(i in 1:n){
    s=0
    t=0

    for(j in 0:x[i]){
        a=choose(x[i],j)*((-1)^j)*(p[2]+j+m-x[i]+p[1])/(p[2]+j+m-x[i])^2
        s=s+a
        j=j+1
    }
    for(j in 0:x[i]){
        b=choose(x[i],j)*((-1)^(j+1))*(p[2]+j+m-x[i]+2*p[1])/(p[2]+j+m-x[i])^3 
        t=t+b
        j=j+1
    }
            
    c=t/s
    dlt=dlt+c
    i=i+1
}
return(dlt) 
}

model<- function(p){
   a=Dla(p)
   b=Dlt(p)
   c(a=a,b=b)
}

m=12
n=1000
B=2000
alpha=5
theta=5
esta<-matrix(0, nrow=B, ncol = 2)
for(b in 1:B){
     for(i in 1:n){
     y[i]<- rslindley(1, theta, alpha, mixture = TRUE)
     x[i]<- rbinom(1, m, 1-exp(-y[i]))
     }
    xstart  <-  c(4.8,5.3)
    esta[b,]= nleqslv(xstart,model)$x
    
}
estb<-c(mean(esta[,1]),mean(esta[,2]))