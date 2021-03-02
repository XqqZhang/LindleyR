library(LindleyR)
##Simulation
#y<- rslindley(100, 5, 5, mixture = TRUE)
#x<- rbinom(100, 12, 1-exp(-y))
set.seed(1)
##MLE    
Lx<- function(alpha,theta,m,n,x){
    
    lx=0
    for(i in 1:n){
            s=0

            for(j in 0:x[i]){
                a=choose(x[i],j)*((-1)^j)*(theta+j+m-x[i]+alpha)/(theta+j+m-x[i])^2
                s=s+a
                j=j+1
            }
            
            l=log(choose(m,x[i]))+2*log(theta)-log(theta+alpha)+log(s)
            lx=lx+l
            i=i+1
        }  
        return(lx)
}

Dla<- function(alpha,theta,m,n,x){
    dla=-n/(theta+alpha)
    for(i in 1:n){
            s=0
            t=0

            for(j in 0:x[i]){
                a=choose(x[i],j)*((-1)^j)*(theta+j+m-x[i]+alpha)/(theta+j+m-x[i])^2
                s=s+a
                j=j+1
            }
            for(j in 0:x[i]){
                b=choose(x[i],j)*((-1)^j)/(theta+j+m-x[i])^2 
                t=t+b
                j=j+1
            }
            
            c=t/s
            dla=dla+c
            i=i+1
        } 
         
        return(dla)
}

Ddla<- function(alpha,theta,m,n,x){
    ddla=(n^2)/(theta+alpha)^2
    for(i in 1:n){
            s=0
            t=0

            for(j in 0:x[i]){
                a=choose(x[i],j)*((-1)^j)*(theta+j+m-x[i]+alpha)/(theta+j+m-x[i])^2
                s=s+a
                j=j+1
            }
            for(j in 0:x[i]){
                b=choose(x[i],j)*((-1)^j)/(theta+j+m-x[i])^2 
                t=t+b
                j=j+1
            }
            
            c=(t/s)^2
            ddla=ddla-c
            i=i+1
        } 
         
        return(ddla)
}

Dlt<- function(alpha,theta,m,n,x){
    dlt=2*n/theta-n/(theta+alpha)
    for(i in 1:n){
            s=0
            t=0

            for(j in 0:x[i]){
                a=choose(x[i],j)*((-1)^j)*(theta+j+m-x[i]+alpha)/(theta+j+m-x[i])^2
                s=s+a
                j=j+1
            }
            for(j in 0:x[i]){
                b=choose(x[i],j)*((-1)^(j+1))*(theta+j+m-x[i]+2*alpha)/(theta+j+m-x[i])^3 
                t=t+b
                j=j+1
            }
            
            c=t/s
            dlt=dlt+c
            i=i+1
        } 
         
        return(dlt)
}

Ddlt<- function(alpha,theta,m,n,x){
    ddlt=-2*n/theta^2+n/(theta+alpha)^2
    for(i in 1:n){
            s=0
            t=0
            u=0  
            for(j in 0:x[i]){
                a=choose(x[i],j)*((-1)^j)*(theta+j+m-x[i]+alpha)/(theta+j+m-x[i])^2
                s=s+a
                j=j+1
            }
            for(j in 0:x[i]){
                b=choose(x[i],j)*((-1)^(j+2))*(theta+j+m-x[i]+3*alpha)/(theta+j+m-x[i])^4 
                t=t+b
                j=j+1
            }
            for(j in 0:x[i]){
                c=choose(x[i],j)*((-1)^(j+1))*(theta+j+m-x[i]+2*alpha)/(theta+j+m-x[i])^3 
                u=u+c
                j=j+1
            }
            
            c=2*t/s-(u/s)^2
            ddlt=ddlt+c
            i=i+1
        } 
         
        return(ddlt)
}

Ddlat<- function(alpha,theta,m,n,x){
    ddlat=n/(theta+alpha)^2
    for(i in 1:n){
            s=0
            t=0
            u=0
            v=0  
            for(j in 0:x[i]){
                a=choose(x[i],j)*((-1)^j)*(theta+j+m-x[i]+alpha)/(theta+j+m-x[i])^2
                s=s+a
                j=j+1
            }
            for(j in 0:x[i]){
                b=choose(x[i],j)*((-1)^(j+1))*2/(theta+j+m-x[i])^3 
                t=t+b
                j=j+1
            }
            for(j in 0:x[i]){
                c=choose(x[i],j)*((-1)^j)/(theta+j+m-x[i])^2 
                u=u+c
                j=j+1
            }
            for(j in 0:x[i]){
                d=choose(x[i],j)*((-1)^(j+1))*(theta+j+m-x[i]+2*alpha)/(theta+j+m-x[i])^3 
                v=v+d
                j=j+1
            }
            c=t/s-u*v/(s^2)
            ddlat=ddlat+c
            i=i+1
        } 
         
        return(ddlat)
}

Ans<- function(alpha,theta,m,n,x){
  dla=Dla(alpha,theta,m,n,x)
  dlt=Dlt(alpha,theta,m,n,x)
  ddla=Ddla(alpha,theta,m,n,x)
  ddlat=Ddlat(alpha,theta,m,n,x)
  ddlt=Ddlt(alpha,theta,m,n,x)
  u=c(dla,dlt)
  Inn=cbind(c(ddla,ddlat),c(ddlat,ddlt))
  inverse=solve(Inn)
  pt<-c(4.9,5.4)
  pt1=pt-inverse%*%u
  while (abs(pt1[1]-pt[1])>0.0001 | abs(pt1[2]-pt[2])>0.0001) {
    pt[1]=pt1[1]
    pt[2]=pt1[2]
    pt1=pt-inverse%*%u
  }
  return(pt1)  
}

alpha=5
theta=5
m=12
n=100
B=2
Z<-matrix(0, nrow=B, ncol = 2)
for (b in 1:B){
     for(i in 1:n){
     y[i]<- rslindley(1, theta, alpha, mixture = TRUE)
     x[i]<- rbinom(1, m, 1-exp(-y[i]))
     }
    Z[b,]<- Ans(alpha,theta,m,n,x)
}
est<-c(mean(Z[,1]),mean(Z[,2]))






