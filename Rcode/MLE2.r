library(LindleyR)
##Simulation
set.seed(1)
y<- rslindley(100, 5, 5, mixture = TRUE)
x<- rbinom(100, 12, 1-exp(-y))

##MLE    
Lx<- function(alpha,theta,m,n){
    
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

Dla<- function(alpha,theta,m,n){
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

Ddla<- function(alpha,theta,m,n){
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

Dlt<- function(alpha,theta,m,n){
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

Ddlt<- function(alpha,theta,m,n){
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

Ddlat<- function(alpha,theta,m,n){
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

Pt<-function(alpha,theta){
    pt<-c(alpha,theta)
    return(pt)
}


In<- function(alpha,theta,m,n){
    ddla=Ddla(alpha,theta,m,n)
    ddlat=Ddlat(alpha,theta,m,n)
    ddlt=Ddlt(alpha,theta,m,n)
    Inn=cbind(c(ddla,ddlat),c(ddlat,ddlt))
    return(Inn)
}

Inverse<- function(alpha,theta,m,n){
    i=In(alpha,theta,m,n)
    inverse=solve(i)
    return(inverse)
}


U<- function(alpha,theta,m,n){
    dla=Dla(alpha,theta,m,n)
    dlt=Dlt(alpha,theta,m,n)
    u=c(dla,dlt)
    return(u)
}

Pt1<- function(alpha,theta,m,n){
  pt=Pt(alpha,theta)
  inverse=Inverse(alpha,theta,m,n)
  u=U(alpha,theta,m,n)   
  pt1=pt-inverse%*%u
  return(pt1)
}

Ans<- function(alpha,theta,m,n){
  pt=Pt(alpha,theta)
  pt1=Pt1(alpha,theta,m,n)
  while (abs(pt1[1]-pt[1])>0.00001 | abs(pt1[2]-pt[2])>0.00001) {
    pt[1]=pt1[1]
    pt[2]=pt1[2]
    pt1=Pt1(pt[1],pt[2],m,n) 
    print(pt1)
  }
  return(pt1)  
}






