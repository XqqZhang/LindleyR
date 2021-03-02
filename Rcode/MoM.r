library(LindleyR)
library(nleqslv)
#y<- rslindley(100, 5, 5, mixture = TRUE)
#x<- rbinom(100, 12, 1-exp(-y))
# model<- function(p){
#    me=mean(x)
#    me2=mean(x^2)
#    a=m*(p[2]^2+p[2]+2*p[1]*p[2]+p[1])/((p[2]+p[1])*(p[2]+1)^2)-me
#    b=m^2-((2*m^2-m)*(p[2]^2)*(p[2]+1+p[1]))/((p[2]+p[1])*(p[2]+1)^2)+((m^2-m)*(p[2]^2)*(p[2]+2+p[1]))/((p[2]+p[1])*(p[2]+2)^2)-me2
#    c(a=a,b=b)
# }
# m=12
# n=100
# B=1
# alpha=5
# theta=5
# esta<-matrix(0, nrow=B, ncol = 2)
# for(b in 1:B){
#      for(i in 1:n){
#      y[i]<- rslindley(1, theta, alpha, mixture = TRUE)
#      x[i]<- rbinom(1, m, 1-exp(-y[i]))
#      }
#     ss<- multiroot(f=model,start=c(4.7,5.3))
#     esta[b,]=ss$root
# }
# estb<-c(mean(esta[,1]),mean(esta[,2])) 
# estb

model<- function(p){
   me=mean(x)
   me2=mean(x^2)
   a=m*(p[2]^2+p[2]+2*p[1]*p[2]+p[1])/((p[2]+p[1])*(p[2]+1)^2)-me
   b=m^2-((2*m^2-m)*(p[2]^2)*(p[2]+1+p[1]))/((p[2]+p[1])*(p[2]+1)^2)+((m^2-m)*(p[2]^2)*(p[2]+2+p[1]))/((p[2]+p[1])*(p[2]+2)^2)-me2
   c(a=a,b=b)
}

m=12
n=2000
B=2000
alpha=5
theta=5
x<-numeric(n)
y<-numeric(n)
momEsta<-matrix(0, nrow=B, ncol = 2)
for(b in 1:B){
     for(i in 1:n){
     y[i]<- rslindley(1, theta, alpha, mixture = TRUE)
     x[i]<- rbinom(1, m, 1-exp(-y[i]))
     }
    xstart  <-  c(4.8,5.3)
    momEsta[b,]= nleqslv(xstart,model)$x
    
}
momEstb<-c(mean(momEsta[,1]),mean(momEsta[,2]))
momSe<-c(sd(momEsta[,1])/length(momEsta[,1]),sd(momEsta[,2])/length(momEsta[,2]))