#sample
library(LindleyR)
set.seed(1)
y<- rslindley(100, 5, 5, mixture = TRUE)
x<- rbinom(100, 12, 1-exp(-y))

#b=(m^2-((2*m^2-m)*(x[2]^2)*(x[2]+1+x[1]))/((x[2]+x[1])*(x[2]+1)^2)+((m^2-m)*(x[2]^2)*(x[2]+2+x[1]))/((x[2]+x[1])*(x[2]+2)^2)-10.95)

library(rootSolve)
m=12
model<- function(x){
   a=m*(x[2]^2+x[2]+2*x[1]*x[2]+x[1])/((x[2]+x[1])*(x[2]+1)^2)-2.43
   
   
   
   c= m^2-((2*m^2-m)*(x[2]^2)*(x[2]+1+x[1]))/((x[2]+x[1])*(x[2]+1)^2)+((m^2-m)*(x[2]^2)*(x[2]+2+x[1]))/((x[2]+x[1])*(x[2]+2)^2)-(m*(x[2]^2+x[2]+2*x[1]*x[2]+x[1])/((x[2]+x[1])*(x[2]+1)^2))^2 - 5.0451
   c(a=a,c=c)
}

ss<-multiroot(f=model,start=c(5.1,5.2))