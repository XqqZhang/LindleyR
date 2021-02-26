library(rootSolve)
m=12
model<- function(x){
   a=m*(x[2]^2+x[2]+2*x[1]*x[2]+x[1])/((x[2]+x[1])*(x[2]+1)^2)-2.43
   
   b=(m^2-((2*m^2-m)*(x[2]^2)*(x[2]+1+x[1]))/((x[2]+x[1])*(x[2]+1)^2)
       +((m^2-m)*(x[2]^2)*(x[2]+2+x[1]))/((x[2]+x[1])*(x[2]+2)^2)-10.95)
   
   c(a=a,b=b)
}

ss<-multiroot(f=model,start=c(5.1,5.2))



