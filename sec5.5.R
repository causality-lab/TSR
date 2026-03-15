#Experiment based on the IHDP (Infant Health and Development Program)

#use the same covariates like Hill (2011)
#package: bartcs
#variables (from the documentation)
#-treatment
#-y_factual: Observed outcome.
#-y_cfactual: Potential outcome given the opposite treatment.
#-mu0: Control conditional means.
#-mu1: Treated conditional means.
#-X1 ~ X6: Confounders with continuous values.
#-X7 ~ X25: Confounders with binary values.

#randomized experiment, but was transformed to an observational one

# We use the original dataset without deleted observation, which can be found here:
# https://www.tandfonline.com/doi/suppl/10.1198/jcgs.2010.08162?scroll=top

load(file="sim.data")
length(colnames(imp1))
#[1] 29

#identify the corresponding variables for treatment X and proxies Z={V1,...,V25}

#scale continuous covariates
Z<-cbind(scale(imp1[,2]),scale(imp1[,3]),scale(imp1[,4]),scale(imp1[,5]), scale(imp1[,6]),
         scale(imp1[,7]),      imp1[,8] ,      imp1[,9] ,      imp1[,10],       imp1[,11], 
               imp1[,12],      imp1[,13],      imp1[,14],      imp1[,15],       imp1[,16],  
               imp1[,17],      imp1[,18],      imp1[,19],      imp1[,20],       imp1[,21],
               imp1[,22],      imp1[,23],      imp1[,24],      imp1[,25], scale(imp1[,26]))

(n<-nrow(Z))
#[1] 985

set.seed(1)

tsr <- NULL
rr  <- NULL
na  <- NULL
t   <- NULL

for(i in 1:1000)
{
  #draw coefficients randomized from a discrete set
  values <- 0:4
  probabilities <- c(0.5, 0.2, 0.15, 0.1, 0.05) 
  discrete <- sample(values, size = 27, replace = TRUE, prob = probabilities)
  
  #randomized treatment
  X<-runif(n,-1,1)

  #generate target Y
  Y<-cbind(rep(1,n),X,Z)%*%discrete+rnorm(n)
  
  #induce selection bias based on X and one variable in Z
  S<-ifelse((X+Z[,3]>0.7),1,0)
  X_S<-X[S==1]
  Y_S<-Y[S==1]
  Z_S<-Z[S==1,]
  
  #ground truth causal effect
  t[i]<-discrete[2]
  
  #RR
  (beta_first<-lm(Y_S~X_S+Z_S)$coefficients)
  Y_second<-cbind(rep(1,n),X,Z)%*%beta_first
  (beta_second<-lm(Y_second~X)$coefficients)  
  g<-function(x){beta_second[1]+beta_second[2]*x}
  rr[i]<-g(1)-g(0)
  
  #TSR
  f<-function(x){beta_first[1]+beta_first[2]*x+sum(beta_first[-c(1,2)]*apply(Z,2,mean))}
  tsr[i]<-f(1)-f(0)
  
  #naive
  (beta_naive<-lm(Y_S~X_S)$coefficients)
  h<-function(x){beta_naive[1]+beta_naive[2]*x}
  na[i]<-h(1)-h(0)
}

par(mfrow=c(1,2),cex.axis = 1.6,cex.lab=1.6,oma=c(0,0.5,0,0),mar=c(7,7,3,3))
boxplot(cbind(tsr-t,rr-t,na-t),names=c("TSR","RR","naive"),col=c("lightgreen","lightblue",adjustcolor("orange", alpha = 0.35)),ylab="Bias",xlab="", main="Selection Bias")
abline(h=0,col="gray")
apply(cbind(tsr-t,rr-t,na-t),2,mean)
#[1]  0.006402254 -0.007051418  0.653637049

#adding confounding

tsr<-NULL
rr<-NULL
na<-NULL
t<-NULL

for(i in 1:1000)
{
  #draw coefficients randomized from a discrete set
  values <- 0:4
  probabilities <- c(0.5, 0.2, 0.15, 0.1, 0.05) 

  #randomized treatment
  X<-runif(n,-1,1)
  discrete <- sample(values, size = 28, replace = TRUE, prob = probabilities)
  
  #add confounding
  C<-rnorm(n,2,2)
  X<-X+C
  
  #generate target Y
  Y<-cbind(rep(1,n),X,Z,C)%*%discrete+rnorm(n)

  #induce selection bias based on X and one variable in Z
  S<-ifelse((X+Z[,3]>0.7),1,0)
  X_S<-X[S==1]
  Y_S<-Y[S==1]
  Z_S<-Z[S==1,]
  C_S<-C[S==1]
  
  #ground truth causal effect
  t[i]<-discrete[2]
  
  #RR
  (beta_first<-lm(Y_S~X_S+Z_S+C_S)$coefficients)
  Y_second<-cbind(rep(1,n),X,Z,C)%*%beta_first
  (beta_second<-lm(Y_second~X)$coefficients)  
  g<-function(x){beta_second[1]+beta_second[2]*x}
  rr[i]<-g(1)-g(0)
  
  #TSR
  f<-function(x){beta_first[1]+beta_first[2]*x+sum(beta_first[-c(1,2)]*apply(cbind(Z,C),2,mean))}
  tsr[i]<-f(1)-f(0)
  
  #naive
  (beta_naive<-lm(Y_S~X_S)$coefficients)
  h<-function(x){beta_naive[1]+beta_naive[2]*x}
  na[i]<-h(1)-h(0)
}

boxplot(cbind(tsr-t,rr-t,na-t),names=c("TSR","RR","naive"),col=c("lightgreen","lightblue",adjustcolor("orange", alpha = 0.35)),ylab="Bias",xlab="",main="Selection Bias and Confounding")
abline(h=0,col="gray")
apply(cbind(tsr-t,rr-t,na-t),2,mean)
#[1] 0.004511279 0.918960668 1.037122255

