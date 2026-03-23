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

#use package for cross-validation
library(glmnet)

load(file="sim.data")
length(colnames(imp1))
#[1] 29

#identify the corresponding variables for treatment X and proxies Z={V1,...,V25}

#scale continuous covariates
Z<-cbind(scale(imp1[,2]),scale(imp1[,3]),scale(imp1[,4]),scale(imp1[,5]),   scale(imp1[,6]),
         scale(imp1[,7]),      imp1[,8] ,      imp1[,9] ,      imp1[,10],         imp1[,11], 
               imp1[,12],      imp1[,13],      imp1[,14],      imp1[,15]-1,       imp1[,16],  
               imp1[,17],      imp1[,18],      imp1[,19],      imp1[,20],         imp1[,21],
               imp1[,22],      imp1[,23],      imp1[,24],      imp1[,25],         imp1[,26])

(n<-nrow(Z))
#[1] 985

set.seed(1)

tsr   <- NULL
tsr_r <- NULL
rr    <- NULL
rr_r  <- NULL
na    <- NULL
t     <- NULL

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
  rr[i]<-beta_second[2]
  
  #RR ridge 
  lambda_seq <- 10^seq(2, -2, by = -.1)
  ridge_cv <- cv.glmnet(cbind(X_S,Z_S), Y_S, alpha = 0, lambda = lambda_seq)
  best_lambda <- ridge_cv$lambda.min
  best_ridge <- glmnet(cbind(X_S,Z_S), Y_S, alpha = 0, lambda = best_lambda)
  beta_first_ridge<-as.vector(coef(best_ridge))
  Y_second_ridge<-cbind(rep(1,n),X,Z)%*%beta_first_ridge
  (beta_second_ridge<-lm(Y_second_ridge~X)$coefficients)  
  rr_r[i]<-beta_second_ridge[2]
  
  #TSR
  tsr[i]<-beta_first[2]
  
  #TSR ridge
  tsr_r[i]<-beta_first_ridge[2]
  
  #naive
  (beta_naive<-lm(Y_S~X_S)$coefficients)
  na[i]<-beta_naive[2]
}

par(mfrow=c(1,1),cex.axis = 2,cex.lab=2,cex.main=2.3,oma=c(0,0.5,0,0),mar=c(7,7,7,3))
boxplot(cbind(tsr-t,tsr_r-t,rr-t,rr_r-t,na-t),names=c("","","","",""),col=c("lightgreen","lightgreen","lightblue","lightblue",adjustcolor("orange", alpha = 0.35)),ylab="Bias",xlab="",main="Selection Bias")
axis(1,at = 1:5,labels = c("TSR\n","TSR \n(ridge)","RR\n", "RR \n(ridge)","naive\n"),tick = FALSE,line = 2.2)   
abline(h=0,col="gray")
apply(cbind(tsr-t,tsr_r-t,rr-t,rr_r-t,na-t),2,mean)
#[1] -0.0006897514 -0.0120191709 -0.0023520992 -0.0135540603  0.6375236977

#adding confounding

tsr   <- NULL
tsr_r <- NULL
rr    <- NULL
rr_r  <- NULL
na    <- NULL
t     <- NULL

for(i in 1:1000)
{
  #draw coefficients randomized from a discrete set
  values <- 0:4
  probabilities <- c(0.5, 0.2, 0.15, 0.1, 0.05) 

  #randomized treatment
  X<-runif(n,-1,1)
  discrete <- sample(values, size = 27, replace = TRUE, prob = probabilities)
  
  #add confounding
  C<-rnorm(n,2,2)
  X<-X+C
  
  #generate target Y
  Y<-cbind(rep(1,n),X,Z)%*%discrete+C+rnorm(n)

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
  rr[i]<-beta_second[2]
  
  #RR ridge 
  lambda_seq <- 10^seq(2, -2, by = -.1)
  ridge_cv <- cv.glmnet(cbind(X_S,Z_S,C_S), Y_S, alpha = 0, lambda = lambda_seq)
  best_lambda <- ridge_cv$lambda.min
  best_ridge <- glmnet(cbind(X_S,Z_S,C_S), Y_S, alpha = 0, lambda = best_lambda)
  beta_first_ridge<-as.vector(coef(best_ridge))
  Y_second_ridge<-cbind(rep(1,n),X,Z,C)%*%beta_first_ridge
  (beta_second_ridge<-lm(Y_second_ridge~X)$coefficients)  
  rr_r[i]<-beta_second_ridge[2]  
    
  #TSR
  tsr[i]<-beta_first[2]
  
  #TSR_ridge
  tsr_r[i]<-beta_first_ridge[2]
  
  #naive
  (beta_naive<-lm(Y_S~X_S)$coefficients)
  na[i]<-beta_naive[2]
}

par(mfrow=c(1,1),cex.axis = 2,cex.lab=2,cex.main=2.3,oma=c(0,0.5,0,0),mar=c(7,7,7,3))
boxplot(cbind(tsr-t,tsr_r-t,rr-t,rr_r-t,na-t),names=c("","","","",""),col=c("lightgreen","lightgreen","lightblue","lightblue",adjustcolor("orange", alpha = 0.35)),ylab="Bias",xlab="",main="Selection Bias\n and Confounding")
axis(1,at = 1:5,labels = c("TSR\n","TSR \n(ridge)","RR\n", "RR \n(ridge)","naive\n"),tick = FALSE,line = 2.2)   
abline(h=0,col="gray")
apply(cbind(tsr-t,tsr_r-t,rr-t,rr_r-t,na-t),2,mean)
#[1] -0.0005705826  0.0085811276  0.9204699508  0.9191743326  1.0486789300