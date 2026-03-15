#Section 5.2

#examples 5 & 6

#load packages
library(glmnet)
library(ggplot2)

#n = sample size
#SIM = number of simulation runs
#bsp = example number
#      note that bsp i corresponds to example i-2 (i.e. bsp3 corresponds to example 1)
#empty = indicates if ((S cap D)=emptyset)
#random = indicates uniformly distributed coefficients

sim56<-function(n,SIM,bsp,empty,random)
{
  #create empty matrices to save results
  
  vars_test <- c("w_test_naive", "w_test_est1", "w_test_est1_ridge", "w_test_est3", "w_test_est3_ridge", 
                 "o_test_naive", "o_test_est1", "o_test_est1_ridge", "o_test_est3", "o_test_est3_ridge")
  lapply(vars_test, function(x) assign(x, NULL, envir = .GlobalEnv))
  
  vars_test_S <- c("w_test_naive_S", "w_test_est1_S", "w_test_est1_ridge_S", "w_test_est3_S", "w_test_est3_ridge_S", 
                   "o_test_naive_S", "o_test_est1_S", "o_test_est1_ridge_S", "o_test_est3_S", "o_test_est3_ridge_S")
  lapply(vars_test_S, function(x) assign(x, NULL, envir = .GlobalEnv))
  
  bias_vars_test <- c("bias_w_test_naive", "bias_w_test_est1", "bias_w_test_est1_ridge", "bias_w_test_est3", "bias_w_test_est3_ridge",
                      "bias_o_test_naive", "bias_o_test_est1", "bias_o_test_est1_ridge", "bias_o_test_est3", "bias_o_test_est3_ridge")
  lapply(bias_vars_test, function(x) assign(x, NULL, envir = .GlobalEnv))
  
  bias_vars_test_S <- c("bias_w_test_naive_S", "bias_w_test_est1_S", "bias_w_test_est1_ridge_S", "bias_w_test_est3_S", "bias_w_test_est3_ridge_S",
                        "bias_o_test_naive_S", "bias_o_test_est1_S", "bias_o_test_est1_ridge_S", "bias_o_test_est3_S", "bias_o_test_est3_ridge_S")
  lapply(bias_vars_test_S, function(x) assign(x, NULL, envir = .GlobalEnv))
  
  if(bsp==5 | bsp==6 | bsp==7 | bsp==8){
    if(empty==F){
      X_complete <- matrix(ncol=SIM,nrow=n)
      Z_complete <- matrix(ncol=SIM,nrow=n)
      W_complete <- matrix(ncol=SIM,nrow=n)
      S_complete <- matrix(ncol=SIM,nrow=n)
      Y_complete <- matrix(ncol=SIM,nrow=n)
    }else{
      X_Sdata_complete <- matrix(ncol=SIM,nrow=n)
      Z_Sdata_complete <- matrix(ncol=SIM,nrow=n)
      W_Sdata_complete <- matrix(ncol=SIM,nrow=n)
      S_Sdata_complete <- matrix(ncol=SIM,nrow=n)
      Y_Sdata_complete <- matrix(ncol=SIM,nrow=n)
      
      X_Ddata_complete <- matrix(ncol=SIM,nrow=n)
      Z_Ddata_complete <- matrix(ncol=SIM,nrow=n)
      W_Ddata_complete <- matrix(ncol=SIM,nrow=n)}}
  
  if(bsp==5 | bsp==6 | bsp==7 | bsp==8){
    betahat_first_complete <- matrix(ncol=SIM,nrow=10)
    betahat_first_ridge_complete <- matrix(ncol=SIM,nrow=10)
    beta_end1_complete <- matrix(ncol=SIM, nrow=4)
    beta_end1_ridge_complete <- matrix(ncol=SIM, nrow=4)
    naive_beta_complete <- matrix(ncol=SIM, nrow=4)
    y0_complete<-NULL
    yx1_complete<-NULL
    yx2_complete<-NULL
    yx3_complete<-NULL
    yz1_complete<-NULL
    yz2_complete<-NULL
    yz3_complete<-NULL
    yw1_complete<-NULL
    yw2_complete<-NULL
    yw3_complete<-NULL
  }
  
  if(bsp==7 | bsp==8){
    beta_ZlmX_complete <- matrix(ncol=SIM, nrow=4)
    beta_Z2lmX_complete <- matrix(ncol=SIM, nrow=4)
    beta_Z3lmX_complete <- matrix(ncol=SIM, nrow=4)
    y0_complete<-NULL
    yx1_complete<-NULL
    yx2_complete<-NULL
    yx3_complete<-NULL
    yz1_complete<-NULL
    yz2_complete<-NULL
    yz3_complete<-NULL
    yw1_complete<-NULL
    yw2_complete<-NULL
    yw3_complete<-NULL
  }
  
  for(sim in 1:SIM)
  {
    print(sim)
    
    #reproducibility
    set.seed(sim)
    
    #data generating processes
      #example 5
      if(bsp==7){
        if(empty==F){
          W<-rnorm(mean=-1,sd=1,n=n)
          X<-rnorm(mean=0,sd=1,n=n)+W
          Z<--2*X+rnorm(mean=0,sd=1,n=n)
          S<-rbinom(p=1/((1+exp(X))*(1+exp(Z))),size=1,n)
          if(random==F){
            y0<--0.2
            yx1<-0
            yx2<-0
            yx3<-0.001
            yz1<-0
            yz2<-0
            yz3<--0.002
            yw1<-0
            yw2<-1
            yw3<-1
          }else{
            y0<-runif(1)-0.7
            yx1<-runif(1)-0.5
            yx2<-runif(1)-0.5
            yx3<-runif(1)-0.5+0.001
            yz1<-runif(1)-0.5
            yz2<-runif(1)-0.5
            yz3<-runif(1)-0.5-0.002
            yw1<-runif(1)-0.5
            yw2<-runif(1)+0.5
            yw3<-runif(1)+0.5
          }
          y0_complete[sim]<-y0
          yx1_complete[sim]<-yx1
          yx2_complete[sim]<-yx2
          yx3_complete[sim]<-yx3
          yz1_complete[sim]<-yz1
          yz2_complete[sim]<-yz2
          yz3_complete[sim]<-yz3
          yw1_complete[sim]<-yw1
          yw2_complete[sim]<-yw2
          yw3_complete[sim]<-yw3
          Y<-y0+yx1*X+yx2*X^2+yx3*X^3+yz1*Z+yz2*Z^2+yz3*Z^3+yw1*W+yw2*W^2+yw3*W^3+rnorm(n)  
          
          Z_S<-Z[S==1]
          W_S<-W[S==1]
          X_S<-X[S==1]
          Y_S<-Y[S==1]
        }else{
          W_Sdata<-rnorm(mean=-1,sd=1,n=n)
          X_Sdata<-rnorm(mean=0,sd=1,n=n)+W_Sdata
          Z_Sdata<--2*X_Sdata+rnorm(mean=0,sd=1,n=n)
          S_Sdata<-rbinom(p=1/((1+exp(X_Sdata))*(1+exp(Z_Sdata))),size=1,n)
          if(random==F){
            y0<--0.2
            yx1<-0
            yx2<-0
            yx3<-0.001
            yz1<-0
            yz2<-0
            yz3<--0.002
            yw1<-0
            yw2<-1
            yw3<-1
          }else{
            y0<-runif(1)-0.7
            yx1<-runif(1)-0.5
            yx2<-runif(1)-0.5
            yx3<-runif(1)-0.5+0.001
            yz1<-runif(1)-0.5
            yz2<-runif(1)-0.5
            yz3<-runif(1)-0.5-0.002
            yw1<-runif(1)-0.5
            yw2<-runif(1)+0.5
            yw3<-runif(1)+0.5
          }
          y0_complete[sim]<-y0
          yx1_complete[sim]<-yx1
          yx2_complete[sim]<-yx2
          yx3_complete[sim]<-yx3
          yz1_complete[sim]<-yz1
          yz2_complete[sim]<-yz2
          yz3_complete[sim]<-yz3
          yw1_complete[sim]<-yw1
          yw2_complete[sim]<-yw2
          yw3_complete[sim]<-yw3
          Y_Sdata<-y0+yx1*X_Sdata+yx2*X_Sdata^2+yx3*X_Sdata^3+yz1*Z_Sdata+yz2*Z_Sdata^2+yz3*Z_Sdata^3+yw1*W_Sdata+yw2*W_Sdata^2+yw3*W_Sdata^3+rnorm(n)
          
          Z_S<-Z_Sdata[S_Sdata==1]
          W_S<-W_Sdata[S_Sdata==1]
          X_S<-X_Sdata[S_Sdata==1]
          Y_S<-Y_Sdata[S_Sdata==1]
          
          W_Ddata<-rnorm(mean=-1,sd=1,n=n)
          X_Ddata<-rnorm(mean=0,sd=1,n=n)+W_Ddata
          Z_Ddata<--2*X_Ddata+rnorm(mean=0,sd=1,n=n)
          
          X<-X_Ddata
          Z<-Z_Ddata
          W<-W_Ddata}}
      #example 6
      if(bsp==8){
        if(empty==F){
          W<-rnorm(mean=2,sd=1,n=n)
          X<-rnorm(mean=0,sd=1,n=n)+W
          Z<-X+rnorm(mean=0,sd=1,n=n)
          S<-ifelse((Z*X)<4&(Z^2*X^2+Z)>0,1,0)
          if(random==F){
            y0<-0
            yx1<-0.5
            yx2<-0
            yx3<-0
            yz1<-0.05
            yz2<-0
            yz3<-0
            yw1<-1
            yw2<-0
            yw3<-0
          }else{
            y0<-runif(1)-0.5
            yx1<-runif(1)
            yx2<-runif(1)-0.5
            yx3<-runif(1)-0.5
            yz1<-runif(1)-0.5+0.05
            yz2<-runif(1)-0.5
            yz3<-runif(1)-0.5
            yw1<-runif(1)+0.5
            yw2<-runif(1)-0.5
            yw3<-runif(1)-0.5
          }
          y0_complete[sim]<-y0
          yx1_complete[sim]<-yx1
          yx2_complete[sim]<-yx2
          yx3_complete[sim]<-yx3
          yz1_complete[sim]<-yz1
          yz2_complete[sim]<-yz2
          yz3_complete[sim]<-yz3
          yw1_complete[sim]<-yw1
          yw2_complete[sim]<-yw2
          yw3_complete[sim]<-yw3
          Y<-y0+yx1*X+yx2*X^2+yx3*X^3+yz1*Z+yz2*Z^2+yz3*Z^3+yw1*W+yw2*W^2+yw3*W^3+rnorm(n)  
          
          Z_S<-Z[S==1]
          W_S<-W[S==1]
          X_S<-X[S==1]
          Y_S<-Y[S==1]
        }else{
          W_Sdata<-rnorm(mean=2,sd=1,n=n)
          X_Sdata<-rnorm(mean=0,sd=1,n=n)+W_Sdata
          Z_Sdata<-X_Sdata+rnorm(mean=0,sd=1,n=n)
          S_Sdata<-ifelse((Z_Sdata*X_Sdata)<4&(Z_Sdata^2*X_Sdata^2+Z_Sdata)>0,1,0)
          if(random==F){
            y0<-0
            yx1<-0.5
            yx2<-0
            yx3<-0
            yz1<-0.05
            yz2<-0
            yz3<-0
            yw1<-1
            yw2<-0
            yw3<-0
          }else{
            y0<-runif(1)-0.5
            yx1<-runif(1)
            yx2<-runif(1)-0.5
            yx3<-runif(1)-0.5
            yz1<-runif(1)-0.5+0.05
            yz2<-runif(1)-0.5
            yz3<-runif(1)-0.5
            yw1<-runif(1)+0.5
            yw2<-runif(1)-0.5
            yw3<-runif(1)-0.5
          }
          y0_complete[sim]<-y0
          yx1_complete[sim]<-yx1
          yx2_complete[sim]<-yx2
          yx3_complete[sim]<-yx3
          yz1_complete[sim]<-yz1
          yz2_complete[sim]<-yz2
          yz3_complete[sim]<-yz3
          yw1_complete[sim]<-yw1
          yw2_complete[sim]<-yw2
          yw3_complete[sim]<-yw3
          Y_Sdata<-y0+yx1*X_Sdata+yx2*X_Sdata^2+yx3*X_Sdata^3+yz1*Z_Sdata+yz2*Z_Sdata^2+yz3*Z_Sdata^3+yw1*W_Sdata+yw2*W_Sdata^2+yw3*W_Sdata^3+rnorm(n)  
          
          Z_S<-Z_Sdata[S_Sdata==1]
          W_S<-W_Sdata[S_Sdata==1]
          X_S<-X_Sdata[S_Sdata==1]
          Y_S<-Y_Sdata[S_Sdata==1]
          
          W_Ddata<-rnorm(mean=2,sd=1,n=n)
          X_Ddata<-rnorm(mean=0,sd=1,n=n)+W_Ddata
          Z_Ddata<-X_Ddata+rnorm(mean=0,sd=1,n=n)
          
          X<-X_Ddata
          Z<-Z_Ddata
          W<-W_Ddata}}
      
      #save realisations from all simulation runs
      if(bsp==5 | bsp==6 | bsp==7 | bsp==8){
        if(empty==F){
          X_complete[,sim]<-X
          Z_complete[,sim]<-Z
          W_complete[,sim]<-W
          S_complete[,sim]<-S
          Y_complete[,sim]<-Y
        }else{
          X_Sdata_complete[,sim] <- X_Sdata
          Z_Sdata_complete[,sim] <- Z_Sdata
          W_Sdata_complete[,sim] <- W_Sdata
          S_Sdata_complete[,sim] <- S_Sdata
          Y_Sdata_complete[,sim] <- Y_Sdata
          
          X_Ddata_complete[,sim] <- X_Ddata
          Z_Ddata_complete[,sim] <- Z_Ddata
          W_Ddata_complete[,sim] <- W_Ddata}}
      
      #first step ridge 
      if(bsp==5 | bsp==6 | bsp==7 | bsp==8){
        lambda_seq <- 10^seq(2, -2, by = -.1)
        ridge_cv <- cv.glmnet(cbind(X_S,I(X_S^2),I(X_S^3),Z_S,I(Z_S^2),I(Z_S^3),W_S,I(W_S^2),I(W_S^3)), Y_S, alpha = 0, lambda = lambda_seq)
        best_lambda <- ridge_cv$lambda.min
        best_ridge <- glmnet(cbind(X_S,I(X_S^2),I(X_S^3),Z_S,I(Z_S^2),I(Z_S^3),W_S,I(W_S^2),I(W_S^3)), Y_S, alpha = 0, lambda = best_lambda)
        betahat_first_ridge<-as.vector(coef(best_ridge))}
      
      #first step OLS
      if(bsp==5 | bsp==6 | bsp==7 | bsp==8){
        lm_fit <- lm(Y_S~X_S+I(X_S^2)+I(X_S^3)+Z_S+I(Z_S^2)+I(Z_S^3)+W_S+I(W_S^2)+I(W_S^3))
        betahat_first<-coef(lm_fit)}
      
      #calculate target variable for second step of RR
      if(bsp==5 | bsp==6 | bsp==7 | bsp==8){
        Y_second_ridge<-(cbind(rep(1,n),X,X^2,X^3,Z,Z^2,Z^3,W,W^2,W^3))%*%as.vector(betahat_first_ridge)
        Y_second<-(cbind(rep(1,n),X,X^2,X^3,Z,Z^2,Z^3,W,W^2,W^3))%*%as.vector(betahat_first)}
      
      #final RR estimates
      beta_end1_ridge<-lm(Y_second_ridge~1+X+I(X^2)+I(X^3))$coef
      beta_end1<-lm(Y_second~1+X+I(X^2)+I(X^3))$coef
      
      #naive estimation only based on S=1
      naive_beta<-lm(Y_S~1+X_S+I(X_S^2)+I(X_S^3))$coef
      
      #second step of TSR
      if(bsp==7 | bsp==8){
        beta_ZlmX<-lm(Z~1+X+I(X^2)+I(X^3))$coefficients
        beta_Z2lmX<-lm(Z^2~1+X+I(X^2)+I(X^3))$coefficients
        beta_Z3lmX<-lm(Z^3~1+X+I(X^2)+I(X^3))$coefficients}
      
      
      #E[Y|do(X)]
      #example 5
      if(bsp==7){
        original<-function(x){
          y0+yx1*x+yx2*x^2+yx3*x^3+yz1*(-2*x)+yz2*(4*x^2+1)+yz3*(-8*x^3-6*x)+yw1*(-1)+yw2*(2)+yw3*(-4)}}
      
      #example 6
      if(bsp==8){
        original<-function(x){
          y0+yx1*x+yx2*x^2+yx3*x^3+yz1*(1.5*x-1)+yz2*(x^2+1)+yz3*(x^3+3*x)+yw1*(2)+yw2*(5)+yw3*(14)}}
      
      
      #E[Y|X]
      #example 5
      if(bsp==7){
        wrong<-function(x){
          y0+yx1*x+yx2*x^2+yx3*x^3+yz1*(-2*x)+yz2*(4*x^2+1)+yz3*(-8*x^3-6*x)+yw1*(-0.5+0.5*x)+yw2*(1-x)+yw3*(-1+3*x)}}
      
      #example 6
      if(bsp==8){
        wrong<-function(x){
          y0+yx1*x+yx2*x^2+yx3*x^3+yz1*(1.5*x-1)+yz2*(x^2+1)+yz3*(x^3+3*x)+yw1*(1+0.5*x)+yw2*(1+2*x)+yw3*(7.5*x-1)}}
      
      #RR ridge
      est1_ridge<-function(x)
      {beta_end1_ridge[1]+beta_end1_ridge[2]*x+beta_end1_ridge[3]*x^2+beta_end1_ridge[4]*x^3}
      
      #RR OLS
      est1<-function(x)
      {beta_end1[1]+beta_end1[2]*x+beta_end1[3]*x^2+beta_end1[4]*x^3}
      
      #naive 
      naive<-function(x)
      {naive_beta[1]+naive_beta[2]*x+naive_beta[3]*x^2+naive_beta[4]*x^3}
      
      #TSR ridge
      if(bsp==7 | bsp==8){
        est3_ridge<-function(x){
          return(betahat_first_ridge[1]+betahat_first_ridge[2]*x+betahat_first_ridge[3]*x^2+betahat_first_ridge[4]*x^3+
                   betahat_first_ridge[5]*(beta_ZlmX[1]+beta_ZlmX[2]*x+beta_ZlmX[3]*x^2+beta_ZlmX[4]*x^3)+
                   betahat_first_ridge[6]*(beta_Z2lmX[1]+beta_Z2lmX[2]*x+beta_Z2lmX[3]*x^2+beta_Z2lmX[4]*x^3)+
                   betahat_first_ridge[7]*(beta_Z3lmX[1]+beta_Z3lmX[2]*x+beta_Z3lmX[3]*x^2+beta_Z3lmX[4]*x^3)+
                   betahat_first_ridge[8]*mean(W)+
                   betahat_first_ridge[9]*mean(W^2)+
                   betahat_first_ridge[10]*mean(W^3))}}
      
      #TSR OLS
      if(bsp==7 | bsp==8){
        est3<-function(x){
          return(betahat_first[1]+betahat_first[2]*x+betahat_first[3]*x^2+betahat_first[4]*x^3+
                   betahat_first[5]*(beta_ZlmX[1]+beta_ZlmX[2]*x+beta_ZlmX[3]*x^2+beta_ZlmX[4]*x^3)+
                   betahat_first[6]*(beta_Z2lmX[1]+beta_Z2lmX[2]*x+beta_Z2lmX[3]*x^2+beta_Z2lmX[4]*x^3)+
                   betahat_first[7]*(beta_Z3lmX[1]+beta_Z3lmX[2]*x+beta_Z3lmX[3]*x^2+beta_Z3lmX[4]*x^3)+
                   betahat_first[8]*mean(W)+
                   betahat_first[9]*mean(W^2)+
                   betahat_first[10]*mean(W^3))}}
      
      
      #save coefficient vectors from all simulation runs
      betahat_first_complete[,sim] <- betahat_first
      betahat_first_ridge_complete[,sim] <- betahat_first_ridge
      beta_end1_complete[,sim] <- beta_end1
      beta_end1_ridge_complete[,sim] <- beta_end1_ridge
      naive_beta_complete[,sim] <- naive_beta
      
      if(bsp==7 | bsp==8){
        beta_ZlmX_complete[,sim]<-beta_ZlmX
        beta_Z2lmX_complete[,sim]<-beta_Z2lmX
        beta_Z3lmX_complete[,sim]<-beta_Z3lmX}
      
      
      if(random==T){
      #test data    
      if(bsp==7){
        W_test<-rnorm(mean=-1,sd=1,n=n)
        X_test<-as.matrix(rnorm(mean=0,sd=1,n=n)+W_test)
        Z_test<--2*X_test+rnorm(mean=0,sd=1,n=n)
        S_test<-rbinom(p=1/((1+exp(X_test))*(1+exp(Z_test))),size=1,n)
        Y_test<-y0+yx1*X_test+yx2*X_test^2+yx3*X_test^3+yz1*Z_test+yz2*Z_test^2+yz3*Z_test^3+yw1*W_test+yw2*W_test^2+yw3*W_test^3    
        X_test_S<-as.matrix(X_test[S_test==1,])
        Y_test_S<-Y_test[S_test==1]}
      if(bsp==8){
        W_test<-rnorm(mean=2,sd=1,n=n)
        X_test<-as.matrix(rnorm(mean=0,sd=1,n=n)+W_test)
        Z_test<-X_test+rnorm(mean=0,sd=1,n=n)
        S_test<-ifelse((Z_test*X_test)<4&(Z_test^2*X_test^2+Z_test)>0,1,0)
        Y_test<-y0+yx1*X_test+yx2*X_test^2+yx3*X_test^3+yz1*Z_test+yz2*Z_test^2+yz3*Z_test^3+yw1*W_test+yw2*W_test^2+yw3*W_test^3  
        X_test_S<-as.matrix(X_test[S_test==1,])
        Y_test_S<-Y_test[S_test==1]}
      
      #MSE (w(rong): E[Y|X], o(riginal): E[Y|do(X)]) 
      
      #in D
      o_test_naive[sim]<-mean((apply(X_test,1,original)-apply(X_test,1,naive))^2)
      o_test_est1[sim]<-mean((apply(X_test,1,original)-apply(X_test,1,est1))^2)
      o_test_est1_ridge[sim]<-mean((apply(X_test,1,original)-apply(X_test,1,est1_ridge))^2)
      o_test_est3[sim]<-mean((apply(X_test,1,original)-apply(X_test,1,est3))^2)
      o_test_est3_ridge[sim]<-mean((apply(X_test,1,original)-apply(X_test,1,est3_ridge))^2)
      
      w_test_naive[sim]<-mean((apply(X_test,1,wrong)-apply(X_test,1,naive))^2)
      w_test_est1[sim]<-mean((apply(X_test,1,wrong)-apply(X_test,1,est1))^2)
      w_test_est1_ridge[sim]<-mean((apply(X_test,1,wrong)-apply(X_test,1,est1_ridge))^2)
      w_test_est3[sim]<-mean((apply(X_test,1,wrong)-apply(X_test,1,est3))^2)
      w_test_est3_ridge[sim]<-mean((apply(X_test,1,wrong)-apply(X_test,1,est3_ridge))^2)
      
      #in S
      o_test_naive_S[sim]<-mean((apply(X_test_S,1,original)-apply(X_test_S,1,naive))^2)
      o_test_est1_S[sim]<-mean((apply(X_test_S,1,original)-apply(X_test_S,1,est1))^2)
      o_test_est1_ridge_S[sim]<-mean((apply(X_test_S,1,original)-apply(X_test_S,1,est1_ridge))^2)
      o_test_est3_S[sim]<-mean((apply(X_test_S,1,original)-apply(X_test_S,1,est3))^2)
      o_test_est3_ridge_S[sim]<-mean((apply(X_test_S,1,original)-apply(X_test_S,1,est3_ridge))^2)
      
      w_test_naive_S[sim]<-mean((apply(X_test_S,1,wrong)-apply(X_test_S,1,naive))^2)
      w_test_est1_S[sim]<-mean((apply(X_test_S,1,wrong)-apply(X_test_S,1,est1))^2)
      w_test_est1_ridge_S[sim]<-mean((apply(X_test_S,1,wrong)-apply(X_test_S,1,est1_ridge))^2)
      w_test_est3_S[sim]<-mean((apply(X_test_S,1,wrong)-apply(X_test_S,1,est3))^2)
      w_test_est3_ridge_S[sim]<-mean((apply(X_test_S,1,wrong)-apply(X_test_S,1,est3_ridge))^2)
      
      #Bias (w: E[Y|X], o: E[Y|do(X)]) 
      
      #in D
      #bias_o_test_naive[sim]<-mean((apply(X_test,1,original)-apply(X_test,1,naive)))
      #bias_o_test_est1[sim]<-mean((apply(X_test,1,original)-apply(X_test,1,est1)))
      #bias_o_test_est1_ridge[sim]<-mean((apply(X_test,1,original)-apply(X_test,1,est1_ridge)))
      #bias_o_test_est3[sim]<-mean((apply(X_test,1,original)-apply(X_test,1,est3)))
      #bias_o_test_est3_ridge[sim]<-mean((apply(X_test,1,original)-apply(X_test,1,est3_ridge)))
      
      #bias_w_test_naive[sim]<-mean((apply(X_test,1,wrong)-apply(X_test,1,naive)))
      #bias_w_test_est1[sim]<-mean((apply(X_test,1,wrong)-apply(X_test,1,est1)))
      #bias_w_test_est1_ridge[sim]<-mean((apply(X_test,1,wrong)-apply(X_test,1,est1_ridge)))
      #bias_w_test_est3[sim]<-mean((apply(X_test,1,wrong)-apply(X_test,1,est3)))
      #bias_w_test_est3_ridge[sim]<-mean((apply(X_test,1,wrong)-apply(X_test,1,est3_ridge)))
      
      #in S
      #bias_o_test_naive_S[sim]<-mean((apply(X_test_S,1,original)-apply(X_test_S,1,naive)))
      #bias_o_test_est1_S[sim]<-mean((apply(X_test_S,1,original)-apply(X_test_S,1,est1)))
      #bias_o_test_est1_ridge_S[sim]<-mean((apply(X_test_S,1,original)-apply(X_test_S,1,est1_ridge)))
      #bias_o_test_est3_S[sim]<-mean((apply(X_test_S,1,original)-apply(X_test_S,1,est3)))
      #bias_o_test_est3_ridge_S[sim]<-mean((apply(X_test_S,1,original)-apply(X_test_S,1,est3_ridge)))
      
      #bias_w_test_naive_S[sim]<-mean((apply(X_test_S,1,wrong)-apply(X_test_S,1,naive)))
      #bias_w_test_est1_S[sim]<-mean((apply(X_test_S,1,wrong)-apply(X_test_S,1,est1)))
      #bias_w_test_est1_ridge_S[sim]<-mean((apply(X_test_S,1,wrong)-apply(X_test_S,1,est1_ridge)))
      #bias_w_test_est3_S[sim]<-mean((apply(X_test_S,1,wrong)-apply(X_test_S,1,est3)))
      #bias_w_test_est3_ridge_S[sim]<-mean((apply(X_test_S,1,wrong)-apply(X_test_S,1,est3_ridge)))
    }
  }
    #return the relevant results  
    if(bsp==7 | bsp==8){      
      if(empty==F){ 
        return(list(SIM = SIM, n = n , bsp = bsp, empty=empty,
                    X_complete = X_complete, Z_complete = Z_complete, W_complete = W_complete, S_complete = S_complete, Y_complete = Y_complete,
                    betahat_first_complete = betahat_first_complete, betahat_first_ridge_complete = betahat_first_ridge_complete,
                    beta_end1_complete = beta_end1_complete, beta_end1_ridge_complete = beta_end1_ridge_complete,
                    naive_beta_complete = naive_beta_complete,
                    o_test_naive=o_test_naive, o_test_est1=o_test_est1, o_test_est1_ridge=o_test_est1_ridge, o_test_est3=o_test_est3, o_test_est3_ridge=o_test_est3_ridge,
                    o_test_naive_S=o_test_naive_S, o_test_est1_S=o_test_est1_S, o_test_est1_ridge_S=o_test_est1_ridge_S, o_test_est3_S=o_test_est3_S, o_test_est3_ridge_S=o_test_est3_ridge_S,
                    w_test_naive=w_test_naive, w_test_est1=w_test_est1, w_test_est1_ridge=w_test_est1_ridge, w_test_est3=w_test_est3, w_test_est3_ridge=w_test_est3_ridge,
                    w_test_naive_S=w_test_naive_S, w_test_est1_S=w_test_est1_S, w_test_est1_ridge_S=w_test_est1_ridge_S, w_test_est3_S=w_test_est3_S, w_test_est3_ridge_S=w_test_est3_ridge_S,
                    bias_o_test_naive=bias_o_test_naive, bias_o_test_est1=bias_o_test_est1, bias_o_test_est1_ridge=bias_o_test_est1_ridge, bias_o_test_est3=bias_o_test_est3, bias_o_test_est3_ridge=bias_o_test_est3_ridge,
                    bias_o_test_naive_S=bias_o_test_naive_S, bias_o_test_est1_S=bias_o_test_est1_S, bias_o_test_est1_ridge_S=bias_o_test_est1_ridge_S, bias_o_test_est3_S=bias_o_test_est3_S, bias_o_test_est3_ridge_S=bias_o_test_est3_ridge_S,
                    bias_w_test_naive=bias_w_test_naive, bias_w_test_est1=bias_w_test_est1, bias_w_test_est1_ridge=bias_w_test_est1_ridge, bias_w_test_est3=bias_w_test_est3, bias_w_test_est3_ridge=bias_w_test_est3_ridge,
                    bias_w_test_naive_S=bias_w_test_naive_S, bias_w_test_est1_S=bias_w_test_est1_S, bias_w_test_est1_ridge_S=bias_w_test_est1_ridge_S, bias_w_test_est3_S=bias_w_test_est3_S, bias_w_test_est3_ridge_S=bias_w_test_est3_ridge_S,
                    beta_ZlmX_complete=beta_ZlmX_complete, beta_Z2lmX_complete=beta_Z2lmX_complete,beta_Z3lmX_complete=beta_Z3lmX_complete,
                    y0_complete=y0_complete,yx1_complete=yx1_complete,yx2_complete=yx2_complete,yx3_complete=yx3_complete,yz1_complete=yz1_complete,yz2_complete=yz2_complete,yz3_complete=yz3_complete,yw1_complete=yw1_complete,yw2_complete=yw2_complete,yw3_complete=yw3_complete))
      }else{
        return(list(SIM = SIM, n = n , bsp = bsp, 
                    X_Sdata_complete = X_Sdata_complete, Z_Sdata_complete = Z_Sdata_complete,  W_Sdata_complete = W_Sdata_complete, S_Sdata_complete = S_Sdata_complete, Y_Sdata_complete = Y_Sdata_complete,
                    X_Ddata_complete = X_Ddata_complete, Z_Ddata_complete = Z_Ddata_complete,  W_Ddata_complete = W_Ddata_complete,
                    betahat_first_complete = betahat_first_complete, betahat_first_ridge_complete = betahat_first_ridge_complete,
                    beta_end1_complete = beta_end1_complete, beta_end1_ridge_complete = beta_end1_ridge_complete,
                    naive_beta_complete = naive_beta_complete,
                    o_test_naive=o_test_naive, o_test_est1=o_test_est1, o_test_est1_ridge=o_test_est1_ridge, o_test_est3=o_test_est3, o_test_est3_ridge=o_test_est3_ridge,
                    o_test_naive_S=o_test_naive_S, o_test_est1_S=o_test_est1_S, o_test_est1_ridge_S=o_test_est1_ridge_S, o_test_est3_S=o_test_est3_S, o_test_est3_ridge_S=o_test_est3_ridge_S,
                    w_test_naive=w_test_naive, w_test_est1=w_test_est1, w_test_est1_ridge=w_test_est1_ridge, w_test_est3=w_test_est3, w_test_est3_ridge=w_test_est3_ridge,
                    w_test_naive_S=w_test_naive_S, w_test_est1_S=w_test_est1_S, w_test_est1_ridge_S=w_test_est1_ridge_S, w_test_est3_S=w_test_est3_S, w_test_est3_ridge_S=w_test_est3_ridge_S,
                    bias_o_test_naive=bias_o_test_naive, bias_o_test_est1=bias_o_test_est1, bias_o_test_est1_ridge=bias_o_test_est1_ridge, bias_o_test_est3=bias_o_test_est3, bias_o_test_est3_ridge=bias_o_test_est3_ridge,
                    bias_o_test_naive_S=bias_o_test_naive_S, bias_o_test_est1_S=bias_o_test_est1_S, bias_o_test_est1_ridge_S=bias_o_test_est1_ridge_S, bias_o_test_est3_S=bias_o_test_est3_S, bias_o_test_est3_ridge_S=bias_o_test_est3_ridge_S,
                    bias_w_test_naive=bias_w_test_naive, bias_w_test_est1=bias_w_test_est1, bias_w_test_est1_ridge=bias_w_test_est1_ridge, bias_w_test_est3=bias_w_test_est3, bias_w_test_est3_ridge=bias_w_test_est3_ridge,
                    bias_w_test_naive_S=bias_w_test_naive_S, bias_w_test_est1_S=bias_w_test_est1_S, bias_w_test_est1_ridge_S=bias_w_test_est1_ridge_S, bias_w_test_est3_S=bias_w_test_est3_S, bias_w_test_est3_ridge_S=bias_w_test_est3_ridge_S,
                    beta_ZlmX_complete=beta_ZlmX_complete, beta_Z2lmX_complete=beta_Z2lmX_complete, beta_Z3lmX_complete=beta_Z3lmX_complete,
                    y0_complete=y0_complete,yx1_complete=yx1_complete,yx2_complete=yx2_complete,yx3_complete=yx3_complete,yz1_complete=yz1_complete,yz2_complete=yz2_complete,yz3_complete=yz3_complete,yw1_complete=yw1_complete,yw2_complete=yw2_complete,yw3_complete=yw3_complete))}}
  }
  
  
  
  
  #####################################################################################################
  #####################################################################################################
  #####################################################################################################
  

#subset

result500_bsp7 <- sim56(n=500,SIM=100,bsp=7,empty=F,random=F)
result1000_bsp7 <- sim56(n=1000,SIM=100,bsp=7,empty=F,random=F)
result5000_bsp7 <- sim56(n=5000,SIM=100,bsp=7,empty=F,random=F)

result500_bsp7_rand <- sim56(n=500,SIM=100,bsp=7,empty=F,random=T)
result1000_bsp7_rand <- sim56(n=1000,SIM=100,bsp=7,empty=F,random=T)
result5000_bsp7_rand <- sim56(n=5000,SIM=100,bsp=7,empty=F,random=T)

result500_bsp8 <- sim56(n=500,SIM=100,bsp=8,empty=F,random=F)
result1000_bsp8 <- sim56(n=1000,SIM=100,bsp=8,empty=F,random=F)
result5000_bsp8 <- sim56(n=5000,SIM=100,bsp=8,empty=F,random=F)

result500_bsp8_rand <- sim56(n=500,SIM=100,bsp=8,empty=F,random=T)
result1000_bsp8_rand <- sim56(n=1000,SIM=100,bsp=8,empty=F,random=T)
result5000_bsp8_rand <- sim56(n=5000,SIM=100,bsp=8,empty=F,random=T)

#empty

empty_result500_bsp7 <- sim56(n=500,SIM=100,bsp=7,empty=T,random=F)
empty_result1000_bsp7 <- sim56(n=1000,SIM=100,bsp=7,empty=T,random=F)
empty_result5000_bsp7 <- sim56(n=5000,SIM=100,bsp=7,empty=T,random=F)

empty_result500_bsp7_rand <- sim56(n=500,SIM=100,bsp=7,empty=T,random=T)
empty_result1000_bsp7_rand <- sim56(n=1000,SIM=100,bsp=7,empty=T,random=T)
empty_result5000_bsp7_rand <- sim56(n=5000,SIM=100,bsp=7,empty=T,random=T)

empty_result500_bsp8 <- sim56(n=500,SIM=100,bsp=8,empty=T,random=F)
empty_result1000_bsp8 <- sim56(n=1000,SIM=100,bsp=8,empty=T,random=F)
empty_result5000_bsp8 <- sim56(n=5000,SIM=100,bsp=8,empty=T,random=F)

empty_result500_bsp8_rand <- sim56(n=500,SIM=100,bsp=8,empty=T,random=T)
empty_result1000_bsp8_rand <- sim56(n=1000,SIM=100,bsp=8,empty=T,random=T)
empty_result5000_bsp8_rand <- sim56(n=5000,SIM=100,bsp=8,empty=T,random=T)


#  save(result500_bsp7,result1000_bsp7,result5000_bsp7,
#       result500_bsp7_rand,result1000_bsp7_rand,result5000_bsp7_rand,
#       result500_bsp8,result1000_bsp8,result5000_bsp8,
#       result500_bsp8_rand,result1000_bsp8_rand,result5000_bsp8_rand,
#       empty_result500_bsp7,empty_result1000_bsp7,empty_result5000_bsp7,
#       empty_result500_bsp7_rand,empty_result1000_bsp7_rand,empty_result5000_bsp7_rand,
#       empty_result500_bsp8,empty_result1000_bsp8,empty_result5000_bsp8,
#       empty_result500_bsp8_rand,empty_result1000_bsp8_rand,empty_result5000_bsp8_rand,
#       file="resultssim78.RData")  

load("resultssim78.RData")


library(ggplot2)



#define middle_green
green_rgb <- col2rgb("green")
darkgreen_rgb <- col2rgb("darkgreen")
middle_green_rgb <- (green_rgb + darkgreen_rgb) / 2
middle_green <- rgb(middle_green_rgb[1], middle_green_rgb[2], middle_green_rgb[3], maxColorValue = 255)



result<-result5000_bsp8
empty<-FALSE
set.seed(1)

SIM<-result$SIM
n<-result$n
bsp<-result$bsp

X_est <- seq(from = -15,to = 10, by = 0.1)    
l_X_est <- length(X_est)
X_est<-as.matrix(X_est)

values_est_naive <- matrix(ncol=SIM,nrow=l_X_est)
values_est_est1 <- matrix(ncol=SIM,nrow=l_X_est)
values_est_est1_ridge <- matrix(ncol=SIM,nrow=l_X_est)
values_est_est3 <- matrix(ncol=SIM,nrow=l_X_est)
values_est_est3_ridge <- matrix(ncol=SIM,nrow=l_X_est)


for (sim in 1:SIM) {
  
  print(sim)
  betahat_first <- as.vector(result$betahat_first_complete[,sim])
  betahat_first_ridge <- as.vector(result$betahat_first_ridge_complete[,sim])
  beta_end1 <- as.vector(result$beta_end1_complete[,sim])
  beta_end1_ridge <- as.vector(result$beta_end1_ridge_complete[,sim])
  naive_beta <- as.vector(result$naive_beta_complete[,sim])
  if(bsp==7|bsp==8){
    beta_ZlmX<-as.vector(result$beta_ZlmX_complete[,sim])
    beta_Z2lmX<-as.vector(result$beta_Z2lmX_complete[,sim])
    beta_Z3lmX<-as.vector(result$beta_Z3lmX_complete[,sim])
  }
  
  if(empty){Z<-result$Z_Ddata_complete[,sim]
  if(bsp==5|bsp==6|bsp==7|bsp==8){
    W<-result$W_Ddata_complete[,sim]}}
  else{  Z<-result$Z_complete[,sim]
  if(bsp==5|bsp==6|bsp==7|bsp==8){
    W<-result$W_complete[,sim]}}
  
  
  #E[Y|do(X)]
  #example 5
  if(bsp==7){
    original<-function(x){
      mean(result$y0_complete)+mean(result$yx1_complete)*x+mean(result$yx2_complete)*x^2+mean(result$yx3_complete)*x^3+mean(result$yz1_complete)*(-2*x)+mean(result$yz2_complete)*(4*x^2+1)+mean(result$yz3_complete)*(-8*x^3-6*x)+mean(result$yw1_complete)*(-1)+mean(result$yw2_complete)*(2)+mean(result$yw3_complete)*(-4)}}
  
  #example 6
  if(bsp==8){
    original<-function(x){
      mean(result$y0_complete)+mean(result$yx1_complete)*x+mean(result$yx2_complete)*x^2+mean(result$yx3_complete)*x^3+mean(result$yz1_complete)*(1.5*x-1)+mean(result$yz2_complete)*(x^2+1)+mean(result$yz3_complete)*(x^3+3*x)+mean(result$yw1_complete)*(2)+mean(result$yw2_complete)*(5)+mean(result$yw3_complete)*(14)}}
  
  #E[Y|X]
  #example 5
  if(bsp==7){
    wrong<-function(x){
      mean(result$y0_complete)+mean(result$yx1_complete)*x+mean(result$yx2_complete)*x^2+mean(result$yx3_complete)*x^3+mean(result$yz1_complete)*(-2*x)+mean(result$yz2_complete)*(4*x^2+1)+mean(result$yz3_complete)*(-8*x^3-6*x)+mean(result$yw1_complete)*(-0.5+0.5*x)+mean(result$yw2_complete)*(1-x)+mean(result$yw3_complete)*(-1+3*x)}}
  
  #example 6
  if(bsp==8){
    wrong<-function(x){
      mean(result$y0_complete)+mean(result$yx1_complete)*x+mean(result$yx2_complete)*x^2+mean(result$yx3_complete)*x^3+mean(result$yz1_complete)*(1.5*x-1)+mean(result$yz2_complete)*(x^2+1)+mean(result$yz3_complete)*(x^3+3*x)+mean(result$yw1_complete)*(1+0.5*x)+mean(result$yw2_complete)*(1+2*x)+mean(result$yw3_complete)*(7.5*x-1)}}
  
  #RR ridge
  est1_ridge<-function(x)
  {beta_end1_ridge[1]+beta_end1_ridge[2]*x+beta_end1_ridge[3]*x^2+beta_end1_ridge[4]*x^3}
  
  #RR OLS
  est1<-function(x)
  {beta_end1[1]+beta_end1[2]*x+beta_end1[3]*x^2+beta_end1[4]*x^3}
  
  #naive 
  naive<-function(x)
  {naive_beta[1]+naive_beta[2]*x+naive_beta[3]*x^2+naive_beta[4]*x^3}
  
  #TSR ridge
  if(bsp==7 | bsp==8){
    est3_ridge<-function(x){
      return(betahat_first_ridge[1]+betahat_first_ridge[2]*x+betahat_first_ridge[3]*x^2+betahat_first_ridge[4]*x^3+
               betahat_first_ridge[5]*(beta_ZlmX[1]+beta_ZlmX[2]*x+beta_ZlmX[3]*x^2+beta_ZlmX[4]*x^3)+
               betahat_first_ridge[6]*(beta_Z2lmX[1]+beta_Z2lmX[2]*x+beta_Z2lmX[3]*x^2+beta_Z2lmX[4]*x^3)+
               betahat_first_ridge[7]*(beta_Z3lmX[1]+beta_Z3lmX[2]*x+beta_Z3lmX[3]*x^2+beta_Z3lmX[4]*x^3)+
               betahat_first_ridge[8]*mean(W)+
               betahat_first_ridge[9]*mean(W^2)+
               betahat_first_ridge[10]*mean(W^3))}}
  
  #TSR OLS
  if(bsp==7 | bsp==8){
    est3<-function(x){
      return(betahat_first[1]+betahat_first[2]*x+betahat_first[3]*x^2+betahat_first[4]*x^3+
               betahat_first[5]*(beta_ZlmX[1]+beta_ZlmX[2]*x+beta_ZlmX[3]*x^2+beta_ZlmX[4]*x^3)+
               betahat_first[6]*(beta_Z2lmX[1]+beta_Z2lmX[2]*x+beta_Z2lmX[3]*x^2+beta_Z2lmX[4]*x^3)+
               betahat_first[7]*(beta_Z3lmX[1]+beta_Z3lmX[2]*x+beta_Z3lmX[3]*x^2+beta_Z3lmX[4]*x^3)+
               betahat_first[8]*mean(W)+
               betahat_first[9]*mean(W^2)+
               betahat_first[10]*mean(W^3))}}
  
  values_est_naive[,sim]<-apply(X_est,1,naive)
  values_est_est1[,sim]<-apply(X_est,1,est1)
  values_est_est1_ridge[,sim]<-apply(X_est,1,est1_ridge)
  values_est_est3[,sim]<-apply(X_est,1,est3)
  values_est_est3_ridge[,sim]<-apply(X_est,1,est3_ridge)
  
}
alpha<-0.05

#calculate lower and upper alpha/2 quantiles
q_lower_naive <- apply(values_est_naive, 1, function(col) quantile(col, alpha/2))
q_lower_est1 <- apply(values_est_est1, 1, function(col) quantile(col, alpha/2))
q_lower_est1_ridge <- apply(values_est_est1_ridge, 1, function(col) quantile(col, alpha/2))
q_lower_est3 <- apply(values_est_est3, 1, function(col) quantile(col, alpha/2))
q_lower_est3_ridge <- apply(values_est_est3_ridge, 1, function(col) quantile(col, alpha/2))

q_upper_naive <- apply(values_est_naive, 1, function(col) quantile(col, 1-alpha/2))
q_upper_est1 <- apply(values_est_est1, 1, function(col) quantile(col, 1-alpha/2))
q_upper_est1_ridge <- apply(values_est_est1_ridge, 1, function(col) quantile(col, 1-alpha/2))
q_upper_est3 <- apply(values_est_est3, 1, function(col) quantile(col, 1-alpha/2))
q_upper_est3_ridge <- apply(values_est_est3_ridge, 1, function(col) quantile(col, 1-alpha/2))


#for the boxplots in figure
if(empty){
  X <- as.vector(result$X_Ddata_complete)
  X_S <- as.vector(result$X_Sdata_complete)[as.vector(result$S_Sdata_complete)==1]
}else{X<-as.vector(result$X_complete)
X_S<-X[as.vector(result$S_complete)==1]}

est1_mean<-apply(values_est_est1,1,mean)
est1_ridge_mean<-apply(values_est_est1_ridge,1,mean)
est3_mean<-apply(values_est_est3,1,mean)
est3_ridge_mean<-apply(values_est_est3_ridge,1,mean)
naive_mean<-apply(values_est_naive,1,mean)




#plots with 95% areas


#example 5

#OLS
ggplot() +
  geom_ribbon(aes(x = X_est, ymin = q_lower_naive, ymax = q_upper_naive,fill="naive"),alpha=0.1) +
  geom_ribbon(aes(x = X_est, ymin = q_lower_est1, ymax = q_upper_est1,fill="RR"),alpha=0.8) +
  geom_ribbon(aes(x = X_est, ymin = q_lower_est3, ymax = q_upper_est3,fill="TSR"),alpha=0.8) +
  geom_line(aes(x = X_est, y = naive_mean, color = "naive"), size = 1.5)+
  geom_line(aes(x = X_est, y = q_lower_naive, color = "naive"),linetype="dashed", size = 1)+
  geom_line(aes(x = X_est, y = q_upper_naive, color = "naive"),linetype="dashed", size = 1)+
  geom_line(aes(x = X_est, y = est1_mean, color = "RR"), size = 1.5)+
  geom_line(aes(x = X_est, y = q_lower_est1, color = "RR"),linetype="dashed", size = 1)+
  geom_line(aes(x = X_est, y = q_upper_est1, color = "RR"),linetype="dashed", size = 1)+
  geom_line(aes(x = X_est, y = est3_mean, color = "TSR"), size = 1.5)+
  geom_line(aes(x = X_est, y = q_lower_est3, color = "TSR"),linetype="dashed", size = 1)+
  geom_line(aes(x = X_est, y = q_upper_est3, color = "TSR"),linetype="dashed", size = 1)+
  geom_line(aes(x = X_est, y = apply(X_est,1,wrong), color = "E[Y|X=x]"), linewidth = 1.3) +
  geom_line(aes(x = X_est, y = apply(X_est,1,original), color = "ground truth"), linewidth = 1.3) +
  theme_minimal() +
  coord_cartesian(xlim = c(-7, 5), ylim = c(-5, 3)) + 
  geom_boxplot(aes(x = X, y = 1.5), width = 0.5, fill = "grey", color = "black", alpha = 0.7) +
  geom_boxplot(aes(x = X_S, y = 2.5), width = 0.5, fill = "grey", color = "black", alpha = 0.7) +
  scale_color_manual(values = c("ground truth" = "black", "naive"="orange", "RR"="blue", "TSR"=middle_green, "E[Y|X=x]" = "purple"), name = NULL,breaks=c("ground truth","naive","RR","TSR","E[Y|X=x]")) +
  scale_fill_manual(values = c("naive" = "orange", "RR" = "lightblue", "TSR" = "lightgreen"), name = NULL) +
  labs(x="x",y=expression(hat(E) * group("[", Y * group("|", do(X == x), ""), "]")))+
  theme(
    axis.title.x = element_text(size = 30),  
    axis.title.y = element_text(size = 20),  
    legend.title = element_text(size = 20),  
    legend.text = element_text(size = 20),  
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.position = "right",
    legend.spacing = unit(1.5, "cm"),
    legend.key.height = unit(1.3, "cm"),  
    legend.key.width = unit(1.3, "cm")
  )+
  guides(
    fill  = guide_legend(order = 2),
    color = guide_legend(order = 1)
  )



#ridge
ggplot() +
  geom_ribbon(aes(x = X_est, ymin = q_lower_naive, ymax = q_upper_naive,fill="naive"),alpha=0.1) +
  geom_ribbon(aes(x = X_est, ymin = q_lower_est1_ridge, ymax = q_upper_est1_ridge,fill="RR (ridge)"),alpha=0.8) +
  geom_ribbon(aes(x = X_est, ymin = q_lower_est3_ridge, ymax = q_upper_est3_ridge,fill="TSR (ridge)"),alpha=0.8) +
  geom_line(aes(x = X_est, y = naive_mean, color = "naive"), size = 1.5)+
  geom_line(aes(x = X_est, y = q_lower_naive, color = "naive"),linetype="dashed", size = 1)+
  geom_line(aes(x = X_est, y = q_upper_naive, color = "naive"),linetype="dashed", size = 1)+
  geom_line(aes(x = X_est, y = est1_ridge_mean, color = "RR (ridge)"), size = 1.5)+
  geom_line(aes(x = X_est, y = q_lower_est1_ridge, color = "RR (ridge)"),linetype="dashed", size = 1)+
  geom_line(aes(x = X_est, y = q_upper_est1_ridge, color = "RR (ridge)"),linetype="dashed", size = 1)+
  geom_line(aes(x = X_est, y = est3_ridge_mean, color = "TSR (ridge)"), size = 1.5)+
  geom_line(aes(x = X_est, y = q_lower_est3_ridge, color = "TSR (ridge)"),linetype="dashed", size = 1)+
  geom_line(aes(x = X_est, y = q_upper_est3_ridge, color = "TSR (ridge)"),linetype="dashed", size = 1)+
  geom_line(aes(x = X_est, y = apply(X_est,1,wrong), color = "E[Y|X=x]"), linewidth = 1.3) +
  geom_line(aes(x = X_est, y = apply(X_est,1,original), color = "ground truth"), linewidth = 1.3) +
  theme_minimal() +
  coord_cartesian(xlim = c(-7, 5), ylim = c(-5, 3)) + 
  geom_boxplot(aes(x = X, y = 1.5), width = 0.5, fill = "grey", color = "black", alpha = 0.7) +
  geom_boxplot(aes(x = X_S, y = 2.5), width = 0.5, fill = "grey", color = "black", alpha = 0.7) +
  scale_color_manual(values = c("ground truth" = "black", "naive"="orange", "RR (ridge)"="blue", "TSR (ridge)"=middle_green, "E[Y|X=x]" = "purple"), name = NULL,breaks=c("ground truth","naive","RR (ridge)","TSR (ridge)","E[Y|X=x]")) +
  scale_fill_manual(values = c("naive" = "orange", "RR (ridge)" = "lightblue", "TSR (ridge)" = "lightgreen"), name = NULL) +
  labs(x="x",y=expression(hat(E) * group("[", Y * group("|", do(X == x), ""), "]")))+
  theme(
    axis.title.x = element_text(size = 30),  
    axis.title.y = element_text(size = 20),  
    legend.title = element_text(size = 20),  
    legend.text = element_text(size = 20),  
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.position = "none",
    legend.spacing = unit(1.5, "cm"),
    legend.key.height = unit(1.3, "cm"),  
    legend.key.width = unit(1.3, "cm")
  )



#example 6 

#OLS
ggplot() +
  geom_ribbon(aes(x = X_est, ymin = q_lower_naive, ymax = q_upper_naive,fill="naive"),alpha=0.1) +
  geom_ribbon(aes(x = X_est, ymin = q_lower_est1, ymax = q_upper_est1,fill="RR"),alpha=0.8) +
  geom_ribbon(aes(x = X_est, ymin = q_lower_est3, ymax = q_upper_est3,fill="TSR"),alpha=0.8) +
  geom_line(aes(x = X_est, y = naive_mean, color = "naive"), size = 1.5)+
  geom_line(aes(x = X_est, y = q_lower_naive, color = "naive"),linetype="dashed", size = 1)+
  geom_line(aes(x = X_est, y = q_upper_naive, color = "naive"),linetype="dashed", size = 1)+
  geom_line(aes(x = X_est, y = est1_mean, color = "RR"), size = 1.5)+
  geom_line(aes(x = X_est, y = q_lower_est1, color = "RR"),linetype="dashed", size = 1)+
  geom_line(aes(x = X_est, y = q_upper_est1, color = "RR"),linetype="dashed", size = 1)+
  geom_line(aes(x = X_est, y = est3_mean, color = "TSR"), size = 1.5)+
  geom_line(aes(x = X_est, y = q_lower_est3, color = "TSR"),linetype="dashed", size = 1)+
  geom_line(aes(x = X_est, y = q_upper_est3, color = "TSR"),linetype="dashed", size = 1)+
  geom_line(aes(x = X_est, y = apply(X_est,1,wrong), color = "E[Y|X=x]"), linewidth = 1.3) +
  geom_line(aes(x = X_est, y = apply(X_est,1,original), color = "ground truth"), linewidth = 1.3) +
  theme_minimal() +
  coord_cartesian(xlim = c(-2.5, 7), ylim = c(-2.5, 7.5)) + 
  geom_boxplot(aes(x = X, y = -2.5), width = 0.5, fill = "grey", color = "black", alpha = 0.7) +
  geom_boxplot(aes(x = X_S, y = -1.5), width = 0.5, fill = "grey", color = "black", alpha = 0.7) +
  scale_color_manual(values = c("ground truth" = "black", "naive"="orange", "RR"="blue", "TSR"=middle_green, "E[Y|X=x]" = "purple"), name = NULL,breaks=c("ground truth","naive","RR","TSR","E[Y|X=x]")) +
  scale_fill_manual(values = c("naive" = "orange", "RR" = "lightblue", "TSR" = "lightgreen"), name = NULL) +
  labs(x="x",y=expression(hat(E) * group("[", Y * group("|", do(X == x), ""), "]")))+
  theme(
    axis.title.x = element_text(size = 30),  
    axis.title.y = element_text(size = 20),  
    legend.title = element_text(size = 20),  
    legend.text = element_text(size = 20),  
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.position = "none",
    legend.spacing = unit(1.5, "cm"),
    legend.key.height = unit(1.3, "cm"),  
    legend.key.width = unit(1.3, "cm")
  )



#ridge
ggplot() +
  geom_ribbon(aes(x = X_est, ymin = q_lower_naive, ymax = q_upper_naive,fill="naive"),alpha=0.1) +
  geom_ribbon(aes(x = X_est, ymin = q_lower_est1_ridge, ymax = q_upper_est1_ridge,fill="RR (ridge)"),alpha=0.8) +
  geom_ribbon(aes(x = X_est, ymin = q_lower_est3_ridge, ymax = q_upper_est3_ridge,fill="TSR (ridge)"),alpha=0.8) +
  geom_line(aes(x = X_est, y = naive_mean, color = "naive"), size = 1.5)+
  geom_line(aes(x = X_est, y = q_lower_naive, color = "naive"),linetype="dashed", size = 1)+
  geom_line(aes(x = X_est, y = q_upper_naive, color = "naive"),linetype="dashed", size = 1)+
  geom_line(aes(x = X_est, y = est1_ridge_mean, color = "RR (ridge)"), size = 1.5)+
  geom_line(aes(x = X_est, y = q_lower_est1_ridge, color = "RR (ridge)"),linetype="dashed", size = 1)+
  geom_line(aes(x = X_est, y = q_upper_est1_ridge, color = "RR (ridge)"),linetype="dashed", size = 1)+
  geom_line(aes(x = X_est, y = est3_ridge_mean, color = "TSR (ridge)"), size = 1.5)+
  geom_line(aes(x = X_est, y = q_lower_est3_ridge, color = "TSR (ridge)"),linetype="dashed", size = 1)+
  geom_line(aes(x = X_est, y = q_upper_est3_ridge, color = "TSR (ridge)"),linetype="dashed", size = 1)+
  geom_line(aes(x = X_est, y = apply(X_est,1,wrong), color = "E[Y|X=x]"), linewidth = 1.3) +
  geom_line(aes(x = X_est, y = apply(X_est,1,original), color = "ground truth"), linewidth = 1.3) +
  theme_minimal() +
  coord_cartesian(xlim = c(-2.5, 7), ylim = c(-2.5, 7.5)) + 
  geom_boxplot(aes(x = X, y = -2.5), width = 0.5, fill = "grey", color = "black", alpha = 0.7) +
  geom_boxplot(aes(x = X_S, y = -1.5), width = 0.5, fill = "grey", color = "black", alpha = 0.7) +
  scale_color_manual(values = c("ground truth" = "black", "naive"="orange", "RR (ridge)"="blue", "TSR (ridge)"=middle_green, "E[Y|X=x]" = "purple"), name = NULL,breaks=c("ground truth","naive","RR (ridge)","TSR (ridge)","E[Y|X=x]")) +
  scale_fill_manual(values = c("naive" = "orange", "RR (ridge)" = "lightblue", "TSR (ridge)" = "lightgreen"), name = NULL) +
  labs(x="x",y=expression(hat(E) * group("[", Y * group("|", do(X == x), ""), "]")))+
  theme(
    axis.title.x = element_text(size = 30),  
    axis.title.y = element_text(size = 20),  
    legend.title = element_text(size = 20),  
    legend.text = element_text(size = 20),  
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.position = "none",
    legend.spacing = unit(1.5, "cm"),
    legend.key.height = unit(1.3, "cm"),  
    legend.key.width = unit(1.3, "cm")
  )






load("resultssim78.RData")  

####################################
#######      bsp7     ########
####################################



#MSE mean
c(mean(result500_bsp7_rand$o_test_naive_S), mean(result500_bsp7_rand$o_test_est1_S), mean(result500_bsp7_rand$o_test_est1_ridge_S), mean(result500_bsp7_rand$o_test_est3_S), mean(result500_bsp7_rand$o_test_est3_ridge_S))
#[1] 29.252106  8.946369  8.939645  3.946392  3.885658
c(mean(result1000_bsp7_rand$o_test_naive_S), mean(result1000_bsp7_rand$o_test_est1_S), mean(result1000_bsp7_rand$o_test_est1_ridge_S), mean(result1000_bsp7_rand$o_test_est3_S), mean(result1000_bsp7_rand$o_test_est3_ridge_S))
#[1] 20.196436  7.346347  7.322737  1.859016  1.832987
c(mean(result5000_bsp7_rand$o_test_naive_S), mean(result5000_bsp7_rand$o_test_est1_S), mean(result5000_bsp7_rand$o_test_est1_ridge_S), mean(result5000_bsp7_rand$o_test_est3_S), mean(result5000_bsp7_rand$o_test_est3_ridge_S))
#[1] 15.5274680  4.5279745  4.5187518  0.2906700  0.2861789

c(mean(result500_bsp7_rand$o_test_naive), mean(result500_bsp7_rand$o_test_est1), mean(result500_bsp7_rand$o_test_est1_ridge), mean(result500_bsp7_rand$o_test_est3), mean(result500_bsp7_rand$o_test_est3_ridge))
#[1] 645.27085  41.67552  42.43579  24.91233  25.90613
c(mean(result1000_bsp7_rand$o_test_naive), mean(result1000_bsp7_rand$o_test_est1), mean(result1000_bsp7_rand$o_test_est1_ridge), mean(result1000_bsp7_rand$o_test_est3), mean(result1000_bsp7_rand$o_test_est3_ridge))
#[1] 458.251444  23.233841  23.071146   8.273629   7.990351
c(mean(result5000_bsp7_rand$o_test_naive), mean(result5000_bsp7_rand$o_test_est1), mean(result5000_bsp7_rand$o_test_est1_ridge), mean(result5000_bsp7_rand$o_test_est3), mean(result5000_bsp7_rand$o_test_est3_ridge))
#[1] 285.395493  13.128223  13.066204   1.633540   1.650222

#MSE sd
c(sd(result500_bsp7_rand$o_test_naive_S), sd(result500_bsp7_rand$o_test_est1_S), sd(result500_bsp7_rand$o_test_est1_ridge_S), sd(result500_bsp7_rand$o_test_est3_S), sd(result500_bsp7_rand$o_test_est3_ridge_S))
#[1] 54.107934  8.364460  8.337654  6.969919  6.736784
c(sd(result1000_bsp7_rand$o_test_naive_S), sd(result1000_bsp7_rand$o_test_est1_S), sd(result1000_bsp7_rand$o_test_est1_ridge_S), sd(result1000_bsp7_rand$o_test_est3_S), sd(result1000_bsp7_rand$o_test_est3_ridge_S))
#[1] 20.623851  7.378685  7.335513  3.900902  3.861521
c(sd(result5000_bsp7_rand$o_test_naive_S), sd(result5000_bsp7_rand$o_test_est1_S), sd(result5000_bsp7_rand$o_test_est1_ridge_S), sd(result5000_bsp7_rand$o_test_est3_S), sd(result5000_bsp7_rand$o_test_est3_ridge_S))
#[1] 12.2115530  3.9490944  3.9376182  0.4233257  0.4158697

c(sd(result500_bsp7_rand$o_test_naive), sd(result500_bsp7_rand$o_test_est1), sd(result500_bsp7_rand$o_test_est1_ridge), sd(result500_bsp7_rand$o_test_est3), sd(result500_bsp7_rand$o_test_est3_ridge))
#[1] 846.81120  64.15996  60.81413  40.01128  40.29705
c(sd(result1000_bsp7_rand$o_test_naive), sd(result1000_bsp7_rand$o_test_est1), sd(result1000_bsp7_rand$o_test_est1_ridge), sd(result1000_bsp7_rand$o_test_est3), sd(result1000_bsp7_rand$o_test_est3_ridge))
#[1] 683.79802  26.25120  25.11390  12.43129  11.58191
c(sd(result5000_bsp7_rand$o_test_naive), sd(result5000_bsp7_rand$o_test_est1), sd(result5000_bsp7_rand$o_test_est1_ridge), sd(result5000_bsp7_rand$o_test_est3), sd(result5000_bsp7_rand$o_test_est3_ridge))
#[1] 349.377241  12.245697  12.114180   2.222432   2.264046


#################
###   empty   ###
#################



#MSE mean
c(mean(empty_result500_bsp7_rand$o_test_naive_S), mean(empty_result500_bsp7_rand$o_test_est1_S), mean(empty_result500_bsp7_rand$o_test_est1_ridge_S), mean(empty_result500_bsp7_rand$o_test_est3_S), mean(empty_result500_bsp7_rand$o_test_est3_ridge_S))
#[1] 26.439543  7.758931  7.582983  2.547173  2.396763
c(mean(empty_result1000_bsp7_rand$o_test_naive_S), mean(empty_result1000_bsp7_rand$o_test_est1_S), mean(empty_result1000_bsp7_rand$o_test_est1_ridge_S), mean(empty_result1000_bsp7_rand$o_test_est3_S), mean(empty_result1000_bsp7_rand$o_test_est3_ridge_S))
#[1] 18.144064  6.404541  6.377790  1.067684  1.051495
c(mean(empty_result5000_bsp7_rand$o_test_naive_S), mean(empty_result5000_bsp7_rand$o_test_est1_S), mean(empty_result5000_bsp7_rand$o_test_est1_ridge_S), mean(empty_result5000_bsp7_rand$o_test_est3_S), mean(empty_result5000_bsp7_rand$o_test_est3_ridge_S))
#[1] 16.2940583  4.3844718  4.3781812  0.2539810  0.2534483

c(mean(empty_result500_bsp7_rand$o_test_naive), mean(empty_result500_bsp7_rand$o_test_est1), mean(empty_result500_bsp7_rand$o_test_est1_ridge), mean(empty_result500_bsp7_rand$o_test_est3), mean(empty_result500_bsp7_rand$o_test_est3_ridge))
#[1] 685.67459  41.06983  40.47980  25.68441  25.37988
c(mean(empty_result1000_bsp7_rand$o_test_naive), mean(empty_result1000_bsp7_rand$o_test_est1), mean(empty_result1000_bsp7_rand$o_test_est1_ridge), mean(empty_result1000_bsp7_rand$o_test_est3), mean(empty_result1000_bsp7_rand$o_test_est3_ridge))
#[1] 447.968835  19.710630  19.280710   7.141670   6.811159
c(mean(empty_result5000_bsp7_rand$o_test_naive), mean(empty_result5000_bsp7_rand$o_test_est1), mean(empty_result5000_bsp7_rand$o_test_est1_ridge), mean(empty_result5000_bsp7_rand$o_test_est3), mean(empty_result5000_bsp7_rand$o_test_est3_ridge))
#[1] 282.292355  12.890830  12.785456   1.429533   1.409451

#MSE sd
c(sd(empty_result500_bsp7_rand$o_test_naive_S), sd(empty_result500_bsp7_rand$o_test_est1_S), sd(empty_result500_bsp7_rand$o_test_est1_ridge_S), sd(empty_result500_bsp7_rand$o_test_est3_S), sd(empty_result500_bsp7_rand$o_test_est3_ridge_S))
#[1] 44.312546  7.946311  7.130953  5.463246  4.596399
c(sd(empty_result1000_bsp7_rand$o_test_naive_S), sd(empty_result1000_bsp7_rand$o_test_est1_S), sd(empty_result1000_bsp7_rand$o_test_est1_ridge_S), sd(empty_result1000_bsp7_rand$o_test_est3_S), sd(empty_result1000_bsp7_rand$o_test_est3_ridge_S))
#[1] 16.378746  4.631317  4.608361  1.634010  1.603251
c(sd(empty_result5000_bsp7_rand$o_test_naive_S), sd(empty_result5000_bsp7_rand$o_test_est1_S), sd(empty_result5000_bsp7_rand$o_test_est1_ridge_S), sd(empty_result5000_bsp7_rand$o_test_est3_S), sd(empty_result5000_bsp7_rand$o_test_est3_ridge_S))
#[1] 14.9919367  3.7623828  3.7551229  0.3074808  0.3091415

c(sd(empty_result500_bsp7_rand$o_test_naive), sd(empty_result500_bsp7_rand$o_test_est1), sd(empty_result500_bsp7_rand$o_test_est1_ridge), sd(empty_result500_bsp7_rand$o_test_est3), sd(empty_result500_bsp7_rand$o_test_est3_ridge))
#[1] 957.57864  66.24976  59.89052  51.83786  46.43091
c(sd(empty_result1000_bsp7_rand$o_test_naive), sd(empty_result1000_bsp7_rand$o_test_est1), sd(empty_result1000_bsp7_rand$o_test_est1_ridge), sd(empty_result1000_bsp7_rand$o_test_est3), sd(empty_result1000_bsp7_rand$o_test_est3_ridge))
#[1] 651.40973  20.32105  20.05145  11.49143  11.30025
c(sd(empty_result5000_bsp7_rand$o_test_naive), sd(empty_result5000_bsp7_rand$o_test_est1), sd(empty_result5000_bsp7_rand$o_test_est1_ridge), sd(empty_result5000_bsp7_rand$o_test_est3), sd(empty_result5000_bsp7_rand$o_test_est3_ridge))
#[1] 344.475262  11.285491  11.082094   1.893172   1.912072






####################################
#######      bsp8     ########
####################################



#MSE mean
c(mean(result500_bsp8_rand$o_test_naive_S), mean(result500_bsp8_rand$o_test_est1_S), mean(result500_bsp8_rand$o_test_est1_ridge_S), mean(result500_bsp8_rand$o_test_est3_S), mean(result500_bsp8_rand$o_test_est3_ridge_S))
#[1] 13.162457  6.351302  6.329774  0.461770  0.441232
c(mean(result1000_bsp8_rand$o_test_naive_S), mean(result1000_bsp8_rand$o_test_est1_S), mean(result1000_bsp8_rand$o_test_est1_ridge_S), mean(result1000_bsp8_rand$o_test_est3_S), mean(result1000_bsp8_rand$o_test_est3_ridge_S))
#[1] 13.1142497  6.3197155  6.3128365  0.2892735  0.2827083
c(mean(result5000_bsp8_rand$o_test_naive_S), mean(result5000_bsp8_rand$o_test_est1_S), mean(result5000_bsp8_rand$o_test_est1_ridge_S), mean(result5000_bsp8_rand$o_test_est3_S), mean(result5000_bsp8_rand$o_test_est3_ridge_S))
#[1] 10.77895187  4.49509829  4.49465821  0.07681068  0.07628328

c(mean(result500_bsp8_rand$o_test_naive), mean(result500_bsp8_rand$o_test_est1), mean(result500_bsp8_rand$o_test_est1_ridge), mean(result500_bsp8_rand$o_test_est3), mean(result500_bsp8_rand$o_test_est3_ridge))
#[1] 261.503710  17.591790  17.357270   4.119865   3.880787
c(mean(result1000_bsp8_rand$o_test_naive), mean(result1000_bsp8_rand$o_test_est1), mean(result1000_bsp8_rand$o_test_est1_ridge), mean(result1000_bsp8_rand$o_test_est3), mean(result1000_bsp8_rand$o_test_est3_ridge))
#[1] 252.639700  14.225791  14.011844   2.498192   2.289504
c(mean(result5000_bsp8_rand$o_test_naive), mean(result5000_bsp8_rand$o_test_est1), mean(result5000_bsp8_rand$o_test_est1_ridge), mean(result5000_bsp8_rand$o_test_est3), mean(result5000_bsp8_rand$o_test_est3_ridge))
#[1] 216.2568038   9.6836539   9.7006043   0.4168671   0.4367547


#MSE sd
c(sd(result500_bsp8_rand$o_test_naive_S), sd(result500_bsp8_rand$o_test_est1_S), sd(result500_bsp8_rand$o_test_est1_ridge_S), sd(result500_bsp8_rand$o_test_est3_S), sd(result500_bsp8_rand$o_test_est3_ridge_S))
#[1] 10.0612663  6.3778463  6.3724054  0.5125688  0.4699172
c(sd(result1000_bsp8_rand$o_test_naive_S), sd(result1000_bsp8_rand$o_test_est1_S), sd(result1000_bsp8_rand$o_test_est1_ridge_S), sd(result1000_bsp8_rand$o_test_est3_S), sd(result1000_bsp8_rand$o_test_est3_ridge_S))
#[1] 8.6230150 6.1118837 6.1099756 0.3522333 0.3475407
c(sd(result5000_bsp8_rand$o_test_naive_S), sd(result5000_bsp8_rand$o_test_est1_S), sd(result5000_bsp8_rand$o_test_est1_ridge_S), sd(result5000_bsp8_rand$o_test_est3_S), sd(result5000_bsp8_rand$o_test_est3_ridge_S))
#[1] 8.5207183 4.8067858 4.8058267 0.1553543 0.1508171

c(sd(result500_bsp8_rand$o_test_naive), sd(result500_bsp8_rand$o_test_est1), sd(result500_bsp8_rand$o_test_est1_ridge), sd(result500_bsp8_rand$o_test_est3), sd(result500_bsp8_rand$o_test_est3_ridge))
#[1] 304.631626  18.277861  19.084103   5.699294   5.872648
c(sd(result1000_bsp8_rand$o_test_naive), sd(result1000_bsp8_rand$o_test_est1), sd(result1000_bsp8_rand$o_test_est1_ridge), sd(result1000_bsp8_rand$o_test_est3), sd(result1000_bsp8_rand$o_test_est3_ridge))
#[1] 258.013131  14.102617  13.759804   3.771183   3.409658
c(sd(result5000_bsp8_rand$o_test_naive), sd(result5000_bsp8_rand$o_test_est1), sd(result5000_bsp8_rand$o_test_est1_ridge), sd(result5000_bsp8_rand$o_test_est3), sd(result5000_bsp8_rand$o_test_est3_ridge))
#[1] 231.0222626  10.0968492  10.1436165   0.5525161   0.5201834



#################
###   empty   ###
#################



#MSE mean
c(mean(empty_result500_bsp8_rand$o_test_naive_S), mean(empty_result500_bsp8_rand$o_test_est1_S), mean(empty_result500_bsp8_rand$o_test_est1_ridge_S), mean(empty_result500_bsp8_rand$o_test_est3_S), mean(empty_result500_bsp8_rand$o_test_est3_ridge_S))
#[1] 12.8477476  6.6133272  6.5865714  0.6155059  0.5973363
c(mean(empty_result1000_bsp8_rand$o_test_naive_S), mean(empty_result1000_bsp8_rand$o_test_est1_S), mean(empty_result1000_bsp8_rand$o_test_est1_ridge_S), mean(empty_result1000_bsp8_rand$o_test_est3_S), mean(empty_result1000_bsp8_rand$o_test_est3_ridge_S))
#[1] 13.4374476  6.3309732  6.3178684  0.2453620  0.2375485
c(mean(empty_result5000_bsp8_rand$o_test_naive_S), mean(empty_result5000_bsp8_rand$o_test_est1_S), mean(empty_result5000_bsp8_rand$o_test_est1_ridge_S), mean(empty_result5000_bsp8_rand$o_test_est3_S), mean(empty_result5000_bsp8_rand$o_test_est3_ridge_S))
#[1] 10.62719146  4.54521498  4.54545781  0.06466543  0.06480277

c(mean(empty_result500_bsp8_rand$o_test_naive), mean(empty_result500_bsp8_rand$o_test_est1), mean(empty_result500_bsp8_rand$o_test_est1_ridge), mean(empty_result500_bsp8_rand$o_test_est3), mean(empty_result500_bsp8_rand$o_test_est3_ridge))
#[1] 245.658188  17.018149  16.786728   4.054963   4.235738
c(mean(empty_result1000_bsp8_rand$o_test_naive), mean(empty_result1000_bsp8_rand$o_test_est1), mean(empty_result1000_bsp8_rand$o_test_est1_ridge), mean(empty_result1000_bsp8_rand$o_test_est3), mean(empty_result1000_bsp8_rand$o_test_est3_ridge))
#[1] 264.854348  14.593849  14.362011   2.168091   2.058049
c(mean(empty_result5000_bsp8_rand$o_test_naive), mean(empty_result5000_bsp8_rand$o_test_est1), mean(empty_result5000_bsp8_rand$o_test_est1_ridge), mean(empty_result5000_bsp8_rand$o_test_est3), mean(empty_result5000_bsp8_rand$o_test_est3_ridge))
#[1] 214.4288354   9.5649833   9.6087282   0.3690637   0.4043703

#MSE sd
c(sd(empty_result500_bsp8_rand$o_test_naive_S), sd(empty_result500_bsp8_rand$o_test_est1_S), sd(empty_result500_bsp8_rand$o_test_est1_ridge_S), sd(empty_result500_bsp8_rand$o_test_est3_S), sd(empty_result500_bsp8_rand$o_test_est3_ridge_S))
#[1] 9.300982 6.381525 6.373992 1.066914 1.025294
c(sd(empty_result1000_bsp8_rand$o_test_naive_S), sd(empty_result1000_bsp8_rand$o_test_est1_S), sd(empty_result1000_bsp8_rand$o_test_est1_ridge_S), sd(empty_result1000_bsp8_rand$o_test_est3_S), sd(empty_result1000_bsp8_rand$o_test_est3_ridge_S))
#[1] 8.7624497 6.0855671 6.0783913 0.2148091 0.2103457
c(sd(empty_result5000_bsp8_rand$o_test_naive_S), sd(empty_result5000_bsp8_rand$o_test_est1_S), sd(empty_result5000_bsp8_rand$o_test_est1_ridge_S), sd(empty_result5000_bsp8_rand$o_test_est3_S), sd(empty_result5000_bsp8_rand$o_test_est3_ridge_S))
#[1] 8.31044002 4.88087158 4.87980806 0.06312979 0.06324453

c(sd(empty_result500_bsp8_rand$o_test_naive), sd(empty_result500_bsp8_rand$o_test_est1), sd(empty_result500_bsp8_rand$o_test_est1_ridge), sd(empty_result500_bsp8_rand$o_test_est3), sd(empty_result500_bsp8_rand$o_test_est3_ridge))
#[1] 259.976668  17.836479  18.680712   5.400978   5.684730
c(sd(empty_result1000_bsp8_rand$o_test_naive), sd(empty_result1000_bsp8_rand$o_test_est1), sd(empty_result1000_bsp8_rand$o_test_est1_ridge), sd(empty_result1000_bsp8_rand$o_test_est3), sd(empty_result1000_bsp8_rand$o_test_est3_ridge))
#[1] 267.889031  15.443119  14.893994   3.187866   3.058300
c(sd(empty_result5000_bsp8_rand$o_test_naive), sd(empty_result5000_bsp8_rand$o_test_est1), sd(empty_result5000_bsp8_rand$o_test_est1_ridge), sd(empty_result5000_bsp8_rand$o_test_est3), sd(empty_result5000_bsp8_rand$o_test_est3_ridge))
#[1] 228.3872055  10.6673252  10.6557855   0.5962631   0.6352241


















































