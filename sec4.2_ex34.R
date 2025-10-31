#Section 4.2

#examples 3 & 4

#load packages
library(glmnet)
library(ggplot2)

#n = sample size
#SIM = number of simulation runs
#bsp = example number
#      note that bsp i corresponds to example i-2 (i.e. bsp3 corresponds to example 1)
#empty = indicates if ((S cap D)=emptyset)
#random = indicates uniformly distributed coefficients

sim34<-function(n,SIM,bsp,empty,random)
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
  
  for(sim in 1:SIM)
  {
    print(sim)
    
    #reproducibility
    set.seed(sim)
    
    #data generating processes
    #example 3
    if(bsp==5){
      if(empty==F){
        W<-rnorm(mean=2,sd=0.3,n=n)
        Z<-rnorm(mean=-0.3,n=n)
        X<-rnorm(mean=0,sd=1,n=n)+W
        S<-ifelse(Z>0 & X<2.5,1,0)
        if(random==F){
          y0<-0
          yx1<-0
          yx2<-0
          yx3<-0.1
          yz1<-0
          yz2<-0.1
          yz3<-0
          yw1<-0.2
          yw2<-0
          yw3<-0.5
        }else{
          y0<-runif(1)-0.5
          yx1<-runif(1)-0.5
          yx2<-runif(1)-0.5
          yx3<-runif(1)-0.4
          yz1<-runif(1)-0.5
          yz2<-runif(1)-0.4
          yz3<-runif(1)-0.5
          yw1<-runif(1)-0.3
          yw2<-runif(1)-0.5
          yw3<-runif(1)
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
        Y<-y0+yx1*X+yx2*X^2+yx3*X^3+yz1*Z+yz2*Z^2+yz3*Z^3+yw1*W+yw2*W^2+yw3*W^3+rnorm(n,sd=0.3)  
        
        Z_S<-Z[S==1]
        W_S<-W[S==1]
        X_S<-X[S==1]
        Y_S<-Y[S==1]
      }else{
        W_Sdata<-rnorm(mean=2,sd=0.3,n=n)
        Z_Sdata<-rnorm(mean=-0.3,n=n)
        X_Sdata<-rnorm(mean=0,sd=1,n=n)+W_Sdata
        S_Sdata<-ifelse(Z_Sdata>0 & X_Sdata<2.5,1,0)
        if(random==F){
          y0<-0
          yx1<-0
          yx2<-0
          yx3<-0.1
          yz1<-0
          yz2<-0.1
          yz3<-0
          yw1<-0.2
          yw2<-0
          yw3<-0.5
        }else{
          y0<-runif(1)-0.5
          yx1<-runif(1)-0.5
          yx2<-runif(1)-0.5
          yx3<-runif(1)-0.4
          yz1<-runif(1)-0.5
          yz2<-runif(1)-0.4
          yz3<-runif(1)-0.5
          yw1<-runif(1)-0.3
          yw2<-runif(1)-0.5
          yw3<-runif(1)
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
        Y_Sdata<-y0+yx1*X_Sdata+yx2*X_Sdata^2+yx3*X_Sdata^3+yz1*Z_Sdata+yz2*Z_Sdata^2+yz3*Z_Sdata^3+yw1*W_Sdata+yw2*W_Sdata^2+yw3*W_Sdata^3+rnorm(n,sd=0.3)  
        
        Z_S<-Z_Sdata[S_Sdata==1]
        W_S<-W_Sdata[S_Sdata==1]
        X_S<-X_Sdata[S_Sdata==1]
        Y_S<-Y_Sdata[S_Sdata==1]
        
        W_Ddata<-rnorm(mean=2,sd=0.3,n=n)
        Z_Ddata<-rnorm(mean=-0.3,n=n)
        X_Ddata<-rnorm(mean=0,sd=1,n=n)+W_Ddata
        
        X<-X_Ddata
        Z<-Z_Ddata
        W<-W_Ddata}}
    #example 4
    if(bsp==6){
      if(empty==F){
        W<-rnorm(mean=2,sd=0.3,n=n)
        X<-rnorm(n=n)+W
        Z<-rnorm(n)
        S<-rbinom(p=1/((1+exp(X))*(1+exp(Z))),size=1,n)
        if(random==F){
          y0<-0
          yx1<-0.5
          yx2<-0
          yx3<-0
          yz1<-1
          yz2<-0
          yz3<-0
          yw1<-3
          yw2<-0
          yw3<-0
        }else{
          y0<-runif(1)-0.5
          yx1<-runif(1)
          yx2<-runif(1)-0.5
          yx3<-runif(1)-0.5
          yz1<-runif(1)+0.5
          yz2<-runif(1)-0.5
          yz3<-runif(1)-0.5
          yw1<-runif(1)+2.5
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
        Y<-y0+yx1*X+yx2*X^2+yx3*X^3+yz1*Z+yz2*Z^2+yz3*Z^3+yw1*W+yw2*W^2+yw2*W^3+rnorm(n)  
        
        Z_S<-Z[S==1]
        W_S<-W[S==1]
        X_S<-X[S==1]
        Y_S<-Y[S==1]
      }else{
        W_Sdata<-rnorm(mean=2,sd=0.3,n=n)
        X_Sdata<-rnorm(n=n)+W_Sdata
        Z_Sdata<-rnorm(n)
        S_Sdata<-rbinom(p=1/((1+exp(X_Sdata))*(1+exp(Z_Sdata))),size=1,n)
          if(random==F){
            y0<-0
            yx1<-0.5
            yx2<-0
            yx3<-0
            yz1<-1
            yz2<-0
            yz3<-0
            yw1<-3
            yw2<-0
            yw3<-0
          }else{
            y0<-runif(1)-0.5
            yx1<-runif(1)
            yx2<-runif(1)-0.5
            yx3<-runif(1)-0.5
            yz1<-runif(1)+0.5
            yz2<-runif(1)-0.5
            yz3<-runif(1)-0.5
            yw1<-runif(1)+2.5
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
          
          W_Ddata<-rnorm(mean=2,sd=0.3,n=n)
          X_Ddata<-rnorm(n=n)+W_Ddata
          Z_Ddata<-rnorm(n)
          
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
      
      #E[Y|do(X)]
      #example 3
      if(bsp==5){
        original<-function(x){
          y0+yx1*x+yx2*x^2+yx3*x^3+yz1*(-0.3)+yz2*(1.09)+yz3*(-0.927)+yw1*(2)+yw2*(4.09)+yw3*(8.54)}}
      
      #example 4
      if(bsp==6){
        original<-function(x){
          y0+yx1*x+yx2*x^2+yx3*x^3+yz1*(0)+yz2*(1)+yz3*(0)+yw1*(2)+yw2*(4.09)+yw3*(8.54)}}
      
      #E[Y|X]
      #example 3
      if(bsp==5){
        wrong<-function(x){
          y0+yx1*x+yx2*x^2+yx3*x^3+yz1*(-0.3)+yz2*(1.09)+yz3*(-0.927)+yw1*(2+0.3^2/(0.3^2+1)*(x-2))+yw2*(4.09+(8.54-8-18)/1.09*(x-2))+yw3*(8.54+(18.1843-17.08)/1.09*(x-2))}}
      
      #example 4
      if(bsp==6){
        wrong<-function(x){
          y0+yx1*x+yx2*x^2+yx3*x^3+yz1*(0)+yz2*(1)+yz3*(0)+yw1*(2+(0.3^2)/(0.3^2+1)*(x-2))+yw2*(4.09+(8.54-8.18)/(0.3^2+1)*(x-2))+yw3*(8.54+(18.1843-17.08)/(0.3^2+1)*(x-2))}}
      
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
      if(bsp==5 | bsp==6){
        est3_ridge<-function(x){
          return(betahat_first_ridge[1]+betahat_first_ridge[2]*x+betahat_first_ridge[3]*x^2+betahat_first_ridge[4]*x^3+
                   betahat_first_ridge[5]*mean(Z)+
                   betahat_first_ridge[6]*mean(Z^2)+
                   betahat_first_ridge[7]*mean(Z^3)+
                   betahat_first_ridge[8]*mean(W)+
                   betahat_first_ridge[9]*mean(W^2)+
                   betahat_first_ridge[10]*mean(W^3))}}
      
      #TSR OLS
      if(bsp==5 |bsp==6){
        est3<-function(x){
          return(betahat_first[1]+betahat_first[2]*x+betahat_first[3]*x^2+betahat_first[4]*x^3+
                   betahat_first[5]*mean(Z)+
                   betahat_first[6]*mean(Z^2)+
                   betahat_first[7]*mean(Z^3)+
                   betahat_first[8]*mean(W)+
                   betahat_first[9]*mean(W^2)+
                   betahat_first[10]*mean(W^3))}}
      
      
      #save coefficient vectors from all simulation runs
      betahat_first_complete[,sim] <- betahat_first
      betahat_first_ridge_complete[,sim] <- betahat_first_ridge
      beta_end1_complete[,sim] <- beta_end1
      beta_end1_ridge_complete[,sim] <- beta_end1_ridge
      naive_beta_complete[,sim] <- naive_beta
  
      if(random==T){
      #MSE (w(rong): E[Y|X], o(riginal): E[Y|do(X)]) 
      
        if(bsp==5){
          W_test<-rnorm(mean=2,sd=0.3,n=n)
          Z_test<-rnorm(mean=-0.3,n=n)
          X_test<-as.matrix(rnorm(mean=0,sd=1,n=n)+W_test)
          S_test<-ifelse(Z_test>0 & X_test<2.5,1,0)
          Y_test<-y0+yx1*X_test+yx2*X_test^2+yx3*X_test^3+yz1*Z_test+yz2*Z_test^2+yz3*Z_test^3+yw1*W_test+yw2*W_test^2+yw3*W_test^3    
          X_test_S<-as.matrix(X_test[S_test==1,])
          Y_test_S<-Y_test[S_test==1]}
        if(bsp==6){
          W_test<-rnorm(mean=2,sd=0.3,n=n)
          X_test<-as.matrix(rnorm(n=n)+W_test)
          Z_test<-rnorm(n)
          S_test<-rbinom(p=1/((1+exp(X_test))*(1+exp(Z_test))),size=1,n)
          Y_test<-y0+yx1*X_test+yx2*X_test^2+yx3*X_test^3+yz1*Z_test+yz2*Z_test^2+yz3*Z_test^3+yw1*W_test+yw2*W_test^2+yw3*W_test^3    
          X_test_S<-as.matrix(X_test[S_test==1,])
          Y_test_S<-Y_test[S_test==1]}
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
    if(bsp==5 | bsp==6){
      if(random==F){
      if(empty==F){
        return(list(SIM = SIM, n = n, bsp = bsp,empty=empty,
                    X_complete = X_complete, Z_complete = Z_complete,W_complete = W_complete, S_complete = S_complete, Y_complete = Y_complete,
                    betahat_first_complete = betahat_first_complete, betahat_first_ridge_complete = betahat_first_ridge_complete,
                    beta_end1_complete = beta_end1_complete, beta_end1_ridge_complete = beta_end1_ridge_complete,
                    naive_beta_complete = naive_beta_complete,
                    y0_complete=y0_complete,yx1_complete=yx1_complete,yx2_complete=yx2_complete,yx3_complete=yx3_complete,yz1_complete=yz1_complete,yz2_complete=yz2_complete,yz3_complete=yz3_complete,yw1_complete=yw1_complete,yw2_complete=yw2_complete,yw3_complete=yw3_complete))
      }else{
        return(list(SIM = SIM, n = n, bsp = bsp,
                    X_Sdata_complete = X_Sdata_complete, Z_Sdata_complete = Z_Sdata_complete, W_Sdata_complete = W_Sdata_complete, S_Sdata_complete = S_Sdata_complete, Y_Sdata_complete = Y_Sdata_complete,
                    X_Ddata_complete = X_Ddata_complete, Z_Ddata_complete = Z_Ddata_complete, W_Ddata_complete = W_Ddata_complete,
                    betahat_first_complete = betahat_first_complete, betahat_first_ridge_complete = betahat_first_ridge_complete,
                    beta_end1_complete = beta_end1_complete, beta_end1_ridge_complete = beta_end1_ridge_complete,
                    naive_beta_complete = naive_beta_complete,
                    y0_complete=y0_complete,yx1_complete=yx1_complete,yx2_complete=yx2_complete,yx3_complete=yx3_complete,yz1_complete=yz1_complete,yz2_complete=yz2_complete,yz3_complete=yz3_complete,yw1_complete=yw1_complete,yw2_complete=yw2_complete,yw3_complete=yw3_complete))
        
      }}else{
        if(empty==F){
          return(list(SIM = SIM, n = n, bsp = bsp,empty=empty,
                      X_complete = X_complete, Z_complete = Z_complete,W_complete = W_complete, S_complete = S_complete, Y_complete = Y_complete,
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
                      y0_complete=y0_complete,yx1_complete=yx1_complete,yx2_complete=yx2_complete,yx3_complete=yx3_complete,yz1_complete=yz1_complete,yz2_complete=yz2_complete,yz3_complete=yz3_complete,yw1_complete=yw1_complete,yw2_complete=yw2_complete,yw3_complete=yw3_complete))
        }else{
          return(list(SIM = SIM, n = n, bsp = bsp,
                      X_Sdata_complete = X_Sdata_complete, Z_Sdata_complete = Z_Sdata_complete, W_Sdata_complete = W_Sdata_complete, S_Sdata_complete = S_Sdata_complete, Y_Sdata_complete = Y_Sdata_complete,
                      X_Ddata_complete = X_Ddata_complete, Z_Ddata_complete = Z_Ddata_complete, W_Ddata_complete = W_Ddata_complete,
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
                      y0_complete=y0_complete,yx1_complete=yx1_complete,yx2_complete=yx2_complete,yx3_complete=yx3_complete,yz1_complete=yz1_complete,yz2_complete=yz2_complete,yz3_complete=yz3_complete,yw1_complete=yw1_complete,yw2_complete=yw2_complete,yw3_complete=yw3_complete))
          
        }
      }
      
      
      }
   
  }

  
  
  
  #####################################################################################################
  #####################################################################################################
  #####################################################################################################
  
  #subset
  
  result500_bsp5 <- sim34(n=500,SIM=100,bsp=5,empty=F,random=F)
  result1000_bsp5 <- sim34(n=1000,SIM=100,bsp=5,empty=F,random=F)
  result5000_bsp5 <- sim34(n=5000,SIM=100,bsp=5,empty=F,random=F)
  
  result500_bsp5_rand <- sim34(n=500,SIM=100,bsp=5,empty=F,random=T)
  result1000_bsp5_rand <- sim34(n=1000,SIM=100,bsp=5,empty=F,random=T)
  result5000_bsp5_rand <- sim34(n=5000,SIM=100,bsp=5,empty=F,random=T)
  
  result500_bsp6 <- sim34(n=500,SIM=100,bsp=6,empty=F,random=F)
  result1000_bsp6 <- sim34(n=1000,SIM=100,bsp=6,empty=F,random=F)
  result5000_bsp6 <- sim34(n=5000,SIM=100,bsp=6,empty=F,random=F)
  
  result500_bsp6_rand <- sim34(n=500,SIM=100,bsp=6,empty=F,random=T)
  result1000_bsp6_rand <- sim34(n=1000,SIM=100,bsp=6,empty=F,random=T)
  result5000_bsp6_rand <- sim34(n=5000,SIM=100,bsp=6,empty=F,random=T)
  
  #empty
  
  empty_result500_bsp5 <- sim34(n=500,SIM=100,bsp=5,empty=T,random=F)
  empty_result1000_bsp5 <- sim34(n=1000,SIM=100,bsp=5,empty=T,random=F)
  empty_result5000_bsp5 <- sim34(n=5000,SIM=100,bsp=5,empty=T,random=F)
  
  empty_result500_bsp5_rand <- sim34(n=500,SIM=100,bsp=5,empty=T,random=T)
  empty_result1000_bsp5_rand <- sim34(n=1000,SIM=100,bsp=5,empty=T,random=T)
  empty_result5000_bsp5_rand <- sim34(n=5000,SIM=100,bsp=5,empty=T,random=T)
  
  empty_result500_bsp6 <- sim34(n=500,SIM=100,bsp=6,empty=T,random=F)
  empty_result1000_bsp6 <- sim34(n=1000,SIM=100,bsp=6,empty=T,random=F)
  empty_result5000_bsp6 <- sim34(n=5000,SIM=100,bsp=6,empty=T,random=F)
  
  empty_result500_bsp6_rand <- sim34(n=500,SIM=100,bsp=6,empty=T,random=T)
  empty_result1000_bsp6_rand <- sim34(n=1000,SIM=100,bsp=6,empty=T,random=T)
  empty_result5000_bsp6_rand <- sim34(n=5000,SIM=100,bsp=6,empty=T,random=T)
  
  
#  save(result500_bsp5,result1000_bsp5,result5000_bsp5,
#       result500_bsp5_rand,result1000_bsp5_rand,result5000_bsp5_rand,
#       result500_bsp6,result1000_bsp6,result5000_bsp6,
#       result500_bsp6_rand,result1000_bsp6_rand,result5000_bsp6_rand,
#       empty_result500_bsp5,empty_result1000_bsp5,empty_result5000_bsp5,
#       empty_result500_bsp5_rand,empty_result1000_bsp5_rand,empty_result5000_bsp5_rand,
#       empty_result500_bsp6,empty_result1000_bsp6,empty_result5000_bsp6,
#       empty_result500_bsp6_rand,empty_result1000_bsp6_rand,empty_result5000_bsp6_rand,
#       file="resultssim56.RData")  
  

load("resultssim56.RData") 



  library(ggplot2)
  
  #define middle_green
  green_rgb <- col2rgb("green")
  darkgreen_rgb <- col2rgb("darkgreen")
  middle_green_rgb <- (green_rgb + darkgreen_rgb) / 2
  middle_green <- rgb(middle_green_rgb[1], middle_green_rgb[2], middle_green_rgb[3], maxColorValue = 255)
  
  
  result<-result5000_bsp6
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
    
    if(empty){Z<-result$Z_Ddata_complete[,sim]
    if(bsp==5|bsp==6|bsp==7|bsp==8){
      W<-result$W_Ddata_complete[,sim]}}
    else{  Z<-result$Z_complete[,sim]
    if(bsp==5|bsp==6|bsp==7|bsp==8){
      W<-result$W_complete[,sim]}}
    
    
    #E[Y|do(X)]
    #example 3
    if(bsp==5){
      original<-function(x){
        mean(result$y0_complete)+mean(result$yx1_complete)*x+mean(result$yx2_complete)*x^2+mean(result$yx3_complete)*x^3+mean(result$yz1_complete)*(-0.3)+mean(result$yz2_complete)*(1.09)+mean(result$yz3_complete)*(-0.927)+mean(result$yw1_complete)*(2)+mean(result$yw2_complete)*(4.09)+mean(result$yw3_complete)*(8.54)}}
    
    #example 4
    if(bsp==6){
      original<-function(x){
        mean(result$y0_complete)+mean(result$yx1_complete)*x+mean(result$yx2_complete)*x^2+mean(result$yx3_complete)*x^3+mean(result$yz1_complete)*(0)+mean(result$yz2_complete)*(1)+mean(result$yz3_complete)*(0)+mean(result$yw1_complete)*(2)+mean(result$yw2_complete)*(4.09)+mean(result$yw3_complete)*(8.54)}}
    
    #E[Y|X]
    #example 3
    if(bsp==5){
      wrong<-function(x){
        mean(result$y0_complete)+mean(result$yx1_complete)*x+mean(result$yx2_complete)*x^2+mean(result$yx3_complete)*x^3+mean(result$yz1_complete)*(-0.3)+mean(result$yz2_complete)*(1.09)+mean(result$yz3_complete)*(-0.927)+mean(result$yw1_complete)*(2+0.3^2/(0.3^2+1)*(x-2))+mean(result$yw2_complete)*(4.09+(8.54-8-18)/1.09*(x-2))+mean(result$yw3_complete)*(8.54+(18.1843-17.08)/1.09*(x-2))}}
    
    #example 4
    if(bsp==6){
      wrong<-function(x){
        mean(result$y0_complete)+mean(result$yx1_complete)*x+mean(result$yx2_complete)*x^2+mean(result$yx3_complete)*x^3+mean(result$yz1_complete)*(0)+mean(result$yz2_complete)*(1)+mean(result$yz3_complete)*(0)+mean(result$yw1_complete)*(2+(0.3^2)/(0.3^2+1)*(x-2))+mean(result$yw2_complete)*(4.09+(8.54-8.18)/(0.3^2+1)*(x-2))+mean(result$yw3_complete)*(8.54+(18.1843-17.08)/(0.3^2+1)*(x-2))}}
    
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
    if(bsp==5 | bsp==6){
      est3_ridge<-function(x){
        return(betahat_first_ridge[1]+betahat_first_ridge[2]*x+betahat_first_ridge[3]*x^2+betahat_first_ridge[4]*x^3+
                 betahat_first_ridge[5]*mean(Z)+
                 betahat_first_ridge[6]*mean(Z^2)+
                 betahat_first_ridge[7]*mean(Z^3)+
                 betahat_first_ridge[8]*mean(W)+
                 betahat_first_ridge[9]*mean(W^2)+
                 betahat_first_ridge[10]*mean(W^3))}}
    
    #TSR OLS
    if(bsp==5 |bsp==6){
      est3<-function(x){
        return(betahat_first[1]+betahat_first[2]*x+betahat_first[3]*x^2+betahat_first[4]*x^3+
                 betahat_first[5]*mean(Z)+
                 betahat_first[6]*mean(Z^2)+
                 betahat_first[7]*mean(Z^3)+
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
  
  
  #example 3
  
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
    geom_line(aes(x = X_est, y = q_lower_est3, color = "TSR"),linetype="dashed", size = 1)+
    geom_line(aes(x = X_est, y = q_upper_est3, color = "TSR"),linetype="dashed", size = 1)+
    geom_line(aes(x = X_est, y = apply(X_est,1,wrong), color = "E[Y|X=x]"), linewidth = 1.3) +
    geom_line(aes(x = X_est, y = apply(X_est,1,original), color = "E[Y|do(X=x)]"), linewidth = 1.3) +
    theme_minimal() +
    coord_cartesian(xlim = c(-2.5, 7), ylim = c(-1, 10)) + 
    geom_boxplot(aes(x = X, y = -1), width = 0.5, fill = "grey", color = "black", alpha = 0.7) +
    geom_boxplot(aes(x = X_S, y = 0), width = 0.5, fill = "grey", color = "black", alpha = 0.7) +
    scale_color_manual(values = c("E[Y|X=x]" = "purple", "E[Y|do(X=x)]" = "black", "naive"="orange", "RR"="blue", "TSR"=middle_green), name = NULL) +
    scale_fill_manual(values = c("naive" = "orange", "RR" = "lightblue", "TSR" = "lightgreen"), name = NULL) +
    labs(x="x",y="")+
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
    geom_line(aes(x = X_est, y = apply(X_est,1,original), color = "E[Y|do(X=x)]"), linewidth = 1.3) +
    theme_minimal() +
    coord_cartesian(xlim = c(-2.5, 7), ylim = c(-1, 10)) + 
    geom_boxplot(aes(x = X, y = -1), width = 0.5, fill = "grey", color = "black", alpha = 0.7) +
    geom_boxplot(aes(x = X_S, y = 0), width = 0.5, fill = "grey", color = "black", alpha = 0.7) +
    scale_color_manual(values = c("E[Y|X=x]" = "purple", "E[Y|do(X=x)]" = "black", "naive"="orange", "RR (ridge)"="blue", "TSR (ridge)"=middle_green), name = NULL) +
    scale_fill_manual(values = c("naive" = "orange", "RR (ridge)" = "lightblue", "TSR (ridge)" = "lightgreen"), name = NULL) +
    labs(x="x",y="")+
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
  
  
  #example 4
  
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
    geom_line(aes(x = X_est, y = apply(X_est,1,original), color = "E[Y|do(X=x)]"), linewidth = 1.3) +
    theme_minimal() +
    coord_cartesian(xlim = c(-2.5, 7), ylim = c(1, 15)) + 
    geom_boxplot(aes(x = X, y = 1), width = 1, fill = "grey", color = "black", alpha = 0.7) +
    geom_boxplot(aes(x = X_S, y = 2.5), width = 1, fill = "grey", color = "black", alpha = 0.7) +
    scale_color_manual(values = c("E[Y|X=x]" = "purple", "E[Y|do(X=x)]" = "black", "naive"="orange", "RR"="blue", "TSR"=middle_green), name = NULL) +
    scale_fill_manual(values = c("naive" = "orange", "RR" = "lightblue", "TSR" = "lightgreen"), name = NULL) +
    labs(x="x",y="")+
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
    geom_line(aes(x = X_est, y = apply(X_est,1,original), color = "E[Y|do(X=x)]"), linewidth = 1.3) +
    theme_minimal() +
    coord_cartesian(xlim = c(-2.5, 7), ylim = c(1, 15)) + 
    geom_boxplot(aes(x = X, y = 1), width = 1, fill = "grey", color = "black", alpha = 0.7) +
    geom_boxplot(aes(x = X_S, y = 2.5), width = 1, fill = "grey", color = "black", alpha = 0.7) +
    scale_color_manual(values = c("E[Y|X=x]" = "purple", "E[Y|do(X=x)]" = "black", "naive"="orange", "RR (ridge)"="blue", "TSR (ridge)"=middle_green), name = NULL) +
    scale_fill_manual(values = c("naive" = "orange", "RR (ridge)" = "lightblue", "TSR (ridge)" = "lightgreen"), name = NULL) +
    labs(x="x",y="")+
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
  
  
  
  
  
load("resultssim56.RData")  

####################################
#######      bsp5     ########
####################################



#MSE mean
c(mean(result500_bsp5_rand$o_test_naive_S), mean(result500_bsp5_rand$o_test_est1_S), mean(result500_bsp5_rand$o_test_est1_ridge_S), mean(result500_bsp5_rand$o_test_est3_S), mean(result500_bsp5_rand$o_test_est3_ridge_S))
#[1] 0.8871010 0.7515690 0.5751506 0.3721857 0.2238777
c(mean(result1000_bsp5_rand$o_test_naive_S), mean(result1000_bsp5_rand$o_test_est1_S), mean(result1000_bsp5_rand$o_test_est1_ridge_S), mean(result1000_bsp5_rand$o_test_est3_S), mean(result1000_bsp5_rand$o_test_est3_ridge_S))
#[1] 0.8462348 0.4610474 0.5168566 0.1494094 0.1875617
c(mean(result5000_bsp5_rand$o_test_naive_S), mean(result5000_bsp5_rand$o_test_est1_S), mean(result5000_bsp5_rand$o_test_est1_ridge_S), mean(result5000_bsp5_rand$o_test_est3_S), mean(result5000_bsp5_rand$o_test_est3_ridge_S))
#[1] 0.74931447 0.29009292 0.38533744 0.02820336 0.09997391

c(mean(result500_bsp5_rand$o_test_naive), mean(result500_bsp5_rand$o_test_est1), mean(result500_bsp5_rand$o_test_est1_ridge), mean(result500_bsp5_rand$o_test_est3), mean(result500_bsp5_rand$o_test_est3_ridge))
#[1] 7.9140620 0.9686047 0.8724697 0.4915395 0.5127053
c(mean(result1000_bsp5_rand$o_test_naive), mean(result1000_bsp5_rand$o_test_est1), mean(result1000_bsp5_rand$o_test_est1_ridge), mean(result1000_bsp5_rand$o_test_est3), mean(result1000_bsp5_rand$o_test_est3_ridge))
#[1] 5.5252574 0.6645752 0.7665563 0.2199803 0.3765903
c(mean(result5000_bsp5_rand$o_test_naive), mean(result5000_bsp5_rand$o_test_est1), mean(result5000_bsp5_rand$o_test_est1_ridge), mean(result5000_bsp5_rand$o_test_est3), mean(result5000_bsp5_rand$o_test_est3_ridge))
#[1] 1.24550157 0.42358102 0.54419732 0.03754944 0.20492636


#MSE sd
c(sd(result500_bsp5_rand$o_test_naive_S), sd(result500_bsp5_rand$o_test_est1_S), sd(result500_bsp5_rand$o_test_est1_ridge_S), sd(result500_bsp5_rand$o_test_est3_S), sd(result500_bsp5_rand$o_test_est3_ridge_S))
#[1] 0.8414769 0.8179506 0.5507757 0.4502169 0.2546193
c(sd(result1000_bsp5_rand$o_test_naive_S), sd(result1000_bsp5_rand$o_test_est1_S), sd(result1000_bsp5_rand$o_test_est1_ridge_S), sd(result1000_bsp5_rand$o_test_est3_S), sd(result1000_bsp5_rand$o_test_est3_ridge_S))
#[1] 0.8954566 0.4378683 0.4022386 0.2003081 0.2104424
c(sd(result5000_bsp5_rand$o_test_naive_S), sd(result5000_bsp5_rand$o_test_est1_S), sd(result5000_bsp5_rand$o_test_est1_ridge_S), sd(result5000_bsp5_rand$o_test_est3_S), sd(result5000_bsp5_rand$o_test_est3_ridge_S))
#[1] 0.78646339 0.24723901 0.30353722 0.04890608 0.13675991

c(sd(result500_bsp5_rand$o_test_naive), sd(result500_bsp5_rand$o_test_est1), sd(result500_bsp5_rand$o_test_est1_ridge), sd(result500_bsp5_rand$o_test_est3), sd(result500_bsp5_rand$o_test_est3_ridge))
#[1] 13.0948361  0.8131654  0.7004721  0.4880217  0.5928714
c(sd(result1000_bsp5_rand$o_test_naive), sd(result1000_bsp5_rand$o_test_est1), sd(result1000_bsp5_rand$o_test_est1_ridge), sd(result1000_bsp5_rand$o_test_est3), sd(result1000_bsp5_rand$o_test_est3_ridge))
#[1] 10.1297568  0.4998519  0.5879058  0.2316439  0.3477438
c(sd(result5000_bsp5_rand$o_test_naive), sd(result5000_bsp5_rand$o_test_est1), sd(result5000_bsp5_rand$o_test_est1_ridge), sd(result5000_bsp5_rand$o_test_est3), sd(result5000_bsp5_rand$o_test_est3_ridge))
#[1] 1.21488162 0.33202446 0.33091559 0.05423269 0.21494682


#################
###   empty   ###
#################



#MSE mean
c(mean(empty_result500_bsp5_rand$o_test_naive_S), mean(empty_result500_bsp5_rand$o_test_est1_S), mean(empty_result500_bsp5_rand$o_test_est1_ridge_S), mean(empty_result500_bsp5_rand$o_test_est3_S), mean(empty_result500_bsp5_rand$o_test_est3_ridge_S))
#[1] 0.9160307 0.7923017 0.6188436 0.3961381 0.2311935
c(mean(empty_result1000_bsp5_rand$o_test_naive_S), mean(empty_result1000_bsp5_rand$o_test_est1_S), mean(empty_result1000_bsp5_rand$o_test_est1_ridge_S), mean(empty_result1000_bsp5_rand$o_test_est3_S), mean(empty_result1000_bsp5_rand$o_test_est3_ridge_S))
#[1] 0.8523185 0.4950933 0.5420266 0.1339101 0.1742297
c(mean(empty_result5000_bsp5_rand$o_test_naive_S), mean(empty_result5000_bsp5_rand$o_test_est1_S), mean(empty_result5000_bsp5_rand$o_test_est1_ridge_S), mean(empty_result5000_bsp5_rand$o_test_est3_S), mean(empty_result5000_bsp5_rand$o_test_est3_ridge_S))
#[1] 0.74732869 0.29318234 0.38442243 0.02921946 0.09794899

c(mean(empty_result500_bsp5_rand$o_test_naive), mean(empty_result500_bsp5_rand$o_test_est1), mean(empty_result500_bsp5_rand$o_test_est1_ridge), mean(empty_result500_bsp5_rand$o_test_est3), mean(empty_result500_bsp5_rand$o_test_est3_ridge))
#[1] 7.5757045 1.0044049 0.9081302 0.5226274 0.5094343
c(mean(empty_result1000_bsp5_rand$o_test_naive), mean(empty_result1000_bsp5_rand$o_test_est1), mean(empty_result1000_bsp5_rand$o_test_est1_ridge), mean(empty_result1000_bsp5_rand$o_test_est3), mean(empty_result1000_bsp5_rand$o_test_est3_ridge))
#[1] 5.6665632 0.6710650 0.7582574 0.2063315 0.3665399
c(mean(empty_result5000_bsp5_rand$o_test_naive), mean(empty_result5000_bsp5_rand$o_test_est1), mean(empty_result5000_bsp5_rand$o_test_est1_ridge), mean(empty_result5000_bsp5_rand$o_test_est3), mean(empty_result5000_bsp5_rand$o_test_est3_ridge))
#[1] 1.22858788 0.42159402 0.54303343 0.03833992 0.20429323

#MSE sd
c(sd(empty_result500_bsp5_rand$o_test_naive_S), sd(empty_result500_bsp5_rand$o_test_est1_S), sd(empty_result500_bsp5_rand$o_test_est1_ridge_S), sd(empty_result500_bsp5_rand$o_test_est3_S), sd(empty_result500_bsp5_rand$o_test_est3_ridge_S))
#[1] 0.8401317 0.8289573 0.5873773 0.4972770 0.2542661
c(sd(empty_result1000_bsp5_rand$o_test_naive_S), sd(empty_result1000_bsp5_rand$o_test_est1_S), sd(empty_result1000_bsp5_rand$o_test_est1_ridge_S), sd(empty_result1000_bsp5_rand$o_test_est3_S), sd(empty_result1000_bsp5_rand$o_test_est3_ridge_S))
#[1] 0.8854529 0.5008554 0.4537136 0.1836580 0.1983837
c(sd(empty_result5000_bsp5_rand$o_test_naive_S), sd(empty_result5000_bsp5_rand$o_test_est1_S), sd(empty_result5000_bsp5_rand$o_test_est1_ridge_S), sd(empty_result5000_bsp5_rand$o_test_est3_S), sd(empty_result5000_bsp5_rand$o_test_est3_ridge_S))
#[1] 0.78814405 0.24506807 0.28652436 0.04561667 0.13814917

c(sd(empty_result500_bsp5_rand$o_test_naive), sd(empty_result500_bsp5_rand$o_test_est1), sd(empty_result500_bsp5_rand$o_test_est1_ridge), sd(empty_result500_bsp5_rand$o_test_est3), sd(empty_result500_bsp5_rand$o_test_est3_ridge))
#[1] 10.7030146  0.8550446  0.7374576  0.5204136  0.5570426
c(sd(empty_result1000_bsp5_rand$o_test_naive), sd(empty_result1000_bsp5_rand$o_test_est1), sd(empty_result1000_bsp5_rand$o_test_est1_ridge), sd(empty_result1000_bsp5_rand$o_test_est3), sd(empty_result1000_bsp5_rand$o_test_est3_ridge))
#[1] 11.5292429  0.5074177  0.5696566  0.2132908  0.3690750
c(sd(empty_result5000_bsp5_rand$o_test_naive), sd(empty_result5000_bsp5_rand$o_test_est1), sd(empty_result5000_bsp5_rand$o_test_est1_ridge), sd(empty_result5000_bsp5_rand$o_test_est3), sd(empty_result5000_bsp5_rand$o_test_est3_ridge))
#[1] 1.18760759 0.32938699 0.32458594 0.05160971 0.20815128





####################################
#######      bsp6     ########
####################################



#MSE mean
c(mean(result500_bsp6_rand$o_test_naive_S), mean(result500_bsp6_rand$o_test_est1_S), mean(result500_bsp6_rand$o_test_est1_ridge_S), mean(result500_bsp6_rand$o_test_est3_S), mean(result500_bsp6_rand$o_test_est3_ridge_S))
#[1] 13.64883 11.89703 11.75123 13.01403 12.75425
c(mean(result1000_bsp6_rand$o_test_naive_S), mean(result1000_bsp6_rand$o_test_est1_S), mean(result1000_bsp6_rand$o_test_est1_ridge_S), mean(result1000_bsp6_rand$o_test_est3_S), mean(result1000_bsp6_rand$o_test_est3_ridge_S))
#[1] 14.63763 13.30937 13.29217 14.44939 14.42298
c(mean(result5000_bsp6_rand$o_test_naive_S), mean(result5000_bsp6_rand$o_test_est1_S), mean(result5000_bsp6_rand$o_test_est1_ridge_S), mean(result5000_bsp6_rand$o_test_est3_S), mean(result5000_bsp6_rand$o_test_est3_ridge_S))
#[1] 11.74511 11.69364 11.69098 13.10261 13.09618

c(mean(result500_bsp6_rand$o_test_naive), mean(result500_bsp6_rand$o_test_est1), mean(result500_bsp6_rand$o_test_est1_ridge), mean(result500_bsp6_rand$o_test_est3), mean(result500_bsp6_rand$o_test_est3_ridge))
#[1] 17.73657 14.94052 13.87487 14.76356 13.61439
c(mean(result1000_bsp6_rand$o_test_naive), mean(result1000_bsp6_rand$o_test_est1), mean(result1000_bsp6_rand$o_test_est1_ridge), mean(result1000_bsp6_rand$o_test_est3), mean(result1000_bsp6_rand$o_test_est3_ridge))
#[1] 17.16311 14.84310 15.07657 14.57342 14.75627
c(mean(result5000_bsp6_rand$o_test_naive), mean(result5000_bsp6_rand$o_test_est1), mean(result5000_bsp6_rand$o_test_est1_ridge), mean(result5000_bsp6_rand$o_test_est3), mean(result5000_bsp6_rand$o_test_est3_ridge))
#[1] 13.31393 13.44088 13.45028 13.19189 13.20349


#MSE sd
c(sd(result500_bsp6_rand$o_test_naive_S), sd(result500_bsp6_rand$o_test_est1_S), sd(result500_bsp6_rand$o_test_est1_ridge_S), sd(result500_bsp6_rand$o_test_est3_S), sd(result500_bsp6_rand$o_test_est3_ridge_S))
#[1] 16.51736 14.15529 14.00174 15.53848 15.22964
c(sd(result1000_bsp6_rand$o_test_naive_S), sd(result1000_bsp6_rand$o_test_est1_S), sd(result1000_bsp6_rand$o_test_est1_ridge_S), sd(result1000_bsp6_rand$o_test_est3_S), sd(result1000_bsp6_rand$o_test_est3_ridge_S))
#[1] 16.62597 15.08471 15.06607 16.88409 16.84145
c(sd(result5000_bsp6_rand$o_test_naive_S), sd(result5000_bsp6_rand$o_test_est1_S), sd(result5000_bsp6_rand$o_test_est1_ridge_S), sd(result5000_bsp6_rand$o_test_est3_S), sd(result5000_bsp6_rand$o_test_est3_ridge_S))
#[1] 14.17333 14.03226 14.02567 15.94039 15.92759

c(sd(result500_bsp6_rand$o_test_naive), sd(result500_bsp6_rand$o_test_est1), sd(result500_bsp6_rand$o_test_est1_ridge), sd(result500_bsp6_rand$o_test_est3), sd(result500_bsp6_rand$o_test_est3_ridge))
#[1] 19.86169 18.83754 15.49859 19.40734 15.34375
c(sd(result1000_bsp6_rand$o_test_naive), sd(result1000_bsp6_rand$o_test_est1), sd(result1000_bsp6_rand$o_test_est1_ridge), sd(result1000_bsp6_rand$o_test_est3), sd(result1000_bsp6_rand$o_test_est3_ridge))
#[1] 18.86864 16.86368 16.93955 16.78994 16.85991
c(sd(result5000_bsp6_rand$o_test_naive), sd(result5000_bsp6_rand$o_test_est1), sd(result5000_bsp6_rand$o_test_est1_ridge), sd(result5000_bsp6_rand$o_test_est3), sd(result5000_bsp6_rand$o_test_est3_ridge))
#[1] 15.99285 16.16629 16.15244 16.04764 16.03939


#################
###   empty   ###
#################



#MSE mean
c(mean(empty_result500_bsp6_rand$o_test_naive_S), mean(empty_result500_bsp6_rand$o_test_est1_S), mean(empty_result500_bsp6_rand$o_test_est1_ridge_S), mean(empty_result500_bsp6_rand$o_test_est3_S), mean(empty_result500_bsp6_rand$o_test_est3_ridge_S))
#[1] 1.9755558 0.7966609 0.5768445 0.4666676 0.2626173
c(mean(empty_result1000_bsp6_rand$o_test_naive_S), mean(empty_result1000_bsp6_rand$o_test_est1_S), mean(empty_result1000_bsp6_rand$o_test_est1_ridge_S), mean(empty_result1000_bsp6_rand$o_test_est3_S), mean(empty_result1000_bsp6_rand$o_test_est3_ridge_S))
#[1] 1.2398684 0.4189283 0.4093819 0.1063182 0.1101837
c(mean(empty_result5000_bsp6_rand$o_test_naive_S), mean(empty_result5000_bsp6_rand$o_test_est1_S), mean(empty_result5000_bsp6_rand$o_test_est1_ridge_S), mean(empty_result5000_bsp6_rand$o_test_est3_S), mean(empty_result5000_bsp6_rand$o_test_est3_ridge_S))
#[1] 0.58965573 0.22439769 0.22530722 0.01377512 0.01541318

c(mean(empty_result500_bsp6_rand$o_test_naive), mean(empty_result500_bsp6_rand$o_test_est1), mean(empty_result500_bsp6_rand$o_test_est1_ridge), mean(empty_result500_bsp6_rand$o_test_est3), mean(empty_result500_bsp6_rand$o_test_est3_ridge))
#[1] 6.226116 2.325862 1.238454 2.070496 1.034836
c(mean(empty_result1000_bsp6_rand$o_test_naive), mean(empty_result1000_bsp6_rand$o_test_est1), mean(empty_result1000_bsp6_rand$o_test_est1_ridge), mean(empty_result1000_bsp6_rand$o_test_est3), mean(empty_result1000_bsp6_rand$o_test_est3_ridge))
#[1] 2.6086369 0.7293561 0.7445886 0.4462002 0.4834856
c(mean(empty_result5000_bsp6_rand$o_test_naive), mean(empty_result5000_bsp6_rand$o_test_est1), mean(empty_result5000_bsp6_rand$o_test_est1_ridge), mean(empty_result5000_bsp6_rand$o_test_est3), mean(empty_result5000_bsp6_rand$o_test_est3_ridge))
#[1] 0.52434315 0.19285638 0.20465903 0.03988480 0.05359514

#MSE sd
c(sd(empty_result500_bsp6_rand$o_test_naive_S), sd(empty_result500_bsp6_rand$o_test_est1_S), sd(empty_result500_bsp6_rand$o_test_est1_ridge_S), sd(empty_result500_bsp6_rand$o_test_est3_S), sd(empty_result500_bsp6_rand$o_test_est3_ridge_S))
#[1] 4.6796716 1.7548839 0.5684180 1.5199896 0.3386209
c(sd(empty_result1000_bsp6_rand$o_test_naive_S), sd(empty_result1000_bsp6_rand$o_test_est1_S), sd(empty_result1000_bsp6_rand$o_test_est1_ridge_S), sd(empty_result1000_bsp6_rand$o_test_est3_S), sd(empty_result1000_bsp6_rand$o_test_est3_ridge_S))
#[1] 1.2870403 0.4634971 0.4243299 0.1935062 0.1660521
c(sd(empty_result5000_bsp6_rand$o_test_naive_S), sd(empty_result5000_bsp6_rand$o_test_est1_S), sd(empty_result5000_bsp6_rand$o_test_est1_ridge_S), sd(empty_result5000_bsp6_rand$o_test_est3_S), sd(empty_result5000_bsp6_rand$o_test_est3_ridge_S))
#[1] 0.57792187 0.25802348 0.25673640 0.01198269 0.01480226

c(sd(empty_result500_bsp6_rand$o_test_naive), sd(empty_result500_bsp6_rand$o_test_est1), sd(empty_result500_bsp6_rand$o_test_est1_ridge), sd(empty_result500_bsp6_rand$o_test_est3), sd(empty_result500_bsp6_rand$o_test_est3_ridge))
#[1] 27.289799 11.343734  1.691686 10.688336  1.552364
c(sd(empty_result1000_bsp6_rand$o_test_naive), sd(empty_result1000_bsp6_rand$o_test_est1), sd(empty_result1000_bsp6_rand$o_test_est1_ridge), sd(empty_result1000_bsp6_rand$o_test_est3), sd(empty_result1000_bsp6_rand$o_test_est3_ridge))
#[1] 4.6950369 1.2022703 1.0558573 0.8643641 0.7485699
c(sd(empty_result5000_bsp6_rand$o_test_naive), sd(empty_result5000_bsp6_rand$o_test_est1), sd(empty_result5000_bsp6_rand$o_test_est1_ridge), sd(empty_result5000_bsp6_rand$o_test_est3), sd(empty_result5000_bsp6_rand$o_test_est3_ridge))
#[1] 0.47048199 0.19866375 0.19976567 0.04746704 0.06614377
































































