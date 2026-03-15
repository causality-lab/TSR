#Section 5.2

#examples 1 & 2

#load packages
library(glmnet)
library(ggplot2)

#n = sample size
#SIM = number of simulation runs
#bsp = example number
#      note that bsp i corresponds to example i-2 (i.e. bsp3 corresponds to example 1)
#empty = indicates if ((S cap D)=emptyset)
#random = indicates uniformly distributed coefficients

sim12<-function(n,SIM,bsp,empty,random)
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
  
  if(bsp==3 | bsp==4){
    if(empty==F){
      X_complete <- matrix(ncol=SIM,nrow=n)
      Z_complete <- matrix(ncol=SIM,nrow=n)
      S_complete <- matrix(ncol=SIM,nrow=n)
      Y_complete <- matrix(ncol=SIM,nrow=n)
    }else{
      X_Sdata_complete <- matrix(ncol=SIM,nrow=n)
      Z_Sdata_complete <- matrix(ncol=SIM,nrow=n)
      S_Sdata_complete <- matrix(ncol=SIM,nrow=n)
      Y_Sdata_complete <- matrix(ncol=SIM,nrow=n)
      
      X_Ddata_complete <- matrix(ncol=SIM,nrow=n)
      Z_Ddata_complete <- matrix(ncol=SIM,nrow=n)}}
  
  if(bsp==3 | bsp==4){
    betahat_first_complete <- matrix(ncol=SIM,nrow=7)
    betahat_first_ridge_complete <- matrix(ncol=SIM,nrow=7)
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
  }
  
  for(sim in 1:SIM)
  {
    print(sim)
    
    #reproducibility
    set.seed(sim)
    
    #data generating processes
    #example 1
    if(bsp==3){
      if(empty==F){
        Z<-rnorm(n,mean=-2)
        X<-rnorm(n)+2*Z
        S<-X+Z < -5
        if(random==F)
        {y0<-0
        yx1<-0
        yx2<-0.2
        yx3<-0.1
        yz1<-2
        yz2<-0
        yz3<-0.1
        }else{
          y0<-runif(1)-0.5
          yx1<-runif(1)-0.5
          yx2<-runif(1)-0.3
          yx3<-runif(1)-0.4
          yz1<-runif(1)+1.5
          yz2<-runif(1)-0.5
          yz3<-runif(1)-0.4
        }
        y0_complete[sim]<-y0
        yx1_complete[sim]<-yx1
        yx2_complete[sim]<-yx2
        yx3_complete[sim]<-yx3
        yz1_complete[sim]<-yz1
        yz2_complete[sim]<-yz2
        yz3_complete[sim]<-yz3
        Y<-y0+yx1*X+yx2*X^2+yx3*X^3+yz1*Z+yz2*Z^2+yz3*Z^3+rnorm(n)  
        Z_S<-Z[S==1]
        X_S<-X[S==1]
        Y_S<-Y[S==1]
      }else{
        Z_Sdata<-rnorm(n,mean=-2)
        X_Sdata<-rnorm(n)+2*Z_Sdata
        S_Sdata<-X_Sdata+Z_Sdata < -5
        if(random==F)
        {y0<-0
        yx1<-0
        yx2<-0.2
        yx3<-0.1
        yz1<-2
        yz2<-0
        yz3<-0.1
        }else{
          y0<-runif(1)-0.5
          yx1<-runif(1)-0.5
          yx2<-runif(1)-0.3
          yx3<-runif(1)-0.4
          yz1<-runif(1)+1.5
          yz2<-runif(1)-0.5
          yz3<-runif(1)-0.4
        }
        y0_complete[sim]<-y0
        yx1_complete[sim]<-yx1
        yx2_complete[sim]<-yx2
        yx3_complete[sim]<-yx3
        yz1_complete[sim]<-yz1
        yz2_complete[sim]<-yz2
        yz3_complete[sim]<-yz3
        Y_Sdata<-y0+yx1*X_Sdata+yx2*X_Sdata^2+yx3*X_Sdata^3+yz1*Z_Sdata+yz2*Z_Sdata^2+yz3*Z_Sdata^3+rnorm(n) 
        
        Z_S<-Z_Sdata[S_Sdata==1]
        X_S<-X_Sdata[S_Sdata==1]
        Y_S<-Y_Sdata[S_Sdata==1]
        
        Z_Ddata<-rnorm(n,mean=-2)
        X_Ddata<-rnorm(n)+2*Z_Ddata
        
        X<-X_Ddata
        Z<-Z_Ddata}}
    #example 2
    if(bsp==4){
      if(empty==F){
        Z<-rnorm(n,mean=-1,sd=2)
        X<-rnorm(n)+Z
        S<-rbinom(p=1/((1+exp(-X))*(1+exp(Z))),size=1,n)
        if(random==F){
          y0<-0
          yx1<-1
          yx2<-0
          yx3<-0
          yz1<-5
          yz2<-0
          yz3<-0
        }else{
          y0<-runif(1)-0.5
          yx1<-runif(1)+0.5
          yx2<-runif(1)-0.5
          yx3<-runif(1)-0.5
          yz1<-runif(1)+4.5
          yz2<-runif(1)-0.5
          yz3<-runif(1)-0.5
        }
        y0_complete[sim]<-y0
        yx1_complete[sim]<-yx1
        yx2_complete[sim]<-yx2
        yx3_complete[sim]<-yx3
        yz1_complete[sim]<-yz1
        yz2_complete[sim]<-yz2
        yz3_complete[sim]<-yz3
        Y<-y0+yx1*X+yx2*X^2+yx3*X^3+yz1*Z+yz2*Z^2+yz3*Z^3+rnorm(n) 
        
        Z_S<-Z[S==1]
        X_S<-X[S==1]
        Y_S<-Y[S==1]
      }else{
        Z_Sdata<-rnorm(n,mean=-1,sd=2)
        X_Sdata<--rnorm(n)+Z_Sdata
        S_Sdata<-rbinom(p=1/((1+exp(-X_Sdata))*(1+exp(Z_Sdata))),size=1,n)
        if(random==F){
          y0<-0
          yx1<-1
          yx2<-0
          yx3<-0
          yz1<-5
          yz2<-0
          yz3<-0
        }else{
          y0<-runif(1)-0.5
          yx1<-runif(1)+0.5
          yx2<-runif(1)-0.5
          yx3<-runif(1)-0.5
          yz1<-runif(1)+4.5
          yz2<-runif(1)-0.5
          yz3<-runif(1)-0.5
        }
        y0_complete[sim]<-y0
        yx1_complete[sim]<-yx1
        yx2_complete[sim]<-yx2
        yx3_complete[sim]<-yx3
        yz1_complete[sim]<-yz1
        yz2_complete[sim]<-yz2
        yz3_complete[sim]<-yz3
        Y_Sdata<-y0+yx1*X_Sdata+yx2*X_Sdata^2+yx3*X_Sdata^3+yz1*Z_Sdata+yz2*Z_Sdata^2+yz3*Z_Sdata^3+rnorm(n)
        
        Z_S<-Z_Sdata[S_Sdata==1]
        X_S<-X_Sdata[S_Sdata==1]
        Y_S<-Y_Sdata[S_Sdata==1]
        
        Z_Ddata<-rnorm(n,mean=-1,sd=2)
        X_Ddata<-rnorm(n)+Z_Ddata
        
        X<-X_Ddata
        Z<-Z_Ddata}}
   
      #save realisations from all simulation runs
      if(bsp==3 | bsp==4){
        if(empty==F){
          X_complete[,sim]<-X
          Z_complete[,sim]<-Z
          S_complete[,sim]<-S
          Y_complete[,sim]<-Y
        }else{
          X_Sdata_complete[,sim] <- X_Sdata
          Z_Sdata_complete[,sim] <- Z_Sdata
          S_Sdata_complete[,sim] <- S_Sdata
          Y_Sdata_complete[,sim] <- Y_Sdata
          
          X_Ddata_complete[,sim] <- X_Ddata
          Z_Ddata_complete[,sim] <- Z_Ddata}}
      
      #first step ridge 
      if(bsp==3 | bsp==4){
        lambda_seq <- 10^seq(2, -2, by = -.1)
        ridge_cv <- cv.glmnet(cbind(X_S,I(X_S^2),I(X_S^3),Z_S,I(Z_S^2),I(Z_S^3)), Y_S, alpha = 0, lambda = lambda_seq)
        best_lambda <- ridge_cv$lambda.min
        best_ridge <- glmnet(cbind(X_S,I(X_S^2),I(X_S^3),Z_S,I(Z_S^2),I(Z_S^3)), Y_S, alpha = 0, lambda = best_lambda)
        betahat_first_ridge<-as.vector(coef(best_ridge))}
      
      #first step OLS
      if(bsp==3 | bsp==4){
        lm_fit <- lm(Y_S~X_S+I(X_S^2)+I(X_S^3)+Z_S+I(Z_S^2)+I(Z_S^3))
        betahat_first<-coef(lm_fit)}
      
      #calculate target variable for second step of RR
      if(bsp==3 | bsp==4){
        Y_second_ridge<-(cbind(rep(1,n),X,X^2,X^3,Z,Z^2,Z^3))%*%as.vector(betahat_first_ridge)
        Y_second<-(cbind(rep(1,n),X,X^2,X^3,Z,Z^2,Z^3))%*%as.vector(betahat_first)}
      
      #final RR estimates
      beta_end1_ridge<-lm(Y_second_ridge~1+X+I(X^2)+I(X^3))$coef
      beta_end1<-lm(Y_second~1+X+I(X^2)+I(X^3))$coef
      
      #naive estimation only based on S=1
      naive_beta<-lm(Y_S~1+X_S+I(X_S^2)+I(X_S^3))$coef
      
      if(random==T){
      #E[Y|do(X)]
      #example 1
      if(bsp==3){
        original<-function(x){
          y0+yx1*x+yx2*x^2+yx3*x^3+yz1*(-2)+yz2*(5)+yz3*(-14)}}
      
      #example 2
      if(bsp==4){
        original<-function(x){
          y0+yx1*x+yx2*x^2+yx3*x^3+yz1*(-1)+yz2*(5)+yz3*(-13)}}
   
      #E[Y|X]
      #example 1
      if(bsp==3){
        wrong<-function(x){y0+yx1*x+yx2*x^2+yx3*x^3+yz1*(-2/5+2/5*x)+yz2*(5-8/5*(x+4))+yz3*(10+6*x)}}
      
      #example 2
      if(bsp==4){
        wrong<-function(x){y0+yx1*x+yx2*x^2+yx3*x^3+yz1*(-1/5+4/5*x)+yz2*(5-8/5*(x+1))+yz3*(-13+60/5+(x+1))}}
      
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
      if(bsp==3 | bsp==4){
        est3_ridge<-function(x){
          return(betahat_first_ridge[1]+betahat_first_ridge[2]*x+betahat_first_ridge[3]*x^2+betahat_first_ridge[4]*x^3+
                   betahat_first_ridge[5]*mean(Z)+
                   betahat_first_ridge[6]*mean(Z^2)+
                   betahat_first_ridge[7]*mean(Z^3))}}
      
      #TSR OLS
      if(bsp==3 | bsp==4){
        est3<-function(x){
          return(betahat_first[1]+betahat_first[2]*x+betahat_first[3]*x^2+betahat_first[4]*x^3+
                   betahat_first[5]*mean(Z)+
                   betahat_first[6]*mean(Z^2)+
                   betahat_first[7]*mean(Z^3))}}
      }
      #save coefficient vectors from all simulation runs
      betahat_first_complete[,sim] <- betahat_first
      betahat_first_ridge_complete[,sim] <- betahat_first_ridge
      beta_end1_complete[,sim] <- beta_end1
      beta_end1_ridge_complete[,sim] <- beta_end1_ridge
      naive_beta_complete[,sim] <- naive_beta
      
      if(random==T){
      #test data    
      if(bsp==3){
        Z_test<-rnorm(n,mean=-2)
        X_test<-as.matrix(rnorm(n)+2*Z_test)
        S_test<-X_test+Z_test < -5   
        Y_test<-y0+yx1*X_test+yx2*X_test^2+yx3*X_test^3+yz1*Z_test+yz2*Z_test^2+yz3*Z_test^3  
        X_test_S<-as.matrix(X_test[S_test==1,])
        Y_test_S<-Y_test[S_test==1]}
      if(bsp==4){
        Z_test<-rnorm(n,mean=-1,sd=2)
        X_test<-as.matrix(rnorm(n)+Z_test)
        S_test<-rbinom(p=1/((1+exp(-X_test))*(1+exp(Z_test))),size=1,n)
        Y_test<-y0+yx1*X_test+yx2*X_test^2+yx3*X_test^3+yz1*Z_test+yz2*Z_test^2+yz3*Z_test^3  
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
    if(bsp==3 | bsp==4){
      
      if(random==F){
      if(empty==F){
        return(list(SIM = SIM, n = n , bsp = bsp, empty=empty,
                    X_complete = X_complete, Z_complete = Z_complete, S_complete = S_complete, Y_complete = Y_complete,
                    betahat_first_complete = betahat_first_complete, betahat_first_ridge_complete = betahat_first_ridge_complete,
                    beta_end1_complete = beta_end1_complete, beta_end1_ridge_complete = beta_end1_ridge_complete,
                    naive_beta_complete = naive_beta_complete,
                    y0_complete=y0_complete,yx1_complete=yx1_complete,yx2_complete=yx2_complete,yx3_complete=yx3_complete,yz1_complete=yz1_complete,yz2_complete=yz2_complete,yz3_complete=yz3_complete))
      }else{
        return(list(SIM = SIM, n = n , bsp = bsp, 
                    X_Sdata_complete = X_Sdata_complete, Z_Sdata_complete = Z_Sdata_complete, S_Sdata_complete = S_Sdata_complete, Y_Sdata_complete = Y_Sdata_complete,
                    X_Ddata_complete = X_Ddata_complete, Z_Ddata_complete = Z_Ddata_complete,
                    betahat_first_complete = betahat_first_complete, betahat_first_ridge_complete = betahat_first_ridge_complete,
                    beta_end1_complete = beta_end1_complete, beta_end1_ridge_complete = beta_end1_ridge_complete,
                    naive_beta_complete = naive_beta_complete,
                    y0_complete=y0_complete,yx1_complete=yx1_complete,yx2_complete=yx2_complete,yx3_complete=yx3_complete,yz1_complete=yz1_complete,yz2_complete=yz2_complete,yz3_complete=yz3_complete))
      }}else{ 
        if(empty==F){
        return(list(SIM = SIM, n = n , bsp = bsp, empty=empty,
                    X_complete = X_complete, Z_complete = Z_complete, S_complete = S_complete, Y_complete = Y_complete,
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
                    y0_complete=y0_complete,yx1_complete=yx1_complete,yx2_complete=yx2_complete,yx3_complete=yx3_complete,yz1_complete=yz1_complete,yz2_complete=yz2_complete,yz3_complete=yz3_complete))
      }else{
        return(list(SIM = SIM, n = n , bsp = bsp, 
                    X_Sdata_complete = X_Sdata_complete, Z_Sdata_complete = Z_Sdata_complete, S_Sdata_complete = S_Sdata_complete, Y_Sdata_complete = Y_Sdata_complete,
                    X_Ddata_complete = X_Ddata_complete, Z_Ddata_complete = Z_Ddata_complete,
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
                    y0_complete=y0_complete,yx1_complete=yx1_complete,yx2_complete=yx2_complete,yx3_complete=yx3_complete,yz1_complete=yz1_complete,yz2_complete=yz2_complete,yz3_complete=yz3_complete))
      }
        
      }
        
        }
  }
  
  
  
  
  #####################################################################################################
  #####################################################################################################
  #####################################################################################################
  
  #subset
  
  result500_bsp3 <- sim12(n=500,SIM=100,bsp=3,empty=F,random=F)
  result1000_bsp3 <- sim12(n=1000,SIM=100,bsp=3,empty=F,random=F)
  result5000_bsp3 <- sim12(n=5000,SIM=100,bsp=3,empty=F,random=F)
  
  result500_bsp3_rand <- sim12(n=500,SIM=100,bsp=3,empty=F,random=T)
  result1000_bsp3_rand <- sim12(n=1000,SIM=100,bsp=3,empty=F,random=T)
  result5000_bsp3_rand <- sim12(n=5000,SIM=100,bsp=3,empty=F,random=T)
  
  result500_bsp4 <- sim12(n=500,SIM=100,bsp=4,empty=F,random=F)
  result1000_bsp4 <- sim12(n=1000,SIM=100,bsp=4,empty=F,random=F)
  result5000_bsp4 <- sim12(n=5000,SIM=100,bsp=4,empty=F,random=F)
  
  result500_bsp4_rand <- sim12(n=500,SIM=100,bsp=4,empty=F,random=T)
  result1000_bsp4_rand <- sim12(n=1000,SIM=100,bsp=4,empty=F,random=T)
  result5000_bsp4_rand <- sim12(n=5000,SIM=100,bsp=4,empty=F,random=T)
  
  #empty
  
  empty_result500_bsp3 <- sim12(n=500,SIM=100,bsp=3,empty=T,random=F)
  empty_result1000_bsp3 <- sim12(n=1000,SIM=100,bsp=3,empty=T,random=F)
  empty_result5000_bsp3 <- sim12(n=5000,SIM=100,bsp=3,empty=T,random=F)
  
  empty_result500_bsp3_rand <- sim12(n=500,SIM=100,bsp=3,empty=T,random=T)
  empty_result1000_bsp3_rand <- sim12(n=1000,SIM=100,bsp=3,empty=T,random=T)
  empty_result5000_bsp3_rand <- sim12(n=5000,SIM=100,bsp=3,empty=T,random=T)
  
  empty_result500_bsp4 <- sim12(n=500,SIM=100,bsp=4,empty=T,random=F)
  empty_result1000_bsp4 <- sim12(n=1000,SIM=100,bsp=4,empty=T,random=F)
  empty_result5000_bsp4 <- sim12(n=5000,SIM=100,bsp=4,empty=T,random=F)
  
  empty_result500_bsp4_rand <- sim12(n=500,SIM=100,bsp=4,empty=T,random=T)
  empty_result1000_bsp4_rand <- sim12(n=1000,SIM=100,bsp=4,empty=T,random=T)
  empty_result5000_bsp4_rand <- sim12(n=5000,SIM=100,bsp=4,empty=T,random=T)

#  save(result500_bsp3,result1000_bsp3,result5000_bsp3,
#       result500_bsp3_rand,result1000_bsp3_rand,result5000_bsp3_rand,
#       result500_bsp4,result1000_bsp4,result5000_bsp4,
#       result500_bsp4_rand,result1000_bsp4_rand,result5000_bsp4_rand,
#       empty_result500_bsp3,empty_result1000_bsp3,empty_result5000_bsp3,
#       empty_result500_bsp3_rand,empty_result1000_bsp3_rand,empty_result5000_bsp3_rand,
#       empty_result500_bsp4,empty_result1000_bsp4,empty_result5000_bsp4,
#       empty_result500_bsp4_rand,empty_result1000_bsp4_rand,empty_result5000_bsp4_rand,
#       file="resultssim34.RData")  
  
load("resultssim34.RData")  
  
  

  
  library(ggplot2)
  
  #define middle_green
  green_rgb <- col2rgb("green")
  darkgreen_rgb <- col2rgb("darkgreen")
  middle_green_rgb <- (green_rgb + darkgreen_rgb) / 2
  middle_green <- rgb(middle_green_rgb[1], middle_green_rgb[2], middle_green_rgb[3], maxColorValue = 255)
  
  
  result<-empty_result5000_bsp4
  empty<-TRUE
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
    }else{ Z<-result$Z_complete[,sim]}
    
    #E[Y|do(X)]
    #example 1
    if(bsp==3){
      original<-function(x){
        mean(result$y0_complete)+mean(result$yx1_complete)*x+mean(result$yx2_complete)*x^2+mean(result$yx3_complete)*x^3+mean(result$yz1_complete)*(-2)+mean(result$yz2_complete)*(5)+mean(result$yz3_complete)*(-14)}}
    
    #example 2
    if(bsp==4){
      original<-function(x){
        mean(result$y0_complete)+mean(result$yx1_complete)*x+mean(result$yx2_complete)*x^2+mean(result$yx3_complete)*x^3+mean(result$yz1_complete)*(-1)+mean(result$yz2_complete)*5+mean(result$yz3_complete)*(-13)}}
    
    #E[Y|X]
    #example 1
    if(bsp==3){
      wrong<-function(x){mean(result$y0_complete)+mean(result$yx1_complete)*x+mean(result$yx2_complete)*x^2+mean(result$yx3_complete)*x^3+mean(result$yz1_complete)*(-2/5+2/5*x)+mean(result$yz2_complete)*(5-8/5*(x+4))+mean(result$yz3_complete)*(10+6*x)}}
    
    #example 2
    if(bsp==4){
      wrong<-function(x){mean(result$y0_complete)+mean(result$yx1_complete)*x+mean(result$yx2_complete)*x^2+mean(result$yx3_complete)*x^3+mean(result$yz1_complete)*(-1/5+4/5*x)+mean(result$yz2_complete)*(5-8/5*(x+1))+mean(result$yz3_complete)*(-13+60/5+(x+1))}}
    
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
    if(bsp==3 | bsp==4){
      est3_ridge<-function(x){
        return(betahat_first_ridge[1]+betahat_first_ridge[2]*x+betahat_first_ridge[3]*x^2+betahat_first_ridge[4]*x^3+
                 betahat_first_ridge[5]*mean(Z)+
                 betahat_first_ridge[6]*mean(Z^2)+
                 betahat_first_ridge[7]*mean(Z^3))}}
    
    #TSR OLS
    if(bsp==3 | bsp==4){
      est3<-function(x){
        return(betahat_first[1]+betahat_first[2]*x+betahat_first[3]*x^2+betahat_first[4]*x^3+
                 betahat_first[5]*mean(Z)+
                 betahat_first[6]*mean(Z^2)+
                 betahat_first[7]*mean(Z^3))}}
    
    
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
  
  #example 1
  
  #OLS
  ggplot() +
    geom_ribbon(aes(x = X_est, ymin = q_lower_naive, ymax = q_upper_naive,fill="naive"),alpha=0.1) +
    geom_ribbon(aes(x = X_est, ymin = q_lower_est1, ymax = q_upper_est1,fill="RR"),alpha=0.8) +
    geom_ribbon(aes(x = X_est, ymin = q_lower_est3, ymax = q_upper_est3,fill="TSR"),alpha=0.8) +
    #geom_line(aes(x = X_est, y = apply(X_est,1,original), color = "E[Y|do(X=x)]"), linewidth = 1.3) +
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
    coord_cartesian(xlim = c(-15, 7.5), ylim = c(-30, 0)) + 
    geom_boxplot(aes(x = X, y = -30), width = 2.5, fill = "grey", color = "black", alpha = 0.7) +
    geom_boxplot(aes(x = X_S, y = -27), width = 2.5, fill = "grey", color = "black", alpha = 0.7) +
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
    #geom_line(aes(x = X_est, y = apply(X_est,1,original), color = "E[Y|do(X=x)]"), linewidth = 1.3) +
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
    coord_cartesian(xlim = c(-15, 7.5), ylim = c(-30, 0)) + 
    geom_boxplot(aes(x = X, y = -30), width = 2.5, fill = "grey", color = "black", alpha = 0.7) +
    geom_boxplot(aes(x = X_S, y = -27), width = 2.5, fill = "grey", color = "black", alpha = 0.7) +
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
  
  
  
  #example 2
  
  #OLS
  ggplot() +
    geom_ribbon(aes(x = X_est, ymin = q_lower_naive, ymax = q_upper_naive,fill="naive"),alpha=0.1) +
    geom_ribbon(aes(x = X_est, ymin = q_lower_est1, ymax = q_upper_est1,fill="RR"),alpha=0.8) +
    geom_ribbon(aes(x = X_est, ymin = q_lower_est3, ymax = q_upper_est3,fill="TSR"),alpha=0.8) +
    #geom_line(aes(x = X_est, y = apply(X_est,1,original), color = "E[Y|do(X=x)]"), linewidth = 1.3) +
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
    coord_cartesian(xlim = c(-15, 10), ylim = c(-80, 50)) + 
    geom_boxplot(aes(x = X, y = -77), width = 5, fill = "grey", color = "black", alpha = 0.7) +
    geom_boxplot(aes(x = X_S, y = -70), width = 5, fill = "grey", color = "black", alpha = 0.7) +
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
    #geom_line(aes(x = X_est, y = apply(X_est,1,original), color = "E[Y|do(X=x)]"), linewidth = 1.3) +
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
    coord_cartesian(xlim = c(-15, 10), ylim = c(-80, 50)) + 
    geom_boxplot(aes(x = X, y = -77), width = 5, fill = "grey", color = "black", alpha = 0.7) +
    geom_boxplot(aes(x = X_S, y = -70), width = 5, fill = "grey", color = "black", alpha = 0.7) +
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
  
  
  
  
  
  
  
  
  
  
  
  
  load("resultssim34.RData")  
  
  ####################################
  #######      bsp3     ########
  ####################################
  
  
  
  #MSE mean
  c(mean(result500_bsp3_rand$o_test_naive_S), mean(result500_bsp3_rand$o_test_est1_S), mean(result500_bsp3_rand$o_test_est1_ridge_S), mean(result500_bsp3_rand$o_test_est3_S), mean(result500_bsp3_rand$o_test_est3_ridge_S))
#  [1] 35.6449442 36.6156269 36.7504607  0.1205925  0.9988038
  c(mean(result1000_bsp3_rand$o_test_naive_S), mean(result1000_bsp3_rand$o_test_est1_S), mean(result1000_bsp3_rand$o_test_est1_ridge_S), mean(result1000_bsp3_rand$o_test_est3_S), mean(result1000_bsp3_rand$o_test_est3_ridge_S))
#  [1] 31.15328192 31.64197615 31.96345256  0.09667763  0.61421332
  c(mean(result5000_bsp3_rand$o_test_naive_S), mean(result5000_bsp3_rand$o_test_est1_S), mean(result5000_bsp3_rand$o_test_est1_ridge_S), mean(result5000_bsp3_rand$o_test_est3_S), mean(result5000_bsp3_rand$o_test_est3_ridge_S))
#  [1] 37.29704859 38.02651997 38.25615494  0.01472025  0.46376699
  
  c(mean(result500_bsp3_rand$o_test_naive), mean(result500_bsp3_rand$o_test_est1), mean(result500_bsp3_rand$o_test_est1_ridge), mean(result500_bsp3_rand$o_test_est3), mean(result500_bsp3_rand$o_test_est3_ridge))
#  [1] 53.5146597 31.6537186 91.8171665  0.9268416 71.8388169
  c(mean(result1000_bsp3_rand$o_test_naive), mean(result1000_bsp3_rand$o_test_est1), mean(result1000_bsp3_rand$o_test_est1_ridge), mean(result1000_bsp3_rand$o_test_est3), mean(result1000_bsp3_rand$o_test_est3_ridge))
#  [1] 34.6604681 27.0270784 77.9731433  0.4737679 61.7326935
  c(mean(result5000_bsp3_rand$o_test_naive), mean(result5000_bsp3_rand$o_test_est1), mean(result5000_bsp3_rand$o_test_est1_ridge), mean(result5000_bsp3_rand$o_test_est3), mean(result5000_bsp3_rand$o_test_est3_ridge))
#  [1] 30.13418824 32.59265494 68.97144740  0.07481121 47.65734062
  
  #MSE sd
  c(sd(result500_bsp3_rand$o_test_naive_S), sd(result500_bsp3_rand$o_test_est1_S), sd(result500_bsp3_rand$o_test_est1_ridge_S), sd(result500_bsp3_rand$o_test_est3_S), sd(result500_bsp3_rand$o_test_est3_ridge_S))
#  [1] 41.5423900 42.6082428 41.9460051  0.1671496  1.6230718
  c(sd(result1000_bsp3_rand$o_test_naive_S), sd(result1000_bsp3_rand$o_test_est1_S), sd(result1000_bsp3_rand$o_test_est1_ridge_S), sd(result1000_bsp3_rand$o_test_est3_S), sd(result1000_bsp3_rand$o_test_est3_ridge_S))
#  [1] 39.7572458 40.5089770 40.8560924  0.2469173  0.6546321
  c(sd(result5000_bsp3_rand$o_test_naive_S), sd(result5000_bsp3_rand$o_test_est1_S), sd(result5000_bsp3_rand$o_test_est1_ridge_S), sd(result5000_bsp3_rand$o_test_est3_S), sd(result5000_bsp3_rand$o_test_est3_ridge_S))
#  [1] 37.75715720 38.56842634 38.44537101  0.02576339  0.68475934
  
  c(sd(result500_bsp3_rand$o_test_naive), sd(result500_bsp3_rand$o_test_est1), sd(result500_bsp3_rand$o_test_est1_ridge), sd(result500_bsp3_rand$o_test_est3), sd(result500_bsp3_rand$o_test_est3_ridge))
#  [1] 80.566016 35.438583 82.686045  1.142754 84.497817
  c(sd(result1000_bsp3_rand$o_test_naive), sd(result1000_bsp3_rand$o_test_est1), sd(result1000_bsp3_rand$o_test_est1_ridge), sd(result1000_bsp3_rand$o_test_est3), sd(result1000_bsp3_rand$o_test_est3_ridge))
#  [1] 48.0693764 33.7713014 71.2139735  0.5997024 74.9253923
  c(sd(result5000_bsp3_rand$o_test_naive), sd(result5000_bsp3_rand$o_test_est1), sd(result5000_bsp3_rand$o_test_est1_ridge), sd(result5000_bsp3_rand$o_test_est3), sd(result5000_bsp3_rand$o_test_est3_ridge))
#  [1] 30.03080394 32.57438263 62.86892344  0.08046236 65.08325328
  
  
  #################
  ###   empty   ###
  #################
  
  
  
  #MSE mean
  c(mean(empty_result500_bsp3_rand$o_test_naive_S), mean(empty_result500_bsp3_rand$o_test_est1_S), mean(empty_result500_bsp3_rand$o_test_est1_ridge_S), mean(empty_result500_bsp3_rand$o_test_est3_S), mean(empty_result500_bsp3_rand$o_test_est3_ridge_S))
#  [1] 35.6744443 36.5361819 36.6901057  0.1670062  0.9643019
  c(mean(empty_result1000_bsp3_rand$o_test_naive_S), mean(empty_result1000_bsp3_rand$o_test_est1_S), mean(empty_result1000_bsp3_rand$o_test_est1_ridge_S), mean(empty_result1000_bsp3_rand$o_test_est3_S), mean(empty_result1000_bsp3_rand$o_test_est3_ridge_S))
#  [1] 31.41031767 31.26722172 32.13499213  0.08513438  0.90782797
  c(mean(empty_result5000_bsp3_rand$o_test_naive_S), mean(empty_result5000_bsp3_rand$o_test_est1_S), mean(empty_result5000_bsp3_rand$o_test_est1_ridge_S), mean(empty_result5000_bsp3_rand$o_test_est3_S), mean(empty_result5000_bsp3_rand$o_test_est3_ridge_S))
#  [1] 37.23974456 37.76166317 37.86313006  0.01502261  0.46462773
  
  c(mean(empty_result500_bsp3_rand$o_test_naive), mean(empty_result500_bsp3_rand$o_test_est1), mean(empty_result500_bsp3_rand$o_test_est1_ridge), mean(empty_result500_bsp3_rand$o_test_est3), mean(empty_result500_bsp3_rand$o_test_est3_ridge))
#  [1] 47.7898994 31.7713384 90.0785818  0.8995274 69.9116680
  c(mean(empty_result1000_bsp3_rand$o_test_naive), mean(empty_result1000_bsp3_rand$o_test_est1), mean(empty_result1000_bsp3_rand$o_test_est1_ridge), mean(empty_result1000_bsp3_rand$o_test_est3), mean(empty_result1000_bsp3_rand$o_test_est3_ridge))
#  [1] 34.7408078 26.7578646 74.1243652  0.4677581 57.1552128
  c(mean(empty_result5000_bsp3_rand$o_test_naive), mean(empty_result5000_bsp3_rand$o_test_est1), mean(empty_result5000_bsp3_rand$o_test_est1_ridge), mean(empty_result5000_bsp3_rand$o_test_est3), mean(empty_result5000_bsp3_rand$o_test_est3_ridge))
#  [1] 30.29553818 32.46810439 69.13614122  0.07139673 47.85212766
  
  #MSE sd
  c(sd(empty_result500_bsp3_rand$o_test_naive_S), sd(empty_result500_bsp3_rand$o_test_est1_S), sd(empty_result500_bsp3_rand$o_test_est1_ridge_S), sd(empty_result500_bsp3_rand$o_test_est3_S), sd(empty_result500_bsp3_rand$o_test_est3_ridge_S))
#  [1] 44.1738445 43.8310978 43.7532661  0.2512945  1.4802920  
  c(sd(empty_result1000_bsp3_rand$o_test_naive_S), sd(empty_result1000_bsp3_rand$o_test_est1_S), sd(empty_result1000_bsp3_rand$o_test_est1_ridge_S), sd(empty_result1000_bsp3_rand$o_test_est3_S), sd(empty_result1000_bsp3_rand$o_test_est3_ridge_S))
#  [1] 39.7455674 38.9143439 39.3206782  0.1347602  1.5097884
  c(sd(empty_result5000_bsp3_rand$o_test_naive_S), sd(empty_result5000_bsp3_rand$o_test_est1_S), sd(empty_result5000_bsp3_rand$o_test_est1_ridge_S), sd(empty_result5000_bsp3_rand$o_test_est3_S), sd(empty_result5000_bsp3_rand$o_test_est3_ridge_S))
#  [1] 37.69850148 38.29836316 38.14130542  0.02085728  0.65052555
  
  c(sd(empty_result500_bsp3_rand$o_test_naive), sd(empty_result500_bsp3_rand$o_test_est1), sd(empty_result500_bsp3_rand$o_test_est1_ridge), sd(empty_result500_bsp3_rand$o_test_est3), sd(empty_result500_bsp3_rand$o_test_est3_ridge))
#  [1] 63.768694 36.686498 80.879262  1.063708 81.750786
  c(sd(empty_result1000_bsp3_rand$o_test_naive), sd(empty_result1000_bsp3_rand$o_test_est1), sd(empty_result1000_bsp3_rand$o_test_est1_ridge), sd(empty_result1000_bsp3_rand$o_test_est3), sd(empty_result1000_bsp3_rand$o_test_est3_ridge))
#  [1] 47.4513032 32.7568960 61.2831525  0.5712976 61.1681450
  c(sd(empty_result5000_bsp3_rand$o_test_naive), sd(empty_result5000_bsp3_rand$o_test_est1), sd(empty_result5000_bsp3_rand$o_test_est1_ridge), sd(empty_result5000_bsp3_rand$o_test_est3), sd(empty_result5000_bsp3_rand$o_test_est3_ridge))
#  [1] 30.32637472 32.50128824 63.63360537  0.07673282 65.57723465
  
  
  
  
  
  ####################################
  #######      bsp4     ########
  ####################################
  
  
  
  #MSE mean
  c(mean(result500_bsp4_rand$o_test_naive_S), mean(result500_bsp4_rand$o_test_est1_S), mean(result500_bsp4_rand$o_test_est1_ridge_S), mean(result500_bsp4_rand$o_test_est3_S), mean(result500_bsp4_rand$o_test_est3_ridge_S))
#  [1] 44.6227494 62.8479430 62.3944877  0.5764024  0.5798602
  c(mean(result1000_bsp4_rand$o_test_naive_S), mean(result1000_bsp4_rand$o_test_est1_S), mean(result1000_bsp4_rand$o_test_est1_ridge_S), mean(result1000_bsp4_rand$o_test_est3_S), mean(result1000_bsp4_rand$o_test_est3_ridge_S))
#  [1] 59.0884032 82.6869440 82.4419531  0.3541259  0.3533019
  c(mean(result5000_bsp4_rand$o_test_naive_S), mean(result5000_bsp4_rand$o_test_est1_S), mean(result5000_bsp4_rand$o_test_est1_ridge_S), mean(result5000_bsp4_rand$o_test_est3_S), mean(result5000_bsp4_rand$o_test_est3_ridge_S))
#  [1] 49.1552437 71.7209196 71.5891267  0.0456946  0.0456709
  
  c(mean(result500_bsp4_rand$o_test_naive), mean(result500_bsp4_rand$o_test_est1), mean(result500_bsp4_rand$o_test_est1_ridge), mean(result500_bsp4_rand$o_test_est3), mean(result500_bsp4_rand$o_test_est3_ridge))
#  [1] 192.911104 180.606911 180.716317   2.425420   2.432447
  c(mean(result1000_bsp4_rand$o_test_naive), mean(result1000_bsp4_rand$o_test_est1), mean(result1000_bsp4_rand$o_test_est1_ridge), mean(result1000_bsp4_rand$o_test_est3), mean(result1000_bsp4_rand$o_test_est3_ridge))
#  [1] 232.6721966 218.9854551 219.2957003   0.8059371   0.8222287
  c(mean(result5000_bsp4_rand$o_test_naive), mean(result5000_bsp4_rand$o_test_est1), mean(result5000_bsp4_rand$o_test_est1_ridge), mean(result5000_bsp4_rand$o_test_est3), mean(result5000_bsp4_rand$o_test_est3_ridge))
#  [1] 178.8953689 183.8282702 183.8122237   0.1156439   0.1236903
  
  
  #MSE sd
  c(sd(result500_bsp4_rand$o_test_naive_S), sd(result500_bsp4_rand$o_test_est1_S), sd(result500_bsp4_rand$o_test_est1_ridge_S), sd(result500_bsp4_rand$o_test_est3_S), sd(result500_bsp4_rand$o_test_est3_ridge_S))
#  [1] 39.0775782 55.6331986 55.1331476  0.7745431  0.7725040
  c(sd(result1000_bsp4_rand$o_test_naive_S), sd(result1000_bsp4_rand$o_test_est1_S), sd(result1000_bsp4_rand$o_test_est1_ridge_S), sd(result1000_bsp4_rand$o_test_est3_S), sd(result1000_bsp4_rand$o_test_est3_ridge_S))
#  [1] 55.4142989 62.0562508 61.8427159  0.6853223  0.6959319
  c(sd(result5000_bsp4_rand$o_test_naive_S), sd(result5000_bsp4_rand$o_test_est1_S), sd(result5000_bsp4_rand$o_test_est1_ridge_S), sd(result5000_bsp4_rand$o_test_est3_S), sd(result5000_bsp4_rand$o_test_est3_ridge_S))
#  [1] 40.60831589 58.77802319 58.67472023  0.05976307  0.05924058
  
  c(sd(result500_bsp4_rand$o_test_naive), sd(result500_bsp4_rand$o_test_est1), sd(result500_bsp4_rand$o_test_est1_ridge), sd(result500_bsp4_rand$o_test_est3), sd(result500_bsp4_rand$o_test_est3_ridge))
#  [1] 277.920941 190.545146 190.174591   3.862193   4.014499
  c(sd(result1000_bsp4_rand$o_test_naive), sd(result1000_bsp4_rand$o_test_est1), sd(result1000_bsp4_rand$o_test_est1_ridge), sd(result1000_bsp4_rand$o_test_est3), sd(result1000_bsp4_rand$o_test_est3_ridge))
#  [1] 314.833494 181.332794 181.574523   1.425579   1.370006
  c(sd(result5000_bsp4_rand$o_test_naive), sd(result5000_bsp4_rand$o_test_est1), sd(result5000_bsp4_rand$o_test_est1_ridge), sd(result5000_bsp4_rand$o_test_est3), sd(result5000_bsp4_rand$o_test_est3_ridge))
#  [1] 231.5467994 172.6455132 172.7924880   0.1073011   0.1189043
  
  
  
  #################
  ###   empty   ###
  #################
  
  
  
  #MSE mean
  c(mean(empty_result500_bsp4_rand$o_test_naive_S), mean(empty_result500_bsp4_rand$o_test_est1_S), mean(empty_result500_bsp4_rand$o_test_est1_ridge_S), mean(empty_result500_bsp4_rand$o_test_est3_S), mean(empty_result500_bsp4_rand$o_test_est3_ridge_S))
#  [1] 43.8826037 61.2624617 60.8088530  0.5540878  0.5643643
  c(mean(empty_result1000_bsp4_rand$o_test_naive_S), mean(empty_result1000_bsp4_rand$o_test_est1_S), mean(empty_result1000_bsp4_rand$o_test_est1_ridge_S), mean(empty_result1000_bsp4_rand$o_test_est3_S), mean(empty_result1000_bsp4_rand$o_test_est3_ridge_S))
#  [1] 60.2507445 85.4867522 85.1975652  0.3918808  0.3962078
  c(mean(empty_result5000_bsp4_rand$o_test_naive_S), mean(empty_result5000_bsp4_rand$o_test_est1_S), mean(empty_result5000_bsp4_rand$o_test_est1_ridge_S), mean(empty_result5000_bsp4_rand$o_test_est3_S), mean(empty_result5000_bsp4_rand$o_test_est3_ridge_S))
#  [1] 48.07364646 70.73879902 70.62531847  0.04364758  0.04401718
  
  c(mean(empty_result500_bsp4_rand$o_test_naive), mean(empty_result500_bsp4_rand$o_test_est1), mean(empty_result500_bsp4_rand$o_test_est1_ridge), mean(empty_result500_bsp4_rand$o_test_est3), mean(empty_result500_bsp4_rand$o_test_est3_ridge))
#  [1] 218.908769 177.403191 176.754937   1.783956   1.870122
  c(mean(empty_result1000_bsp4_rand$o_test_naive), mean(empty_result1000_bsp4_rand$o_test_est1), mean(empty_result1000_bsp4_rand$o_test_est1_ridge), mean(empty_result1000_bsp4_rand$o_test_est3), mean(empty_result1000_bsp4_rand$o_test_est3_ridge))
#  [1] 205.7929574 228.0138896 228.1654573   0.9859473   1.0800008
  c(mean(empty_result5000_bsp4_rand$o_test_naive), mean(empty_result5000_bsp4_rand$o_test_est1), mean(empty_result5000_bsp4_rand$o_test_est1_ridge), mean(empty_result5000_bsp4_rand$o_test_est3), mean(empty_result5000_bsp4_rand$o_test_est3_ridge))
#  [1] 158.3365402 185.2388034 185.2908046   0.1186099   0.1235730
  
  #MSE sd
  c(sd(empty_result500_bsp4_rand$o_test_naive_S), sd(empty_result500_bsp4_rand$o_test_est1_S), sd(empty_result500_bsp4_rand$o_test_est1_ridge_S), sd(empty_result500_bsp4_rand$o_test_est3_S), sd(empty_result500_bsp4_rand$o_test_est3_ridge_S))
#  [1] 40.2648141 54.9784906 54.5550490  0.8853187  0.8762496  
  c(sd(empty_result1000_bsp4_rand$o_test_naive_S), sd(empty_result1000_bsp4_rand$o_test_est1_S), sd(empty_result1000_bsp4_rand$o_test_est1_ridge_S), sd(empty_result1000_bsp4_rand$o_test_est3_S), sd(empty_result1000_bsp4_rand$o_test_est3_ridge_S))
#  [1] 50.0424325 66.6186601 66.4067656  0.7588129  0.7555894
  c(sd(empty_result5000_bsp4_rand$o_test_naive_S), sd(empty_result5000_bsp4_rand$o_test_est1_S), sd(empty_result5000_bsp4_rand$o_test_est1_ridge_S), sd(empty_result5000_bsp4_rand$o_test_est3_S), sd(empty_result5000_bsp4_rand$o_test_est3_ridge_S))
#  [1] 38.90925425 58.99898156 58.91292711  0.05930191  0.05922655
  
  c(sd(empty_result500_bsp4_rand$o_test_naive), sd(empty_result500_bsp4_rand$o_test_est1), sd(empty_result500_bsp4_rand$o_test_est1_ridge), sd(empty_result500_bsp4_rand$o_test_est3), sd(empty_result500_bsp4_rand$o_test_est3_ridge))
#  [1] 362.722719 180.897580 179.734133   2.115034   2.448291
  c(sd(empty_result1000_bsp4_rand$o_test_naive), sd(empty_result1000_bsp4_rand$o_test_est1), sd(empty_result1000_bsp4_rand$o_test_est1_ridge), sd(empty_result1000_bsp4_rand$o_test_est3), sd(empty_result1000_bsp4_rand$o_test_est3_ridge))
#  [1] 257.185461 190.308922 190.809707   1.171618   1.229523
  c(sd(empty_result5000_bsp4_rand$o_test_naive), sd(empty_result5000_bsp4_rand$o_test_est1), sd(empty_result5000_bsp4_rand$o_test_est1_ridge), sd(empty_result5000_bsp4_rand$o_test_est3), sd(empty_result5000_bsp4_rand$o_test_est3_ridge))
#  [1] 167.9747360 178.3857117 178.5123567   0.1146284   0.1227108
    
