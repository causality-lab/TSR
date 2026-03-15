

#Section 5.1 randomized coefficients

#load packages
library(glmnet)
library(ggplot2)

#n = sample size
#SIM = number of simulation runs
#empty = indicates if ((S cap D)=emptyset)
#linear = indicates linear case

sim1<-function(n,SIM,empty,linear)
{
  #create empty matrices to save results
  
  vars_test <- c("o_test_naive", "o_test_est1", "o_test_est1_ridge", "o_test_est3", "o_test_est3_ridge")
  lapply(vars_test, function(x) assign(x, NULL, envir = .GlobalEnv))
  
  vars_test_S <- c("o_test_naive_S", "o_test_est1_S", "o_test_est1_ridge_S", "o_test_est3_S", "o_test_est3_ridge_S")
  lapply(vars_test_S, function(x) assign(x, NULL, envir = .GlobalEnv))
  
  bias_vars_test <- c("bias_o_test_naive", "bias_o_test_est1", "bias_o_test_est1_ridge", "bias_o_test_est3", "bias_o_test_est3_ridge")
  lapply(bias_vars_test, function(x) assign(x, NULL, envir = .GlobalEnv))
  
  bias_vars_test_S <- c("bias_o_test_naive_S", "bias_o_test_est1_S", "bias_o_test_est1_ridge_S", "bias_o_test_est3_S", "bias_o_test_est3_ridge_S")
  lapply(bias_vars_test_S, function(x) assign(x, NULL, envir = .GlobalEnv))
  
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
    Z_Ddata_complete <- matrix(ncol=SIM,nrow=n)}
  
 
    betahat_first_complete <- matrix(ncol=SIM,nrow=5)
    betahat_first_ridge_complete <- matrix(ncol=SIM,nrow=5)
    beta_end1_complete <- matrix(ncol=SIM, nrow=3)
    beta_end1_ridge_complete <- matrix(ncol=SIM, nrow=3)
    naive_beta_complete <- matrix(ncol=SIM, nrow=3)

  y0_complete<-NULL
  yx1_complete<-NULL
  yx2_complete<-NULL
  yz1_complete<-NULL
  yz2_complete<-NULL
  
  for(sim in 1:SIM)
  {
    print(sim)
    
    #reproducibility
    set.seed(sim)
    
    #data generating processes
    if(empty==F){
      X<-rnorm(n)
      Z<-rnorm(n,mean=-2)
      S<-X+Z < -2
      if(linear==F){
        y0<-runif(1)-0.5
        yx1<-runif(1)-0.5
        yx2<-runif(1)+2.5
        yz1<-runif(1)+4.5
        yz2<-runif(1)-0.5
        Y<-y0+yx1*X+yx2*X^2+yz1*Z+yz2*Z^2+rnorm(n)
      }else{
        y0<-runif(1)-0.5
        yx1<-runif(1)+2.5
        yx2<-runif(1)-0.5
        yz1<-runif(1)+4.5
        yz2<-runif(1)-0.5
        Y<-y0+yx1*X+yx2*X^2+yz1*Z+yz2*Z^2+rnorm(n)
      }
      Z_S<-Z[S==1]
      X_S<-X[S==1]
      Y_S<-Y[S==1]
    }else{
      X_Sdata<-rnorm(n)
      Z_Sdata<-rnorm(n,mean=-2)
      S_Sdata<-X_Sdata+Z_Sdata < -2
      if(linear==F){
        y0<-runif(1)-0.5
        yx1<-runif(1)-0.5
        yx2<-runif(1)+2.5
        yz1<-runif(1)+4.5
        yz2<-runif(1)-0.5
        Y_Sdata<-y0+yx1*X_Sdata+yx2*X_Sdata^2+yz1*Z_Sdata+yz2*Z_Sdata^2+rnorm(n)
      }else{
        y0<-runif(1)-0.5
        yx1<-runif(1)+2.5
        yx2<-runif(1)-0.5
        yz1<-runif(1)+4.5
        yz2<-runif(1)-0.5
        Y_Sdata<-y0+yx1*X_Sdata+yx2*X_Sdata^2+yz1*Z_Sdata+yz2*Z_Sdata^2+rnorm(n)
      }
      Z_S<-Z_Sdata[S_Sdata==1]
      X_S<-X_Sdata[S_Sdata==1]
      Y_S<-Y_Sdata[S_Sdata==1]
      
      X_Ddata<-rnorm(n)
      Z_Ddata<-rnorm(n,mean=-2)
      
      X<-X_Ddata
      Z<-Z_Ddata}
    
    #save realisations from all simulation runs
    if(empty==F){
      X_complete[,sim]<-X
      Z_complete[,sim]<-Z
      S_complete[,sim]<-S
      Y_complete[,sim]<-Y
    }else{ X_Sdata_complete[,sim] <- X_Sdata
    Z_Sdata_complete[,sim] <- Z_Sdata
    S_Sdata_complete[,sim] <- S_Sdata
    Y_Sdata_complete[,sim] <- Y_Sdata
    
    X_Ddata_complete[,sim] <- X_Ddata
    Z_Ddata_complete[,sim] <- Z_Ddata}
    
     
      #first step ridge 
      lambda_seq <- 10^seq(2, -2, by = -.1)
      ridge_cv <- cv.glmnet(cbind(X_S,I(X_S^2),Z_S,I(Z_S^2)), Y_S, alpha = 0, lambda = lambda_seq)
      best_lambda <- ridge_cv$lambda.min
      best_ridge <- glmnet(cbind(X_S,I(X_S^2),Z_S,I(Z_S^2)), Y_S, alpha = 0, lambda = best_lambda)
      betahat_first_ridge<-as.vector(coef(best_ridge))
      
      #first step OLS
      lm_fit <- lm(Y_S~X_S+I(X_S^2)+Z_S+I(Z_S^2))
      betahat_first<-coef(lm_fit)
      
      #calculate target variable for second step of RR
      Y_second_ridge<-(cbind(rep(1,n),X,X^2,Z,Z^2))%*%as.vector(betahat_first_ridge)
      Y_second<-(cbind(rep(1,n),X,X^2,Z,Z^2))%*%as.vector(betahat_first)
      
      #final RR estimates
      beta_end1_ridge<-lm(Y_second_ridge~1+X+I(X^2))$coef
      beta_end1<-lm(Y_second~1+X+I(X^2))$coef
      
      #naive estimation only based on S=1
      naive_beta<-lm(Y_S~1+X_S+I(X_S^2))$coef
      
    
      #E[Y|do(X)] = E[Y|X]
      original<-function(x)
      {y0+yx1*x+yx2*x^2+yz1*(-2)+yz2*5}
      
      #RR ridge
      est1_ridge<-function(x)
      {beta_end1_ridge[1]+beta_end1_ridge[2]*x+beta_end1_ridge[3]*x^2}
      
      #RR OLS
      est1<-function(x)
      {beta_end1[1]+beta_end1[2]*x+beta_end1[3]*x^2}
      
      #naive
      naive<-function(x)
      {naive_beta[1]+naive_beta[2]*x+naive_beta[3]*x^2}
      
      #TSR ridge
      est3_ridge<-function(x)
      {return(betahat_first_ridge[1]+betahat_first_ridge[2]*x+betahat_first_ridge[3]*x^2+
                betahat_first_ridge[4]*mean(Z)+
                betahat_first_ridge[5]*mean(Z^2))}
      
      #TSR OLS
      est3<-function(x)
      {return(betahat_first[1]+betahat_first[2]*x+betahat_first[3]*x^2+
                betahat_first[4]*mean(Z)+
                betahat_first[5]*mean(Z^2))}
      
    #save coefficient vectors from all simulation runs
    betahat_first_complete[,sim] <- betahat_first
    betahat_first_ridge_complete[,sim] <- betahat_first_ridge
    beta_end1_complete[,sim] <- beta_end1
    beta_end1_ridge_complete[,sim] <- beta_end1_ridge
    naive_beta_complete[,sim] <- naive_beta
    
    #test data  
    X_test<-as.matrix(rnorm(n=n))
    Z_test<-rnorm(n,mean=-2)
    S_test<-X_test+Z_test< -2
    if(linear==F){
      Y_test<-y0+yx1*X_test+yx2*X_test^2+yz1*Z_test+yz2*Z_test^2
    }else{
      Y_test<-y0+yx1*X_test+yx2*X_test^2+yz1*Z_test+yz2*Z_test^2
    }
    X_test_S<-as.matrix(X_test[S_test==1,])
    Y_test_S<-Y_test[S_test==1]
    
    #MSE (o: E[Y|do(X)] = E[Y|X]) 
    
    #in D
    o_test_naive[sim]<-mean((apply(X_test,1,original)-apply(X_test,1,naive))^2)
    o_test_est1[sim]<-mean((apply(X_test,1,original)-apply(X_test,1,est1))^2)
    o_test_est1_ridge[sim]<-mean((apply(X_test,1,original)-apply(X_test,1,est1_ridge))^2)
    o_test_est3[sim]<-mean((apply(X_test,1,original)-apply(X_test,1,est3))^2)
    o_test_est3_ridge[sim]<-mean((apply(X_test,1,original)-apply(X_test,1,est3_ridge))^2)
    
    #in S
    o_test_naive_S[sim]<-mean((apply(X_test_S,1,original)-apply(X_test_S,1,naive))^2)
    o_test_est1_S[sim]<-mean((apply(X_test_S,1,original)-apply(X_test_S,1,est1))^2)
    o_test_est1_ridge_S[sim]<-mean((apply(X_test_S,1,original)-apply(X_test_S,1,est1_ridge))^2)
    o_test_est3_S[sim]<-mean((apply(X_test_S,1,original)-apply(X_test_S,1,est3))^2)
    o_test_est3_ridge_S[sim]<-mean((apply(X_test_S,1,original)-apply(X_test_S,1,est3_ridge))^2)
    
    #Bias (o: E[Y|do(X)] = E[Y|X]) 
    
    #in D
    bias_o_test_naive[sim]<-mean((apply(X_test,1,original)-apply(X_test,1,naive)))
    bias_o_test_est1[sim]<-mean((apply(X_test,1,original)-apply(X_test,1,est1)))
    bias_o_test_est1_ridge[sim]<-mean((apply(X_test,1,original)-apply(X_test,1,est1_ridge)))
    bias_o_test_est3[sim]<-mean((apply(X_test,1,original)-apply(X_test,1,est3)))
    bias_o_test_est3_ridge[sim]<-mean((apply(X_test,1,original)-apply(X_test,1,est3_ridge)))
    
    #in S
    bias_o_test_naive_S[sim]<-mean((apply(X_test_S,1,original)-apply(X_test_S,1,naive)))
    bias_o_test_est1_S[sim]<-mean((apply(X_test_S,1,original)-apply(X_test_S,1,est1)))
    bias_o_test_est1_ridge_S[sim]<-mean((apply(X_test_S,1,original)-apply(X_test_S,1,est1_ridge)))
    bias_o_test_est3_S[sim]<-mean((apply(X_test_S,1,original)-apply(X_test_S,1,est3)))
    bias_o_test_est3_ridge_S[sim]<-mean((apply(X_test_S,1,original)-apply(X_test_S,1,est3_ridge)))
  }
  
  y0_complete[sim]<-y0
  yx1_complete[sim]<-yx1
  yx2_complete[sim]<-yx2
  yz1_complete[sim]<-yz1
  yz2_complete[sim]<-yz2
  
  #return the relevant results  
  if(empty==F){
    return(list(SIM = SIM, n = n, empty = empty, linear = linear,
                X_complete = X_complete, Z_complete = Z_complete, S_complete = S_complete, Y_complete = Y_complete,
                betahat_first_complete = betahat_first_complete, betahat_first_ridge_complete = betahat_first_ridge_complete,
                beta_end1_complete = beta_end1_complete, beta_end1_ridge_complete = beta_end1_ridge_complete,
                naive_beta_complete = naive_beta_complete,
                o_test_naive=o_test_naive, o_test_est1=o_test_est1, o_test_est1_ridge=o_test_est1_ridge, o_test_est3=o_test_est3, o_test_est3_ridge=o_test_est3_ridge,
                o_test_naive_S=o_test_naive_S, o_test_est1_S=o_test_est1_S, o_test_est1_ridge_S=o_test_est1_ridge_S, o_test_est3_S=o_test_est3_S, o_test_est3_ridge_S=o_test_est3_ridge_S,
                bias_o_test_naive=bias_o_test_naive, bias_o_test_est1=bias_o_test_est1, bias_o_test_est1_ridge=bias_o_test_est1_ridge, bias_o_test_est3=bias_o_test_est3, bias_o_test_est3_ridge=bias_o_test_est3_ridge,
                bias_o_test_naive_S=bias_o_test_naive_S, bias_o_test_est1_S=bias_o_test_est1_S, bias_o_test_est1_ridge_S=bias_o_test_est1_ridge_S, bias_o_test_est3_S=bias_o_test_est3_S, bias_o_test_est3_ridge_S=bias_o_test_est3_ridge_S,
                y0=y0,yx1=yx1,yx2=yx2,yz1=yz1,yz2=yz2))
  }else{
    return(list(SIM = SIM, n = n, empty = empty, linear = linear,
                X_Sdata_complete = X_Sdata_complete, Z_Sdata_complete = Z_Sdata_complete, S_Sdata_complete = S_Sdata_complete, Y_Sdata_complete = Y_Sdata_complete,
                X_Ddata_complete = X_Ddata_complete, Z_Ddata_complete = Z_Ddata_complete,
                betahat_first_complete = betahat_first_complete, betahat_first_ridge_complete = betahat_first_ridge_complete,
                beta_end1_complete = beta_end1_complete, beta_end1_ridge_complete = beta_end1_ridge_complete,
                naive_beta_complete = naive_beta_complete,
                o_test_naive=o_test_naive, o_test_est1=o_test_est1, o_test_est1_ridge=o_test_est1_ridge, o_test_est3=o_test_est3, o_test_est3_ridge=o_test_est3_ridge,
                o_test_naive_S=o_test_naive_S, o_test_est1_S=o_test_est1_S, o_test_est1_ridge_S=o_test_est1_ridge_S, o_test_est3_S=o_test_est3_S, o_test_est3_ridge_S=o_test_est3_ridge_S,
                bias_o_test_naive=bias_o_test_naive, bias_o_test_est1=bias_o_test_est1, bias_o_test_est1_ridge=bias_o_test_est1_ridge, bias_o_test_est3=bias_o_test_est3, bias_o_test_est3_ridge=bias_o_test_est3_ridge,
                bias_o_test_naive_S=bias_o_test_naive_S, bias_o_test_est1_S=bias_o_test_est1_S, bias_o_test_est1_ridge_S=bias_o_test_est1_ridge_S, bias_o_test_est3_S=bias_o_test_est3_S, bias_o_test_est3_ridge_S=bias_o_test_est3_ridge_S,
                y0=y0,yx1=yx1,yx2=yx2,yz1=yz1,yz2=yz2))}
  
}



#####################################################################################################
#####################################################################################################
#####################################################################################################

result500_2 <- sim1(n=500,SIM=100,empty=F,linear=F)
result1000_2 <- sim1(n=1000,SIM=100,empty=F,linear=F)
result5000_2 <- sim1(n=5000,SIM=100,empty=F,linear=F)

empty_result500_2 <- sim1(n=500,SIM=100,empty=T,linear=F)
empty_result1000_2 <- sim1(n=1000,SIM=100,empty=T,linear=F)
empty_result5000_2 <- sim1(n=5000,SIM=100,empty=T,linear=F)

linear_result500_2 <- sim1(n=500,SIM=100,empty=F,linear=T)
linear_result1000_2 <- sim1(n=1000,SIM=100,empty=F,linear=T)
linear_result5000_2 <- sim1(n=5000,SIM=100,empty=F,linear=T)

empty_linear_result500_2 <- sim1(n=500,SIM=100,empty=T,linear=T)
empty_linear_result1000_2 <- sim1(n=1000,SIM=100,empty=T,linear=T)
empty_linear_result5000_2 <- sim1(n=5000,SIM=100,empty=T,linear=T)



result500_plot<-result500_2 
result1000_plot<-result1000_2 
result5000_plot<-result5000_2 

result500_plot<-empty_result500_2 
result1000_plot<-empty_result1000_2 
result5000_plot<-empty_result5000_2 

result500_plot<-linear_result500_2 
result1000_plot<-linear_result1000_2 
result5000_plot<-linear_result5000_2 

result500_plot<-empty_linear_result500_2 
result1000_plot<-empty_linear_result1000_2 
result5000_plot<-empty_linear_result5000_2 

#MSE

par(mfrow=c(1,1),cex.axis = 1.4,cex.lab=2.5,oma=c(0,0.5,0,0),mar=c(7,7,3,3))
groups=c("n=500","n=500","n=1000","n=1000","n=5000","n=5000")
boxplot(ylab="MSE",cbind(result500_plot$o_test_est1,result500_plot$o_test_est3,result1000_plot$o_test_est1,result1000_plot$o_test_est3,result5000_plot$o_test_est1,result5000_plot$o_test_est3),ylim=c(0,0.9),pch=20,col=c("lightblue","lightgreen"),xaxt = "n")
axis(1, at = c(1.5, 3.5, 5.5), labels = c("n=500", "n=1000", "n=5000"))
legend(legend=c("RR","TSR"),col=c("lightblue","lightgreen"),"topright",pch=15,cex=1.7)
boxplot(ylab="MSE",cbind(result500_plot$o_test_est1_S,result500_plot$o_test_est3_S,result1000_plot$o_test_est1_S,result1000_plot$o_test_est3_S,result5000_plot$o_test_est1_S,result5000_plot$o_test_est3_S),ylim=c(0,0.9),pch=20,col=c("lightblue","lightgreen"),xaxt = "n")
axis(1, at = c(1.5, 3.5, 5.5), labels = c("n=500", "n=1000", "n=5000"))
legend(legend=c("RR","TSR"),col=c("lightblue","lightgreen"),"topright",pch=15,cex=1.7)

groups=c("n=500","n=500","n=1000","n=1000","n=5000","n=5000")
boxplot(ylab="MSE",cbind(result500_plot$o_test_est1_ridge,result500_plot$o_test_est3_ridge,result1000_plot$o_test_est1_ridge,result1000_plot$o_test_est3_ridge,result5000_plot$o_test_est1_ridge,result5000_plot$o_test_est3_ridge),ylim=c(0,0.9),pch=20,col=c("lightblue","lightgreen"),xaxt = "n")
axis(1, at = c(1.5, 3.5, 5.5), labels = c("n=500", "n=1000", "n=5000"))
legend(legend=c("RR (ridge)","TSR (ridge)"),col=c("lightblue","lightgreen"),"topright",pch=15,cex=1.7)
boxplot(ylab="MSE",cbind(result500_plot$o_test_est1_ridge_S,result500_plot$o_test_est3_ridge_S,result1000_plot$o_test_est1_ridge_S,result1000_plot$o_test_est3_ridge_S,result5000_plot$o_test_est1_ridge_S,result5000_plot$o_test_est3_ridge_S),ylim=c(0,0.9),pch=20,col=c("lightblue","lightgreen"),xaxt = "n")
axis(1, at = c(1.5, 3.5, 5.5), labels = c("n=500", "n=1000", "n=5000"))
legend(legend=c("RR (ridge)","TSR (ridge)"),col=c("lightblue","lightgreen"),"topright",pch=15,cex=1.7)


#save(
#  
#  result500_2, 
#  result1000_2, 
#  result5000_2, 
#  
#  empty_result500_2, 
#  empty_result1000_2, 
#  empty_result5000_2, 
#  
#  linear_result500_2, 
#  linear_result1000_2, 
#  linear_result5000_2, 
#  
#  empty_linear_result500_2, 
#  empty_linear_result1000_2, 
#  empty_linear_result5000_2, 
#  
#  file="boxplots_sec41_2.Rdata") 


load("boxplots_sec41_2.Rdata")




sim1<-function(n,SIM,empty,linear)
{
  #create empty matrices to save results
  
  vars_test <- c("o_test_naive", "o_test_est1", "o_test_est1_ridge", "o_test_est3", "o_test_est3_ridge")
  lapply(vars_test, function(x) assign(x, NULL, envir = .GlobalEnv))
  
  vars_test_S <- c("o_test_naive_S", "o_test_est1_S", "o_test_est1_ridge_S", "o_test_est3_S", "o_test_est3_ridge_S")
  lapply(vars_test_S, function(x) assign(x, NULL, envir = .GlobalEnv))
  
  bias_vars_test <- c("bias_o_test_naive", "bias_o_test_est1", "bias_o_test_est1_ridge", "bias_o_test_est3", "bias_o_test_est3_ridge")
  lapply(bias_vars_test, function(x) assign(x, NULL, envir = .GlobalEnv))
  
  bias_vars_test_S <- c("bias_o_test_naive_S", "bias_o_test_est1_S", "bias_o_test_est1_ridge_S", "bias_o_test_est3_S", "bias_o_test_est3_ridge_S")
  lapply(bias_vars_test_S, function(x) assign(x, NULL, envir = .GlobalEnv))
  
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
    Z_Ddata_complete <- matrix(ncol=SIM,nrow=n)}
  
  
  betahat_first_complete <- matrix(ncol=SIM,nrow=5)
  betahat_first_ridge_complete <- matrix(ncol=SIM,nrow=5)
  beta_end1_complete <- matrix(ncol=SIM, nrow=3)
  beta_end1_ridge_complete <- matrix(ncol=SIM, nrow=3)
  naive_beta_complete <- matrix(ncol=SIM, nrow=3)
  
  y0_complete<-NULL
  yx1_complete<-NULL
  yx2_complete<-NULL
  yz1_complete<-NULL
  yz2_complete<-NULL
  
  for(sim in 1:SIM)
  {
    print(sim)
    
    #reproducibility
    set.seed(sim)
    
    #data generating processes
    if(empty==F){
      X<-rnorm(n)
      Z<-rnorm(n,mean=-2)
      S<-X+Z < -2
      if(linear==F){
        y0<-0.5-0.5
        yx1<-0.5-0.5
        yx2<-0.5+2.5
        yz1<-0.5+4.5
        yz2<-0.5-0.5
        Y<-y0+yx1*X+yx2*X^2+yz1*Z+yz2*Z^2+rnorm(n)
      }else{
        y0<-0.5-0.5
        yx1<-0.5+2.5
        yx2<-0.5-0.5
        yz1<-0.5+4.5
        yz2<-0.5-0.5
        Y<-y0+yx1*X+yx2*X^2+yz1*Z+yz2*Z^2+rnorm(n)
      }
      Z_S<-Z[S==1]
      X_S<-X[S==1]
      Y_S<-Y[S==1]
    }else{
      X_Sdata<-rnorm(n)
      Z_Sdata<-rnorm(n,mean=-2)
      S_Sdata<-X_Sdata+Z_Sdata < -2
      if(linear==F){
        y0<-0.5-0.5
        yx1<-0.5-0.5
        yx2<-0.5+2.5
        yz1<-0.5+4.5
        yz2<-0.5-0.5
        Y_Sdata<-y0+yx1*X_Sdata+yx2*X_Sdata^2+yz1*Z_Sdata+yz2*Z_Sdata^2+rnorm(n)
      }else{
        y0<-0.5-0.5
        yx1<-0.5+2.5
        yx2<-0.5-0.5
        yz1<-0.5+4.5
        yz2<-0.5-0.5
        Y_Sdata<-y0+yx1*X_Sdata+yx2*X_Sdata^2+yz1*Z_Sdata+yz2*Z_Sdata^2+rnorm(n)
      }
      Z_S<-Z_Sdata[S_Sdata==1]
      X_S<-X_Sdata[S_Sdata==1]
      Y_S<-Y_Sdata[S_Sdata==1]
      
      X_Ddata<-rnorm(n)
      Z_Ddata<-rnorm(n,mean=-2)
      
      X<-X_Ddata
      Z<-Z_Ddata}
    
    #save realisations from all simulation runs
    if(empty==F){
      X_complete[,sim]<-X
      Z_complete[,sim]<-Z
      S_complete[,sim]<-S
      Y_complete[,sim]<-Y
    }else{ X_Sdata_complete[,sim] <- X_Sdata
    Z_Sdata_complete[,sim] <- Z_Sdata
    S_Sdata_complete[,sim] <- S_Sdata
    Y_Sdata_complete[,sim] <- Y_Sdata
    
    X_Ddata_complete[,sim] <- X_Ddata
    Z_Ddata_complete[,sim] <- Z_Ddata}
    
    
    #first step ridge 
    lambda_seq <- 10^seq(2, -2, by = -.1)
    ridge_cv <- cv.glmnet(cbind(X_S,I(X_S^2),Z_S,I(Z_S^2)), Y_S, alpha = 0, lambda = lambda_seq)
    best_lambda <- ridge_cv$lambda.min
    best_ridge <- glmnet(cbind(X_S,I(X_S^2),Z_S,I(Z_S^2)), Y_S, alpha = 0, lambda = best_lambda)
    betahat_first_ridge<-as.vector(coef(best_ridge))
    
    #first step OLS
    lm_fit <- lm(Y_S~X_S+I(X_S^2)+Z_S+I(Z_S^2))
    betahat_first<-coef(lm_fit)
    
    #calculate target variable for second step of RR
    Y_second_ridge<-(cbind(rep(1,n),X,X^2,Z,Z^2))%*%as.vector(betahat_first_ridge)
    Y_second<-(cbind(rep(1,n),X,X^2,Z,Z^2))%*%as.vector(betahat_first)
    
    #final RR estimates
    beta_end1_ridge<-lm(Y_second_ridge~1+X+I(X^2))$coef
    beta_end1<-lm(Y_second~1+X+I(X^2))$coef
    
    #naive estimation only based on S=1
    naive_beta<-lm(Y_S~1+X_S+I(X_S^2))$coef
    
    
    #E[Y|do(X)] = E[Y|X]
    original<-function(x)
    {y0+yx1*x+yx2*x^2+yz1*(-2)+yz2*5}
    
    #RR ridge
    est1_ridge<-function(x)
    {beta_end1_ridge[1]+beta_end1_ridge[2]*x+beta_end1_ridge[3]*x^2}
    
    #RR OLS
    est1<-function(x)
    {beta_end1[1]+beta_end1[2]*x+beta_end1[3]*x^2}
    
    #naive
    naive<-function(x)
    {naive_beta[1]+naive_beta[2]*x+naive_beta[3]*x^2}
    
    #TSR ridge
    est3_ridge<-function(x)
    {return(betahat_first_ridge[1]+betahat_first_ridge[2]*x+betahat_first_ridge[3]*x^2+
              betahat_first_ridge[4]*mean(Z)+
              betahat_first_ridge[5]*mean(Z^2))}
    
    #TSR OLS
    est3<-function(x)
    {return(betahat_first[1]+betahat_first[2]*x+betahat_first[3]*x^2+
              betahat_first[4]*mean(Z)+
              betahat_first[5]*mean(Z^2))}
    
    #save coefficient vectors from all simulation runs
    betahat_first_complete[,sim] <- betahat_first
    betahat_first_ridge_complete[,sim] <- betahat_first_ridge
    beta_end1_complete[,sim] <- beta_end1
    beta_end1_ridge_complete[,sim] <- beta_end1_ridge
    naive_beta_complete[,sim] <- naive_beta
    
    #test data  
    X_test<-as.matrix(rnorm(n=n))
    Z_test<-rnorm(n,mean=-2)
    S_test<-X_test+Z_test< -2
    if(linear==F){
      Y_test<-y0+yx1*X_test+yx2*X_test^2+yz1*Z_test+yz2*Z_test^2
    }else{
      Y_test<-y0+yx1*X_test+yx2*X_test^2+yz1*Z_test+yz2*Z_test^2
    }
    X_test_S<-as.matrix(X_test[S_test==1,])
    Y_test_S<-Y_test[S_test==1]
    
    #MSE (o: E[Y|do(X)] = E[Y|X]) 
    
    #in D
    o_test_naive[sim]<-mean((apply(X_test,1,original)-apply(X_test,1,naive))^2)
    o_test_est1[sim]<-mean((apply(X_test,1,original)-apply(X_test,1,est1))^2)
    o_test_est1_ridge[sim]<-mean((apply(X_test,1,original)-apply(X_test,1,est1_ridge))^2)
    o_test_est3[sim]<-mean((apply(X_test,1,original)-apply(X_test,1,est3))^2)
    o_test_est3_ridge[sim]<-mean((apply(X_test,1,original)-apply(X_test,1,est3_ridge))^2)
    
    #in S
    o_test_naive_S[sim]<-mean((apply(X_test_S,1,original)-apply(X_test_S,1,naive))^2)
    o_test_est1_S[sim]<-mean((apply(X_test_S,1,original)-apply(X_test_S,1,est1))^2)
    o_test_est1_ridge_S[sim]<-mean((apply(X_test_S,1,original)-apply(X_test_S,1,est1_ridge))^2)
    o_test_est3_S[sim]<-mean((apply(X_test_S,1,original)-apply(X_test_S,1,est3))^2)
    o_test_est3_ridge_S[sim]<-mean((apply(X_test_S,1,original)-apply(X_test_S,1,est3_ridge))^2)
    
    #Bias (o: E[Y|do(X)] = E[Y|X]) 
    
    #in D
    bias_o_test_naive[sim]<-mean((apply(X_test,1,original)-apply(X_test,1,naive)))
    bias_o_test_est1[sim]<-mean((apply(X_test,1,original)-apply(X_test,1,est1)))
    bias_o_test_est1_ridge[sim]<-mean((apply(X_test,1,original)-apply(X_test,1,est1_ridge)))
    bias_o_test_est3[sim]<-mean((apply(X_test,1,original)-apply(X_test,1,est3)))
    bias_o_test_est3_ridge[sim]<-mean((apply(X_test,1,original)-apply(X_test,1,est3_ridge)))
    
    #in S
    bias_o_test_naive_S[sim]<-mean((apply(X_test_S,1,original)-apply(X_test_S,1,naive)))
    bias_o_test_est1_S[sim]<-mean((apply(X_test_S,1,original)-apply(X_test_S,1,est1)))
    bias_o_test_est1_ridge_S[sim]<-mean((apply(X_test_S,1,original)-apply(X_test_S,1,est1_ridge)))
    bias_o_test_est3_S[sim]<-mean((apply(X_test_S,1,original)-apply(X_test_S,1,est3)))
    bias_o_test_est3_ridge_S[sim]<-mean((apply(X_test_S,1,original)-apply(X_test_S,1,est3_ridge)))
  }
  
  y0_complete[sim]<-y0
  yx1_complete[sim]<-yx1
  yx2_complete[sim]<-yx2
  yz1_complete[sim]<-yz1
  yz2_complete[sim]<-yz2
  
  #return the relevant results  
  if(empty==F){
    return(list(SIM = SIM, n = n, empty = empty, linear = linear,
                X_complete = X_complete, Z_complete = Z_complete, S_complete = S_complete, Y_complete = Y_complete,
                betahat_first_complete = betahat_first_complete, betahat_first_ridge_complete = betahat_first_ridge_complete,
                beta_end1_complete = beta_end1_complete, beta_end1_ridge_complete = beta_end1_ridge_complete,
                naive_beta_complete = naive_beta_complete,
                o_test_naive=o_test_naive, o_test_est1=o_test_est1, o_test_est1_ridge=o_test_est1_ridge, o_test_est3=o_test_est3, o_test_est3_ridge=o_test_est3_ridge,
                o_test_naive_S=o_test_naive_S, o_test_est1_S=o_test_est1_S, o_test_est1_ridge_S=o_test_est1_ridge_S, o_test_est3_S=o_test_est3_S, o_test_est3_ridge_S=o_test_est3_ridge_S,
                bias_o_test_naive=bias_o_test_naive, bias_o_test_est1=bias_o_test_est1, bias_o_test_est1_ridge=bias_o_test_est1_ridge, bias_o_test_est3=bias_o_test_est3, bias_o_test_est3_ridge=bias_o_test_est3_ridge,
                bias_o_test_naive_S=bias_o_test_naive_S, bias_o_test_est1_S=bias_o_test_est1_S, bias_o_test_est1_ridge_S=bias_o_test_est1_ridge_S, bias_o_test_est3_S=bias_o_test_est3_S, bias_o_test_est3_ridge_S=bias_o_test_est3_ridge_S,
                y0=y0,yx1=yx1,yx2=yx2,yz1=yz1,yz2=yz2))
  }else{
    return(list(SIM = SIM, n = n, empty = empty, linear = linear,
                X_Sdata_complete = X_Sdata_complete, Z_Sdata_complete = Z_Sdata_complete, S_Sdata_complete = S_Sdata_complete, Y_Sdata_complete = Y_Sdata_complete,
                X_Ddata_complete = X_Ddata_complete, Z_Ddata_complete = Z_Ddata_complete,
                betahat_first_complete = betahat_first_complete, betahat_first_ridge_complete = betahat_first_ridge_complete,
                beta_end1_complete = beta_end1_complete, beta_end1_ridge_complete = beta_end1_ridge_complete,
                naive_beta_complete = naive_beta_complete,
                o_test_naive=o_test_naive, o_test_est1=o_test_est1, o_test_est1_ridge=o_test_est1_ridge, o_test_est3=o_test_est3, o_test_est3_ridge=o_test_est3_ridge,
                o_test_naive_S=o_test_naive_S, o_test_est1_S=o_test_est1_S, o_test_est1_ridge_S=o_test_est1_ridge_S, o_test_est3_S=o_test_est3_S, o_test_est3_ridge_S=o_test_est3_ridge_S,
                bias_o_test_naive=bias_o_test_naive, bias_o_test_est1=bias_o_test_est1, bias_o_test_est1_ridge=bias_o_test_est1_ridge, bias_o_test_est3=bias_o_test_est3, bias_o_test_est3_ridge=bias_o_test_est3_ridge,
                bias_o_test_naive_S=bias_o_test_naive_S, bias_o_test_est1_S=bias_o_test_est1_S, bias_o_test_est1_ridge_S=bias_o_test_est1_ridge_S, bias_o_test_est3_S=bias_o_test_est3_S, bias_o_test_est3_ridge_S=bias_o_test_est3_ridge_S,
                y0=y0,yx1=yx1,yx2=yx2,yz1=yz1,yz2=yz2))}
  
}



#####################################################################################################
#####################################################################################################
#####################################################################################################

result500_1 <- sim1(n=500,SIM=100,empty=F,linear=F)
result1000_1 <- sim1(n=1000,SIM=100,empty=F,linear=F)
result5000_1 <- sim1(n=5000,SIM=100,empty=F,linear=F)

empty_result500_1 <- sim1(n=500,SIM=100,empty=T,linear=F)
empty_result1000_1 <- sim1(n=1000,SIM=100,empty=T,linear=F)
empty_result5000_1 <- sim1(n=5000,SIM=100,empty=T,linear=F)

linear_result500_1 <- sim1(n=500,SIM=100,empty=F,linear=T)
linear_result1000_1 <- sim1(n=1000,SIM=100,empty=F,linear=T)
linear_result5000_1 <- sim1(n=5000,SIM=100,empty=F,linear=T)

empty_linear_result500_1 <- sim1(n=500,SIM=100,empty=T,linear=T)
empty_linear_result1000_1 <- sim1(n=1000,SIM=100,empty=T,linear=T)
empty_linear_result5000_1 <- sim1(n=5000,SIM=100,empty=T,linear=T)



result500<-result500_1 
result1000<-result1000_1 
result5000<-result5000_1 

result500<-empty_result500_1 
result1000<-empty_result1000_1 
result5000<-empty_result5000_1 

result500<-linear_result500_1 
result1000<-linear_result1000_1 
result5000<-linear_result5000_1 

result500<-empty_linear_result500_1 
result1000<-empty_linear_result1000_1 
result5000<-empty_linear_result5000_1 



#Plot with curves for one specific case (n=500, quadratic, S \subset D)

result <- result500
SIM <- result$SIM
n <- result$n
empty <- result$empty
linear <- result$linear

X_est <- seq(from = -10,to = 10, by = 0.1)    
l_X_est <- length(X_est)
X_est <- as.matrix(X_est)

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
  
  if(empty==F){
    X <- as.vector(result$X_complete)
    Z <- as.vector(result$Z_complete)
    S <- as.vector(result$S_complete)
    Y <- as.vector(result$Y_complete)
    X_S<-X[S==1]
  }else{
    X <- as.vector(result$X_Ddata_complete)
    X_S <- as.vector(result$X_Sdata_complete[result$S_Sdata_complete==1])
    
    Z <- as.vector(result$Z_Ddata_complete)
    S <- as.vector(result$S_Ddata_complete)
    Y <- as.vector(result$Y_Ddata_complete)}
  
  #E[Y|do(X)] = E[Y|X]
  original<-function(x)
  {result$y0+result$yx1*x+result$yx2*x^2+result$yz1*(-2)+result$yz2*5}
  
  #RR ridge
  est1_ridge<-function(x)
  {beta_end1_ridge[1]+beta_end1_ridge[2]*x+beta_end1_ridge[3]*x^2}
  
  #RR OLS
  est1<-function(x)
  {beta_end1[1]+beta_end1[2]*x+beta_end1[3]*x^2}
  
  #naive
  naive<-function(x)
  {naive_beta[1]+naive_beta[2]*x+naive_beta[3]*x^2}
  
  #TSR ridge
  est3_ridge<-function(x)
  {return(betahat_first_ridge[1]+betahat_first_ridge[2]*x+betahat_first_ridge[3]*x^2+
            betahat_first_ridge[4]*mean(Z)+
            betahat_first_ridge[5]*mean(Z^2))}
  
  #TSR OLS
  est3<-function(x)
  {return(betahat_first[1]+betahat_first[2]*x+betahat_first[3]*x^2+
            betahat_first[4]*mean(Z)+
            betahat_first[5]*mean(Z^2))}
 
  
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

#calculate the means of the estimates
est1_mean<-apply(values_est_est1,1,mean)
est1_ridge_mean<-apply(values_est_est1_ridge,1,mean)
est3_mean<-apply(values_est_est3,1,mean)
est3_ridge_mean<-apply(values_est_est3_ridge,1,mean)
naive_mean<-apply(values_est_naive,1,mean)

green_rgb <- col2rgb("green")
darkgreen_rgb <- col2rgb("darkgreen")
middle_green_rgb <- (green_rgb + darkgreen_rgb) / 2
middle_green <- rgb(middle_green_rgb[1], middle_green_rgb[2], middle_green_rgb[3], maxColorValue = 255)

#quadratic

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
  geom_line(aes(x = X_est, y = apply(X_est,1,original), color = "ground truth"), linewidth = 1.3) +
  theme_minimal() +
  coord_cartesian(xlim = c(-4.6, 4.6), ylim = c(-30, 40)) + 
  geom_boxplot(aes(x = X, y = -27), width = 5, fill = "grey", color = "black", alpha = 0.7) +
  geom_boxplot(aes(x = X_S, y = -20), width = 5, fill = "grey", color = "black", alpha = 0.7) +
  scale_color_manual(values = c("ground truth" = "black", "naive"="orange", "RR"="blue", "TSR"=middle_green), name = NULL) +
  scale_fill_manual(values = c("naive" = "orange", "RR" = "lightblue", "TSR" = "lightgreen"), name = NULL) +
  theme(legend.position = "top")+
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
  geom_line(aes(x = X_est, y = apply(X_est,1,original), color = "ground truth"), linewidth = 1.3) +
  theme_minimal() +
  coord_cartesian(xlim = c(-4.6, 4.6), ylim = c(-30, 40)) + 
  geom_boxplot(aes(x = X, y = -27), width = 5, fill = "grey", color = "black", alpha = 0.7) +
  geom_boxplot(aes(x = X_S, y = -20), width = 5, fill = "grey", color = "black", alpha = 0.7) +
  scale_color_manual(values = c("ground truth" = "black", "naive"="orange", "RR (ridge)"="blue", "TSR (ridge)"=middle_green), name = NULL) +
  scale_fill_manual(values = c("naive" = "orange", "RR (ridge)" = "lightblue", "TSR (ridge)" = "lightgreen"), name = NULL) +
  theme(legend.position = "top")+
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



#linear


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
  geom_line(aes(x = X_est, y = apply(X_est,1,original), color = "ground truth"), linewidth = 1.3) +
  theme_minimal() +
  coord_cartesian(xlim = c(-4.6, 4.6), ylim = c(-50, 10)) + 
  geom_boxplot(aes(x = X, y = -47), width = 5, fill = "grey", color = "black", alpha = 0.7) +
  geom_boxplot(aes(x = X_S, y = -40), width = 5, fill = "grey", color = "black", alpha = 0.7) +
  scale_color_manual(values = c("ground truth" = "black", "naive"="orange", "RR"="blue", "TSR"=middle_green), name = NULL) +
  scale_fill_manual(values = c("naive" = "orange", "RR" = "lightblue", "TSR" = "lightgreen"), name = NULL) +
  theme(legend.position = "top")+
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
  geom_line(aes(x = X_est, y = apply(X_est,1,original), color = "ground truth"), linewidth = 1.3) +
  theme_minimal() +
  coord_cartesian(xlim = c(-4.6, 4.6), ylim = c(-50, 10)) + 
  geom_boxplot(aes(x = X, y = -47), width = 5, fill = "grey", color = "black", alpha = 0.7) +
  geom_boxplot(aes(x = X_S, y = -40), width = 5, fill = "grey", color = "black", alpha = 0.7) +
  scale_color_manual(values = c("ground truth" = "black", "naive"="orange", "RR (ridge)"="blue", "TSR (ridge)"=middle_green), name = NULL) +
  scale_fill_manual(values = c("naive" = "orange", "RR (ridge)" = "lightblue", "TSR (ridge)" = "lightgreen"), name = NULL) +
  theme(legend.position = "top")+
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




#sec 4.1 results -randomized coefficients

####################################
#######      quadratic      ########
####################################



#MSE mean
c(mean(result500_2$o_test_naive_S), mean(result500_2$o_test_est1_S), mean(result500_2$o_test_est1_ridge_S), mean(result500_2$o_test_est3_S), mean(result500_2$o_test_est3_ridge_S))
#[1] 13.12729337  0.18435783  0.18053287  0.07594283  0.07420850
c(mean(result1000_2$o_test_naive_S), mean(result1000_2$o_test_est1_S), mean(result1000_2$o_test_est1_ridge_S), mean(result1000_2$o_test_est3_S), mean(result1000_2$o_test_est3_ridge_S))
#[1] 12.85392091  0.09129182  0.09093798  0.03459136  0.03438119
c(mean(result5000_2$o_test_naive_S), mean(result5000_2$o_test_est1_S), mean(result5000_2$o_test_est1_ridge_S), mean(result5000_2$o_test_est3_S), mean(result5000_2$o_test_est3_ridge_S))
#[1] 13.271128505  0.018520669  0.020037311  0.006500773  0.007995492

c(mean(result500_2$o_test_naive), mean(result500_2$o_test_est1), mean(result500_2$o_test_est1_ridge), mean(result500_2$o_test_est3), mean(result500_2$o_test_est3_ridge))
#[1] 35.28125200  0.19697496  0.19423682  0.10021599  0.09947769
c(mean(result1000_2$o_test_naive), mean(result1000_2$o_test_est1), mean(result1000_2$o_test_est1_ridge), mean(result1000_2$o_test_est3), mean(result1000_2$o_test_est3_ridge))
#[1] 33.04541112  0.10298463  0.10143157  0.04855360  0.04883651
c(mean(result5000_2$o_test_naive), mean(result5000_2$o_test_est1), mean(result5000_2$o_test_est1_ridge), mean(result5000_2$o_test_est3), mean(result5000_2$o_test_est3_ridge))
#[1] 34.919597889  0.021852740  0.023927974  0.009162473  0.011438996


#MSE sd
c(sd(result500_2$o_test_naive_S), sd(result500_2$o_test_est1_S), sd(result500_2$o_test_est1_ridge_S), sd(result500_2$o_test_est3_S), sd(result500_2$o_test_est3_ridge_S))
#[1] 6.62121597 0.20781937 0.20325910 0.07841582 0.07691176
c(sd(result1000_2$o_test_naive_S), sd(result1000_2$o_test_est1_S), sd(result1000_2$o_test_est1_ridge_S), sd(result1000_2$o_test_est3_S), sd(result1000_2$o_test_est3_ridge_S))
#[1] 6.10119809 0.08475962 0.08523829 0.04206452 0.04176733
c(sd(result5000_2$o_test_naive_S), sd(result5000_2$o_test_est1_S), sd(result5000_2$o_test_est1_ridge_S), sd(result5000_2$o_test_est3_S), sd(result5000_2$o_test_est3_ridge_S))
#[1] 6.016373865 0.018511941 0.018688187 0.006655763 0.008053742

c(sd(result500_2$o_test_naive), sd(result500_2$o_test_est1), sd(result500_2$o_test_est1_ridge), sd(result500_2$o_test_est3), sd(result500_2$o_test_est3_ridge))
#[1] 19.4606750  0.1715773  0.1703883  0.0962441  0.1001531
c(sd(result1000_2$o_test_naive), sd(result1000_2$o_test_est1), sd(result1000_2$o_test_est1_ridge), sd(result1000_2$o_test_est3), sd(result1000_2$o_test_est3_ridge))
#[1] 16.31962070  0.08453729  0.08654949  0.05664407  0.05932058
c(sd(result5000_2$o_test_naive), sd(result5000_2$o_test_est1), sd(result5000_2$o_test_est1_ridge), sd(result5000_2$o_test_est3), sd(result5000_2$o_test_est3_ridge))
#[1] 17.294978771  0.020709634  0.020249562  0.007977598  0.010152760



#################
###   empty   ###
#################



#MSE mean
c(mean(empty_result500_2$o_test_naive_S), mean(empty_result500_2$o_test_est1_S), mean(empty_result500_2$o_test_est1_ridge_S), mean(empty_result500_2$o_test_est3_S), mean(empty_result500_2$o_test_est3_ridge_S))
#[1] 13.39166109  0.16455281  0.16300262  0.06199675  0.06163348
c(mean(empty_result1000_2$o_test_naive_S), mean(empty_result1000_2$o_test_est1_S), mean(empty_result1000_2$o_test_est1_ridge_S), mean(empty_result1000_2$o_test_est3_S), mean(empty_result1000_2$o_test_est3_ridge_S))
#[1] 12.76928671  0.08312488  0.08222381  0.03103097  0.03083180
c(mean(empty_result5000_2$o_test_naive_S), mean(empty_result5000_2$o_test_est1_S), mean(empty_result5000_2$o_test_est1_ridge_S), mean(empty_result5000_2$o_test_est3_S), mean(empty_result5000_2$o_test_est3_ridge_S))
#[1] 13.218647518  0.021143450  0.022092247  0.009017970  0.009994856

c(mean(empty_result500_2$o_test_naive), mean(empty_result500_2$o_test_est1), mean(empty_result500_2$o_test_est1_ridge), mean(empty_result500_2$o_test_est3), mean(empty_result500_2$o_test_est3_ridge))
#[1] 35.19287009  0.22509545  0.22371851  0.09390262  0.09489599
c(mean(empty_result1000_2$o_test_naive), mean(empty_result1000_2$o_test_est1), mean(empty_result1000_2$o_test_est1_ridge), mean(empty_result1000_2$o_test_est3), mean(empty_result1000_2$o_test_est3_ridge))
#[1] 33.41623021  0.09603687  0.09568151  0.04270291  0.04329694
c(mean(empty_result5000_2$o_test_naive), mean(empty_result5000_2$o_test_est1), mean(empty_result5000_2$o_test_est1_ridge), mean(empty_result5000_2$o_test_est3), mean(empty_result5000_2$o_test_est3_ridge))
#[1] 35.06929702  0.02430621  0.02554498  0.01276720  0.01425024

#MSE sd
c(sd(empty_result500_2$o_test_naive_S), sd(empty_result500_2$o_test_est1_S), sd(empty_result500_2$o_test_est1_ridge_S), sd(empty_result500_2$o_test_est3_S), sd(empty_result500_2$o_test_est3_ridge_S))
#[1] 6.79828963 0.18063415 0.18049279 0.07330028 0.07490283
c(sd(empty_result1000_2$o_test_naive_S), sd(empty_result1000_2$o_test_est1_S), sd(empty_result1000_2$o_test_est1_ridge_S), sd(empty_result1000_2$o_test_est3_S), sd(empty_result1000_2$o_test_est3_ridge_S))
#[1] 5.84977863 0.08909894 0.08925492 0.04693764 0.04694088
c(sd(empty_result5000_2$o_test_naive_S), sd(empty_result5000_2$o_test_est1_S), sd(empty_result5000_2$o_test_est1_ridge_S), sd(empty_result5000_2$o_test_est3_S), sd(empty_result5000_2$o_test_est3_ridge_S))
#[1] 5.97278581 0.02154902 0.02121897 0.01199304 0.01234643

c(sd(empty_result500_2$o_test_naive), sd(empty_result500_2$o_test_est1), sd(empty_result500_2$o_test_est1_ridge), sd(empty_result500_2$o_test_est3), sd(empty_result500_2$o_test_est3_ridge))
#[1] 19.8857301  0.2343205  0.2349753  0.1069215  0.1114006
c(sd(empty_result1000_2$o_test_naive), sd(empty_result1000_2$o_test_est1), sd(empty_result1000_2$o_test_est1_ridge), sd(empty_result1000_2$o_test_est3), sd(empty_result1000_2$o_test_est3_ridge))
#[1] 16.65001887  0.08800525  0.09031048  0.05180145  0.05636804
c(sd(empty_result5000_2$o_test_naive), sd(empty_result5000_2$o_test_est1), sd(empty_result5000_2$o_test_est1_ridge), sd(empty_result5000_2$o_test_est3), sd(empty_result5000_2$o_test_est3_ridge))
#[1] 17.51899866  0.02436812  0.02532257  0.01528951  0.01701702




####################################
#######       linear        ########
####################################



#MSE mean
c(mean(linear_result500_2$o_test_naive_S), mean(linear_result500_2$o_test_est1_S), mean(linear_result500_2$o_test_est1_ridge_S), mean(linear_result500_2$o_test_est3_S), mean(linear_result500_2$o_test_est3_ridge_S))
#[1] 13.12729337  0.18435783  0.18121136  0.07594283  0.07553143
c(mean(linear_result1000_2$o_test_naive_S), mean(linear_result1000_2$o_test_est1_S), mean(linear_result1000_2$o_test_est1_ridge_S), mean(linear_result1000_2$o_test_est3_S), mean(linear_result1000_2$o_test_est3_ridge_S))
#[1] 12.85392091  0.09129182  0.09217859  0.03459136  0.03606434
c(mean(linear_result5000_2$o_test_naive_S), mean(linear_result5000_2$o_test_est1_S), mean(linear_result5000_2$o_test_est1_ridge_S), mean(linear_result5000_2$o_test_est3_S), mean(linear_result5000_2$o_test_est3_ridge_S))
#[1] 13.271128505  0.018520669  0.022594427  0.006500773  0.010677841

c(mean(linear_result500_2$o_test_naive), mean(linear_result500_2$o_test_est1), mean(linear_result500_2$o_test_est1_ridge), mean(linear_result500_2$o_test_est3), mean(linear_result500_2$o_test_est3_ridge))
#[1] 35.2812520  0.1969750  0.1965480  0.1002160  0.1017934
c(mean(linear_result1000_2$o_test_naive), mean(linear_result1000_2$o_test_est1), mean(linear_result1000_2$o_test_est1_ridge), mean(linear_result1000_2$o_test_est3), mean(linear_result1000_2$o_test_est3_ridge))
#[1] 33.04541112  0.10298463  0.10298769  0.04855360  0.05148678
c(mean(linear_result5000_2$o_test_naive), mean(linear_result5000_2$o_test_est1), mean(linear_result5000_2$o_test_est1_ridge), mean(linear_result5000_2$o_test_est3), mean(linear_result5000_2$o_test_est3_ridge))
#[1] 34.919597889  0.021852740  0.027660256  0.009162473  0.015384191

#MSE sd
c(sd(linear_result500_2$o_test_naive_S), sd(linear_result500_2$o_test_est1_S), sd(linear_result500_2$o_test_est1_ridge_S), sd(linear_result500_2$o_test_est3_S), sd(linear_result500_2$o_test_est3_ridge_S))
#[1] 6.62121597 0.20781937 0.20200752 0.07841582 0.07754779
c(sd(linear_result1000_2$o_test_naive_S), sd(linear_result1000_2$o_test_est1_S), sd(linear_result1000_2$o_test_est1_ridge_S), sd(linear_result1000_2$o_test_est3_S), sd(linear_result1000_2$o_test_est3_ridge_S))
#[1] 6.10119809 0.08475962 0.08718952 0.04206452 0.04395317
c(sd(linear_result5000_2$o_test_naive_S), sd(linear_result5000_2$o_test_est1_S), sd(linear_result5000_2$o_test_est1_ridge_S), sd(linear_result5000_2$o_test_est3_S), sd(linear_result5000_2$o_test_est3_ridge_S))
#[1] 6.016373865 0.018511941 0.019026565 0.006655763 0.010483599

c(sd(linear_result500_2$o_test_naive), sd(linear_result500_2$o_test_est1), sd(linear_result500_2$o_test_est1_ridge), sd(linear_result500_2$o_test_est3), sd(linear_result500_2$o_test_est3_ridge))
#[1] 19.4606750  0.1715773  0.1710995  0.0962441  0.1035901
c(sd(linear_result1000_2$o_test_naive), sd(linear_result1000_2$o_test_est1), sd(linear_result1000_2$o_test_est1_ridge), sd(linear_result1000_2$o_test_est3), sd(linear_result1000_2$o_test_est3_ridge))
#[1] 16.31962070  0.08453729  0.09017366  0.05664407  0.06382551
c(sd(linear_result5000_2$o_test_naive), sd(linear_result5000_2$o_test_est1), sd(linear_result5000_2$o_test_est1_ridge), sd(linear_result5000_2$o_test_est3), sd(linear_result5000_2$o_test_est3_ridge))
#[1] 17.294978771  0.020709634  0.021142322  0.007977598  0.013953557


#################
###   empty   ###
#################



#MSE mean
c(mean(empty_linear_result500_2$o_test_naive_S), mean(empty_linear_result500_2$o_test_est1_S), mean(empty_linear_result500_2$o_test_est1_ridge_S), mean(empty_linear_result500_2$o_test_est3_S), mean(empty_linear_result500_2$o_test_est3_ridge_S))
#[1] 13.39166109  0.16455281  0.16253843  0.06199675  0.06274075
c(mean(empty_linear_result1000_2$o_test_naive_S), mean(empty_linear_result1000_2$o_test_est1_S), mean(empty_linear_result1000_2$o_test_est1_ridge_S), mean(empty_linear_result1000_2$o_test_est3_S), mean(empty_linear_result1000_2$o_test_est3_ridge_S))
#[1] 12.76928671  0.08312488  0.08345375  0.03103097  0.03270526
c(mean(empty_linear_result5000_2$o_test_naive_S), mean(empty_linear_result5000_2$o_test_est1_S), mean(empty_linear_result5000_2$o_test_est1_ridge_S), mean(empty_linear_result5000_2$o_test_est3_S), mean(empty_linear_result5000_2$o_test_est3_ridge_S))
#[1] 13.21864752  0.02114345  0.02429830  0.00901797  0.01209800

c(mean(empty_linear_result500_2$o_test_naive), mean(empty_linear_result500_2$o_test_est1), mean(empty_linear_result500_2$o_test_est1_ridge), mean(empty_linear_result500_2$o_test_est3), mean(empty_linear_result500_2$o_test_est3_ridge))
#[1] 35.19287009  0.22509545  0.22365549  0.09390262  0.09680743
c(mean(empty_linear_result1000_2$o_test_naive), mean(empty_linear_result1000_2$o_test_est1), mean(empty_linear_result1000_2$o_test_est1_ridge), mean(empty_linear_result1000_2$o_test_est3), mean(empty_linear_result1000_2$o_test_est3_ridge))
#[1] 33.41623021  0.09603687  0.09768017  0.04270291  0.04605673
c(mean(empty_linear_result5000_2$o_test_naive), mean(empty_linear_result5000_2$o_test_est1), mean(empty_linear_result5000_2$o_test_est1_ridge), mean(empty_linear_result5000_2$o_test_est3), mean(empty_linear_result5000_2$o_test_est3_ridge))
#[1] 35.06929702  0.02430621  0.02873780  0.01276720  0.01745335

#MSE sd
c(sd(empty_linear_result500_2$o_test_naive_S), sd(empty_linear_result500_2$o_test_est1_S), sd(empty_linear_result500_2$o_test_est1_ridge_S), sd(empty_linear_result500_2$o_test_est3_S), sd(empty_linear_result500_2$o_test_est3_ridge_S))
#[1] 6.79828963 0.18063415 0.17981544 0.07330028 0.07528662
c(sd(empty_linear_result1000_2$o_test_naive_S), sd(empty_linear_result1000_2$o_test_est1_S), sd(empty_linear_result1000_2$o_test_est1_ridge_S), sd(empty_linear_result1000_2$o_test_est3_S), sd(empty_linear_result1000_2$o_test_est3_ridge_S))
#[1] 5.84977863 0.08909894 0.09037558 0.04693764 0.04901321
c(sd(empty_linear_result5000_2$o_test_naive_S), sd(empty_linear_result5000_2$o_test_est1_S), sd(empty_linear_result5000_2$o_test_est1_ridge_S), sd(empty_linear_result5000_2$o_test_est3_S), sd(empty_linear_result5000_2$o_test_est3_ridge_S))
#[1] 5.97278581 0.02154902 0.02158416 0.01199304 0.01345607

c(sd(empty_linear_result500_2$o_test_naive), sd(empty_linear_result500_2$o_test_est1), sd(empty_linear_result500_2$o_test_est1_ridge), sd(empty_linear_result500_2$o_test_est3), sd(empty_linear_result500_2$o_test_est3_ridge))
#[1] 19.8857301  0.2343205  0.2346895  0.1069215  0.1135098
c(sd(empty_linear_result1000_2$o_test_naive), sd(empty_linear_result1000_2$o_test_est1), sd(empty_linear_result1000_2$o_test_est1_ridge), sd(empty_linear_result1000_2$o_test_est3), sd(empty_linear_result1000_2$o_test_est3_ridge))
#[1] 16.65001887  0.08800525  0.09274316  0.05180145  0.06083609
c(sd(empty_linear_result5000_2$o_test_naive), sd(empty_linear_result5000_2$o_test_est1), sd(empty_linear_result5000_2$o_test_est1_ridge), sd(empty_linear_result5000_2$o_test_est3), sd(empty_linear_result5000_2$o_test_est3_ridge))
#[1] 17.51899866  0.02436812  0.02684439  0.01528951  0.01971154



