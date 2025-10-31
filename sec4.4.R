#Robustness to misspecification

#load packages
library(glmnet)
library(ggplot2)

#n = sample size
#SIM = number of simulation runs
#robust takes values (0,1,2,3,4,5), values explained below

sim1<-function(n,SIM,robust)
{
  #create empty matrices to save results

  X_Sdata_complete <- matrix(ncol=SIM,nrow=n)
  Z_Sdata_complete <- matrix(ncol=SIM,nrow=n)
  U_Sdata_complete <- matrix(ncol=SIM,nrow=n)
  S_Sdata_complete <- matrix(ncol=SIM,nrow=n)
  Y_Sdata_complete <- matrix(ncol=SIM,nrow=n)
  
  X_Ddata_complete <- matrix(ncol=SIM,nrow=n)
  Z_Ddata_complete <- matrix(ncol=SIM,nrow=n)
  U_Ddata_complete <- matrix(ncol=SIM,nrow=n)
  
  betahat_first_complete <- matrix(ncol=SIM,nrow=5)
  betahat_first_ridge_complete <- matrix(ncol=SIM,nrow=5)
  beta_end1_complete <- matrix(ncol=SIM, nrow=3)
  beta_end1_ridge_complete <- matrix(ncol=SIM, nrow=3)
  naive_beta_complete <- matrix(ncol=SIM, nrow=3)
  
  for(sim in 1:SIM)
  {
    print(sim)
    
    #reproducibility
    set.seed(sim)
    
    #data generating processes
      
      if(robust==0){            #without U
        
        Z_Sdata<-rnorm(n,mean=-2)
        X_Sdata<-rnorm(n)+2*rnorm(n,mean=-2)
        S_Sdata<-X_Sdata+Z_Sdata< -5
        Y_Sdata<-0.2*X_Sdata^2+5*Z_Sdata+rnorm(n)
        
        Z_S<-Z_Sdata[S_Sdata==1]
        X_S<-X_Sdata[S_Sdata==1]
        Y_S<-Y_Sdata[S_Sdata==1]
        
        Z_Ddata<-rnorm(n,mean=-2)
        X_Ddata<-rnorm(n)+2*rnorm(n,mean=-2)
        
        X<-X_Ddata
        Z<-Z_Ddata
        }
    
      if(robust==1){           #missing covariate
        
        Z_Sdata<-rnorm(n,mean=-2)
        X_Sdata<-rnorm(n)+2*rnorm(n,mean=-2)
        U_Sdata<-rnorm(n,mean=3)
        S_Sdata<-X_Sdata+rnorm(n,mean=-2)< -5
        Y_Sdata<-0.2*X_Sdata^2+5*Z_Sdata+4*U_Sdata+rnorm(n)
        
        Z_S<-Z_Sdata[S_Sdata==1]
        X_S<-X_Sdata[S_Sdata==1]
        Y_S<-Y_Sdata[S_Sdata==1]
        U_S<-U_Sdata[S_Sdata==1]
        
        Z_Ddata<-rnorm(n,mean=-2)
        X_Ddata<-rnorm(n)+2*rnorm(n,mean=-2)
        U_Ddata<-rnorm(n,mean=3)
        
        X<-X_Ddata
        Z<-Z_Ddata
        U<-U_Ddata
        
        U_Sdata_complete[,sim] <- U_Sdata
        U_Ddata_complete[,sim] <- U_Ddata
        }
        
      if(robust==2){           #U as confounder (type 1)
        
        Z_Sdata<-rnorm(n,mean=-2)
        X_Sdata<-rnorm(n)+2*rnorm(n,mean=-2)
        U_Sdata<-rnorm(n,mean=3)
        S_Sdata<-X_Sdata+Z_Sdata+0.1*U_Sdata< -5
        Y_Sdata<-0.2*X_Sdata^2+5*Z_Sdata+4*U_Sdata+rnorm(n)
        
        Z_S<-Z_Sdata[S_Sdata==1]
        X_S<-X_Sdata[S_Sdata==1]
        Y_S<-Y_Sdata[S_Sdata==1]
        U_S<-U_Sdata[S_Sdata==1]
        
        Z_Ddata<-rnorm(n,mean=-2)
        X_Ddata<-rnorm(n)+2*rnorm(n,mean=-2)
        U_Ddata<-rnorm(n,mean=3)
        
        X<-X_Ddata
        Z<-Z_Ddata
        U<-U_Ddata
        
        U_Sdata_complete[,sim] <- U_Sdata
        U_Ddata_complete[,sim] <- U_Ddata
        }
  
      if(robust==3){           #U as cause of S (type 1)
      
        Z_Sdata<-rnorm(n,mean=-2)
        X_Sdata<-rnorm(n)+2*rnorm(n,mean=-2)
        U_Sdata<-rnorm(n,mean=3)
        S_Sdata<-X_Sdata+Z_Sdata+0.1*U_Sdata< -5
        Y_Sdata<-0.2*X_Sdata^2+5*Z_Sdata+rnorm(n)
      
        Z_S<-Z_Sdata[S_Sdata==1]
        X_S<-X_Sdata[S_Sdata==1]
        Y_S<-Y_Sdata[S_Sdata==1]
        U_S<-U_Sdata[S_Sdata==1]
      
        Z_Ddata<-rnorm(n,mean=-2)
        X_Ddata<-rnorm(n)+2*rnorm(n,mean=-2)
        U_Ddata<-rnorm(n,mean=3)
      
        X<-X_Ddata
        Z<-Z_Ddata
        U<-U_Ddata
        
        U_Sdata_complete[,sim] <- U_Sdata
        U_Ddata_complete[,sim] <- U_Ddata
        }
    
    if(robust==4){           #U as confounder (type 2)
      
      Z_Sdata<-rnorm(n,mean=-2)
      X_Sdata<-rnorm(n)+2*rnorm(n,mean=-2)
      U_Sdata<-rnorm(n,mean=3)
      S_Sdata<-(X_Sdata+Z_Sdata < -5) & (U_Sdata < 1.5)
      Y_Sdata<-0.2*X_Sdata^2+5*Z_Sdata+4*U_Sdata+rnorm(n)
      
      Z_S<-Z_Sdata[S_Sdata==1]
      X_S<-X_Sdata[S_Sdata==1]
      Y_S<-Y_Sdata[S_Sdata==1]
      U_S<-U_Sdata[S_Sdata==1]
      
      Z_Ddata<-rnorm(n,mean=-2)
      X_Ddata<-rnorm(n)+2*rnorm(n,mean=-2)
      U_Ddata<-rnorm(n,mean=3)
      
      X<-X_Ddata
      Z<-Z_Ddata
      U<-U_Ddata
      
      U_Sdata_complete[,sim] <- U_Sdata
      U_Ddata_complete[,sim] <- U_Ddata
      }
      
    if(robust==5){           #U as cause of S (type 2)
      
      Z_Sdata<-rnorm(n,mean=-2)
      X_Sdata<-rnorm(n)+2*rnorm(n,mean=-2)
      U_Sdata<-rnorm(n,mean=3)
      S_Sdata<-(X_Sdata+Z_Sdata < -5) & (U_Sdata < 1.5)
      Y_Sdata<-0.2*X_Sdata^2+5*Z_Sdata+rnorm(n)
      
      Z_S<-Z_Sdata[S_Sdata==1]
      X_S<-X_Sdata[S_Sdata==1]
      Y_S<-Y_Sdata[S_Sdata==1]
      U_S<-U_Sdata[S_Sdata==1]
      
      Z_Ddata<-rnorm(n,mean=-2)
      X_Ddata<-rnorm(n)+2*rnorm(n,mean=-2)
      U_Ddata<-rnorm(n,mean=3)
      
      X<-X_Ddata
      Z<-Z_Ddata
      U<-U_Ddata
      
      U_Sdata_complete[,sim] <- U_Sdata
      U_Ddata_complete[,sim] <- U_Ddata
      }
    
    #save realisations from all simulation runs
    X_Sdata_complete[,sim] <- X_Sdata
    Z_Sdata_complete[,sim] <- Z_Sdata
    S_Sdata_complete[,sim] <- S_Sdata
    Y_Sdata_complete[,sim] <- Y_Sdata

    X_Ddata_complete[,sim] <- X_Ddata
    Z_Ddata_complete[,sim] <- Z_Ddata

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
    
    #RR ridge
    est1_ridge<-function(x)
    {
      beta_end1_ridge[1]+beta_end1_ridge[2]*x+beta_end1_ridge[3]*x^2
    }
    
    #RR OLS
    est1<-function(x)
    {
      beta_end1[1]+beta_end1[2]*x+beta_end1[3]*x^2
    }
    
    #naive
    naive<-function(x)
    {
      naive_beta[1]+naive_beta[2]*x+naive_beta[3]*x^2
    }
    
    #TSR ridge
    est3_ridge<-function(x)
    {
      return(betahat_first_ridge[1]+betahat_first_ridge[2]*x+betahat_first_ridge[3]*x^2+
             betahat_first_ridge[4]*mean(Z)+
             betahat_first_ridge[5]*mean(Z^2)
      )
    }
    
    #TSR OLS
    est3<-function(x)
    {
      return(betahat_first[1]+betahat_first[2]*x+betahat_first[3]*x^2+
             betahat_first[4]*mean(Z)+
             betahat_first[5]*mean(Z^2)
      )
    }
    
    #save coefficient vectors from all simulation runs
    betahat_first_complete[,sim] <- betahat_first
    betahat_first_ridge_complete[,sim] <- betahat_first_ridge
    beta_end1_complete[,sim] <- beta_end1
    beta_end1_ridge_complete[,sim] <- beta_end1_ridge
    naive_beta_complete[,sim] <- naive_beta
    
  }
  
  
  #return the relevant results
  return(list(SIM = SIM, n = n , robust = robust,
              X_Sdata_complete = X_Sdata_complete, Z_Sdata_complete = Z_Sdata_complete, U_Sdata_complete = U_Sdata_complete, S_Sdata_complete = S_Sdata_complete, Y_Sdata_complete = Y_Sdata_complete,
              X_Ddata_complete = X_Ddata_complete, Z_Ddata_complete = Z_Ddata_complete, U_Ddata_complete = U_Ddata_complete,
              betahat_first_complete = betahat_first_complete, betahat_first_ridge_complete = betahat_first_ridge_complete,
              beta_end1_complete = beta_end1_complete, beta_end1_ridge_complete = beta_end1_ridge_complete,
              naive_beta_complete = naive_beta_complete))
  }


#####################################################################################################
#####################################################################################################
#####################################################################################################

result1000_robust0 <- sim1(n=1000,SIM=100,robust=0)
result1000_robust1 <- sim1(n=1000,SIM=100,robust=1)
result1000_robust2 <- sim1(n=1000,SIM=100,robust=2)
result1000_robust3 <- sim1(n=1000,SIM=100,robust=3)
result1000_robust4 <- sim1(n=1000,SIM=100,robust=4)
result1000_robust5 <- sim1(n=1000,SIM=100,robust=5)

result5000_robust0 <- sim1(n=5000,SIM=100,robust=0)
result5000_robust1 <- sim1(n=5000,SIM=100,robust=1)
result5000_robust2 <- sim1(n=5000,SIM=100,robust=2)
result5000_robust3 <- sim1(n=5000,SIM=100,robust=3)
result5000_robust4 <- sim1(n=5000,SIM=100,robust=4)
result5000_robust5 <- sim1(n=5000,SIM=100,robust=5)

#save(result1000_robust0,result1000_robust1,result1000_robust2,result1000_robust3,result1000_robust4,result1000_robust5,
#     result5000_robust0,result5000_robust1,result5000_robust2,result5000_robust3,result5000_robust4,result5000_robust5,
#     file="robust_5cases.Rdata") 

load("robust_5cases.Rdata")


#plots with 95% areas

result <- result5000_robust5
set.seed(1)

SIM <- result$SIM
n <- result$n
robust <- result$robust

X_est <- seq(from = -15,to = 10, by = 0.1)    
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
  
  Z <- result$Z_Ddata_complete[,sim]
  
  #E[Y|do(X)] 
    
  if(robust==0 | robust==3 | robust==5){
    original <- function(x)
    {
      0.2*x^2-10
    }}
  if(robust==1 | robust==2 | robust==4){
    original <- function(x)
    {
      0.3*x^2-10+12
    }}
  est1_ridge <- function(x)
  {
    beta_end1_ridge[1]+beta_end1_ridge[2]*x+beta_end1_ridge[3]*x^2
  }
  
  est1 <- function(x)
  {
    beta_end1[1]+beta_end1[2]*x+beta_end1[3]*x^2
  }
  
  naive <- function(x)
  {
    naive_beta[1]+naive_beta[2]*x+naive_beta[3]*x^2
  }
  
  est3_ridge <- function(x)
  {
    return(betahat_first_ridge[1]+betahat_first_ridge[2]*x+betahat_first_ridge[3]*x^2+
           betahat_first_ridge[4]*mean(Z)+
           betahat_first_ridge[5]*mean(Z^2)
    )
  }
  
  est3 <- function(x)
  {
    return(betahat_first[1]+betahat_first[2]*x+betahat_first[3]*x^2+
           betahat_first[4]*mean(Z)+
           betahat_first[5]*mean(Z^2)
    )
  }
  
  values_est_naive[,sim] <- apply(X_est,1,naive)
  values_est_est1[,sim] <- apply(X_est,1,est1)
  values_est_est1_ridge[,sim] <- apply(X_est,1,est1_ridge)
  values_est_est3[,sim] <- apply(X_est,1,est3)
  values_est_est3_ridge[,sim] <- apply(X_est,1,est3_ridge)
  
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

#calculate means of estimations
est1_mean <- apply(values_est_est1,1,mean)
est1_ridge_mean <- apply(values_est_est1_ridge,1,mean)
est3_mean <- apply(values_est_est3,1,mean)
est3_ridge_mean <- apply(values_est_est3_ridge,1,mean)
naive_mean <- apply(values_est_naive,1,mean)

#define middle_green
green_rgb <- col2rgb("green")
darkgreen_rgb <- col2rgb("darkgreen")
middle_green_rgb <- (green_rgb + darkgreen_rgb) / 2
middle_green <- rgb(middle_green_rgb[1], middle_green_rgb[2], middle_green_rgb[3], maxColorValue = 255)

#for the boxplots in the figure
X <- as.vector(result$X_Ddata_complete)
X_S <- as.vector(result$X_Sdata_complete)[as.vector(result$S_Sdata_complete)==1]

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
  geom_line(aes(x = X_est, y = apply(X_est,1,original), color = "E[Y|do(X=x)]"), linewidth = 1.3) +
  theme_minimal() +
  coord_cartesian(xlim = c(-15, 7.5), ylim = c(-28, 30)) + 
  geom_boxplot(aes(x = X, y = -27), width = 5, fill = "grey", color = "black", alpha = 0.7) +
  geom_boxplot(aes(x = X_S, y = -20), width = 5, fill = "grey", color = "black", alpha = 0.7) +
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
    legend.position = "right",
    legend.spacing = unit(1.5, "cm"),
    legend.key.height = unit(1.3, "cm"),  
    legend.key.width = unit(1.3, "cm")
  )+
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(order = 2)
  )



