
#Section 4.3

#load packages
library(glmnet)
library(ggplot2)

#n = sample size
#SIM = simulation runs
#empty = indicates if ((S cap D)=emptyset)

sim1<-function(n,SIM,empty){
  
  #create empty matrices to save results 
  
  if(empty==F){
    X_complete <- matrix(ncol=SIM,nrow=n)
    Zm_complete <- matrix(ncol=SIM,nrow=n)
    Zp_complete <- matrix(ncol=SIM,nrow=n)
    S_complete <- matrix(ncol=SIM,nrow=n)
    Y_complete <- matrix(ncol=SIM,nrow=n) 
    }else{  
      X_Sdata_complete <- matrix(ncol=SIM,nrow=n)
      Zm_Sdata_complete <- matrix(ncol=SIM,nrow=n)
      Zp_Sdata_complete <- matrix(ncol=SIM,nrow=n)
      S_Sdata_complete <- matrix(ncol=SIM,nrow=n)
      Y_Sdata_complete <- matrix(ncol=SIM,nrow=n)
  
      X_Ddata_complete <- matrix(ncol=SIM,nrow=n)
      Zm_Ddata_complete <- matrix(ncol=SIM,nrow=n)
      Zp_Ddata_complete <- matrix(ncol=SIM,nrow=n)}
  
  betahat_first_complete <- matrix(ncol=SIM,nrow=7)
  betahat_first_ridge_complete <- matrix(ncol=SIM,nrow=7)
  beta_end1_complete <- matrix(ncol=SIM, nrow=3)
  beta_end1_ridge_complete <- matrix(ncol=SIM, nrow=3)
  naive_beta_complete <- matrix(ncol=SIM, nrow=3)
  beta_ZlmX_complete <- matrix(ncol=SIM, nrow=5)
  beta_Z2lmX_complete <- matrix(ncol=SIM, nrow=5)
  beta_ZlmX_ridge_complete <- matrix(ncol=SIM, nrow=5)
  beta_Z2lmX_ridge_complete <- matrix(ncol=SIM, nrow=5)
  
  for(sim in 1:SIM)
  {
    print(sim)
    
    #reproducibility
    set.seed(sim)
    
    #data generating processes
    if(empty==F){
      U<-rnorm(n)
      Zp<-2*U+rnorm(n)
      X<-Zp+rnorm(n=n)
      Zm<-X+2*rnorm(n)+2*U
      S<-Zm+X>5
      Y<-0.5*X^2+2*Zm+3*rnorm(n)+2*U
    
      Zm_S<-Zm[S==1]
      Zp_S<-Zp[S==1]
      X_S<-X[S==1]
      Y_S<-Y[S==1]
      }else{
        U_Sdata<-rnorm(n)
        Zp_Sdata<-2*U_Sdata+rnorm(n)
        X_Sdata<-Zp_Sdata+rnorm(n=n)
        Zm_Sdata<-X_Sdata+2*rnorm(n)+2*U_Sdata
        S_Sdata<-Zm_Sdata+X_Sdata>5
        Y_Sdata<-0.5*X_Sdata^2+2*Zm_Sdata+3*rnorm(n)+2*U_Sdata
      
        Zm_S<-Zm_Sdata[S_Sdata==1]
        Zp_S<-Zp_Sdata[S_Sdata==1]
        X_S<-X_Sdata[S_Sdata==1]
        Y_S<-Y_Sdata[S_Sdata==1]
      
        U_Ddata<-rnorm(n)
        Zp_Ddata<-2*U_Ddata+rnorm(n)
        X_Ddata<-Zp_Ddata+rnorm(n=n)
        Zm_Ddata<-X_Ddata+2*rnorm(n)+2*U_Ddata
      
        Zp<-Zp_Ddata
        X<-X_Ddata
        Zm<-Zm_Ddata}
    
    #save realisations from all simulation runs
    if(empty==F){
      X_complete[,sim]<-X
      Zp_complete[,sim]<-Zp
      Zm_complete[,sim]<-Zm
      S_complete[,sim]<-S
      Y_complete[,sim]<-Y
      }else{
        X_Sdata_complete[,sim]<-X_Sdata
        Zp_Sdata_complete[,sim]<-Zp_Sdata
        Zm_Sdata_complete[,sim]<-Zm_Sdata
        S_Sdata_complete[,sim]<-S_Sdata
        Y_Sdata_complete[,sim]<-Y_Sdata
      
        X_Ddata_complete[,sim]<-X_Ddata
        Zp_Ddata_complete[,sim]<-Zp_Ddata
        Zm_Ddata_complete[,sim]<-Zm_Ddata}
    
    #first step ridge
    lambda_seq <- 10^seq(2, -2, by = -.1)
    ridge_cv <- cv.glmnet(cbind(X_S,I(X_S^2),Zp_S,I(Zp_S^2),Zm_S,I(Zm_S^2)), Y_S, alpha = 0, lambda = lambda_seq)
    best_lambda <- ridge_cv$lambda.min
    best_ridge <- glmnet(cbind(X_S,I(X_S^2),Zp_S,I(Zp_S^2),Zm_S,I(Zm_S^2)), Y_S, alpha = 0, lambda = best_lambda)
    betahat_first_ridge<-as.vector(coef(best_ridge))
    
    #first step OLS
    lm_fit <- lm(Y_S~X_S+I(X_S^2)+Zp_S+I(Zp_S^2)+Zm_S+I(Zm_S^2))
    betahat_first<-coef(lm_fit)
    
    #calculate target variable for second step of RR
    Y_second_ridge<-(cbind(rep(1,n),X,X^2,Zp,Zp^2,Zm,Zm^2))%*%as.vector(betahat_first_ridge)
    Y_second<-(cbind(rep(1,n),X,X^2,Zp,Zp^2,Zm,Zm^2))%*%as.vector(betahat_first)
    
    #final RR estimates
    beta_end1_ridge<-lm(Y_second_ridge~1+X+I(X^2))$coef
    beta_end1<-lm(Y_second~1+X+I(X^2))$coef
    
    #naive estimation only based on S=1 
    naive_beta<-lm(Y_S~1+X_S+I(X_S^2))$coef
    
    #second step of TSR OLS
    beta_ZlmX<-lm(Zm~1+X+I(X^2)+Zp+I(Zp^2))$coef
    beta_Z2lmX<-lm(Zm^2~1+X+I(X^2)+Zp+I(Zp^2))$coef
    
    #second step of TSR ridge
    lambda_seq <- 10^seq(2, -2, by = -.1)
    ridge_cv <- cv.glmnet(cbind(X,I(X^2),Zp,I(Zp^2)), Zm, alpha = 0, lambda = lambda_seq)
    best_lambda <- ridge_cv$lambda.min
    best_ridge <- glmnet(cbind(X,I(X^2),Zp,I(Zp^2)), Zm, alpha = 0, lambda = best_lambda)
    beta_ZlmX_ridge<-as.vector(coef(best_ridge))
    
    lambda_seq <- 10^seq(2, -2, by = -.1)
    ridge_cv <- cv.glmnet(cbind(X,I(X^2),Zp,I(Zp^2)), Zm^2, alpha = 0, lambda = lambda_seq)
    best_lambda <- ridge_cv$lambda.min
    best_ridge <- glmnet(cbind(X,I(X^2),Zp,I(Zp^2)), Zm^2, alpha = 0, lambda = best_lambda)
    beta_Z2lmX_ridge<-as.vector(coef(best_ridge))
    
    #save coefficient vectors from all simulation runs
    betahat_first_complete[,sim] <- betahat_first
    betahat_first_ridge_complete[,sim] <- betahat_first_ridge
    beta_end1_complete[,sim] <- beta_end1
    beta_end1_ridge_complete[,sim] <- beta_end1_ridge
    naive_beta_complete[,sim] <- naive_beta
    beta_ZlmX_complete[,sim] <- beta_ZlmX
    beta_Z2lmX_complete[,sim] <- beta_Z2lmX
    beta_ZlmX_ridge_complete[,sim] <- beta_ZlmX_ridge
    beta_Z2lmX_ridge_complete[,sim] <- beta_Z2lmX_ridge
  } 
  
  #return the relevant results  
  if(empty==F){
    return(list(X_complete = X_complete, Zp_complete=Zp_complete, Zm_complete = Zm_complete, S_complete = S_complete, Y_complete = Y_complete,
              betahat_first_complete = betahat_first_complete, betahat_first_ridge_complete = betahat_first_ridge_complete,
              beta_end1_complete = beta_end1_complete, beta_end1_ridge_complete = beta_end1_ridge_complete,
              naive_beta_complete = naive_beta_complete, 
              beta_ZlmX_complete = beta_ZlmX_complete, beta_Z2lmX_complete = beta_Z2lmX_complete,
              beta_ZlmX_ridge_complete = beta_ZlmX_ridge_complete, beta_Z2lmX_ridge_complete = beta_Z2lmX_ridge_complete,
              SIM=SIM,n=n,empty=empty))
    }else{
      return(list(X_Sdata_complete = X_Sdata_complete, Zp_Sdata_complete=Zp_Sdata_complete, Zm_Sdata_complete = Zm_Sdata_complete, S_Sdata_complete = S_Sdata_complete, Y_Sdata_complete = Y_Sdata_complete,
                X_Ddata_complete = X_Ddata_complete, Zp_Ddata_complete=Zp_Ddata_complete, Zm_Ddata_complete = Zm_Ddata_complete,
                betahat_first_complete = betahat_first_complete, betahat_first_ridge_complete = betahat_first_ridge_complete,
                beta_end1_complete = beta_end1_complete, beta_end1_ridge_complete = beta_end1_ridge_complete,
                naive_beta_complete = naive_beta_complete, 
                beta_ZlmX_complete = beta_ZlmX_complete, beta_Z2lmX_complete = beta_Z2lmX_complete,
                beta_ZlmX_ridge_complete = beta_ZlmX_ridge_complete, beta_Z2lmX_ridge_complete = beta_Z2lmX_ridge_complete,
                SIM=SIM,n=n,empty=empty))}
}


result500 <- sim1(n=500,SIM=100,empty=F)
result2000 <- sim1(n=2000,SIM=100,empty=F)

empty_result500 <- sim1(n=500,SIM=100,empty=T)
empty_result2000 <- sim1(n=2000,SIM=100,empty=T)


#plot with curves for one case (n=2000, S \cap D= empyset)
result<-empty_result2000
SIM<-result$SIM
n<-result$n
empty<-result$empty

if(empty==F){
  X <- as.vector(result$X_complete)
  Zp <- as.vector(result$Zp_complete)
  Zm <- as.vector(result$Zm_complete)
  S <- as.vector(result$S_complete)
  Y <- as.vector(result$Y_complete)
  X_S<-X[S==1]
  }else{
    X <- as.vector(result$X_Ddata_complete)
    X_S <- as.vector(result$X_Sdata_complete[result$S_Sdata_complete==1])
  
    Zp <- as.vector(result$Zp_Ddata_complete)
    Zm <- as.vector(result$Zm_Ddata_complete)
    S <- as.vector(result$S_Ddata_complete)
    Y <- as.vector(result$Y_Ddata_complete)}

#x values for evaluating 95% areas
X_est <- seq(from = -15,to = 15, by = 0.2 )   
l_X_est <- length(X_est)

values_est_naive <- matrix(ncol=SIM,nrow=l_X_est)
values_est_est1 <- matrix(ncol=SIM,nrow=l_X_est)
values_est_est1_ridge <- matrix(ncol=SIM,nrow=l_X_est)
values_est_est3 <- matrix(ncol=SIM,nrow=l_X_est)
values_est_est3_ridge <- matrix(ncol=SIM,nrow=l_X_est)


X_est<-as.matrix(X_est)

for (sim in 1:SIM) {
  print(sim)
  
  betahat_first <- as.vector(result$betahat_first_complete[,sim])
  betahat_first_ridge <- as.vector(result$betahat_first_ridge_complete[,sim])
  beta_end1 <- as.vector(result$beta_end1_complete[,sim])
  beta_end1_ridge <- as.vector(result$beta_end1_ridge_complete[,sim])
  naive_beta <- as.vector(result$naive_beta_complete[,sim])
  beta_ZlmX <- as.vector(result$beta_ZlmX_complete[,sim])
  beta_Z2lmX <- as.vector(result$beta_Z2lmX_complete[,sim])
  beta_ZlmX_ridge <- as.vector(result$beta_ZlmX_ridge_complete[,sim])
  beta_Z2lmX_ridge <- as.vector(result$beta_Z2lmX_ridge_complete[,sim])
  
  if(empty==F){Zp<-result$Zp_complete[,sim]
  }else{Zp <- as.vector(result$Zp_Ddata_complete)}
  
  #E[Y|do(X)] 
  original<-function(x)
  {1/2*x^2+2*x}
  
  #RR ridge
  est1_ridge<-function(x)
  {beta_end1_ridge[1]+beta_end1_ridge[2]*x+beta_end1_ridge[3]*x^2}
  
  #RR
  est1<-function(x)
  {beta_end1[1]+beta_end1[2]*x+beta_end1[3]*x^2}    
  
  #naive
  naive<-function(x)
  {naive_beta[1]+naive_beta[2]*x+naive_beta[3]*x^2}
  
  #TSR ridge
  est3_ridge<-function(x)
  {return(betahat_first_ridge[1]+betahat_first_ridge[2]*x+betahat_first_ridge[3]*x^2+
             betahat_first_ridge[4]*mean(Zp)+
             betahat_first_ridge[5]*mean(Zp^2)+
             betahat_first_ridge[6]*(beta_ZlmX_ridge[1]+beta_ZlmX_ridge[2]*x+beta_ZlmX_ridge[3]*x^2+beta_ZlmX_ridge[4]*mean(Zp)+beta_ZlmX_ridge[5]*mean(Zp^2))+
             betahat_first_ridge[7]*(beta_Z2lmX_ridge[1]+beta_Z2lmX_ridge[2]*x+beta_Z2lmX_ridge[3]*x^2+beta_Z2lmX_ridge[4]*mean(Zp)+beta_Z2lmX_ridge[5]*mean(Zp^2)))} 
  
  #TSR
  est3<-function(x)
  {return(betahat_first[1]+betahat_first[2]*x+betahat_first[3]*x^2+
             betahat_first[4]*mean(Zp)+
             betahat_first[5]*mean(Zp^2)+
             betahat_first[6]*(beta_ZlmX[1]+beta_ZlmX[2]*x+beta_ZlmX[3]*x^2+beta_ZlmX[4]*mean(Zp)+beta_ZlmX[5]*mean(Zp^2))+
             betahat_first[7]*(beta_Z2lmX[1]+beta_Z2lmX[2]*x+beta_Z2lmX[3]*x^2+beta_Z2lmX[4]*mean(Zp)+beta_Z2lmX[5]*mean(Zp^2)))} 
  values_est_naive[,sim]<-apply(X_est,1,naive)
  values_est_est1[,sim]<-apply(X_est,1,est1)
  values_est_est1_ridge[,sim]<-apply(X_est,1,est1_ridge)
  values_est_est3[,sim]<-apply(X_est,1,est3)
  values_est_est3_ridge[,sim]<-apply(X_est,1,est3_ridge)
  
}

alpha<-0.05

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


est1_mean<-apply(values_est_est1,1,mean)
est1_ridge_mean<-apply(values_est_est1_ridge,1,mean)
est3_mean<-apply(values_est_est3,1,mean)
est3_ridge_mean<-apply(values_est_est3_ridge,1,mean)

green_rgb <- col2rgb("green")
darkgreen_rgb <- col2rgb("darkgreen")
middle_green_rgb <- (green_rgb + darkgreen_rgb) / 2
middle_green <- rgb(middle_green_rgb[1], middle_green_rgb[2], middle_green_rgb[3], maxColorValue = 255)

#OLS
ggplot() +
  geom_ribbon(aes(x = X_est, ymin = q_lower_est1, ymax = q_upper_est1, fill = "RR"), alpha = 0.8) +
  geom_ribbon(aes(x = X_est, ymin = q_lower_est3, ymax = q_upper_est3, fill = "TSR"), alpha = 0.8) +
  geom_boxplot(aes(x = X, y = -32), width = 7.5, fill = "grey", alpha = 0.7) +
  geom_boxplot(aes(x = X_S, y = -20), width = 7.5, fill = "grey", alpha = 0.7) +
  #geom_line(aes(x = X_est, y = apply(X_est,1,original), color = "E[Y|do(X=x)]"), linewidth = 1.3) +
  geom_line(aes(x = X_est, y = est1_mean, color = "RR"), size = 1.5)+
  geom_line(aes(x = X_est, y = q_lower_est1, color = "RR"),linetype="dashed", size = 1)+
  geom_line(aes(x = X_est, y = q_upper_est1, color = "RR"),linetype="dashed", size = 1)+
  geom_line(aes(x = X_est, y = est3_mean, color = "TSR"), size = 1.5)+
  geom_line(aes(x = X_est, y = q_lower_est3, color = "TSR"),linetype="dashed", size = 1)+
  geom_line(aes(x = X_est, y = q_upper_est3, color = "TSR"),linetype="dashed", size = 1)+
  geom_line(aes(x = X_est, y = apply(X_est,1,original), color = "ground truth"), linewidth = 1.3) +
  theme_minimal() +
  coord_cartesian(ylim = c(-40, 100),xlim=c(-10,10))+
  scale_fill_manual(values = c( "TSR" = "lightgreen", "RR" = "lightblue"), name = NULL) +
  scale_color_manual(values = c("ground truth" = "black", "RR"="blue","TSR"=middle_green), name = NULL) +
  theme(legend.position = "top")+
  labs(x="x",y=expression(hat(E) * group("[", Y * group("|", do(X == x), ""), "]")))+
  theme(
    axis.title.x = element_text(size = 20),  
    axis.title.y = element_text(size = 20),  
    legend.title = element_text(size = 20),  
    legend.text = element_text(size = 20),  
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.position = "right",
    legend.spacing = unit(2, "cm"),
    legend.key.height = unit(1, "cm"),  
    legend.key.width = unit(1, "cm")
  )+
  guides(fill = guide_legend(order = 2),
         color = guide_legend(order = 2))


#ridge
ggplot() +
  geom_ribbon(aes(x = X_est, ymin = q_lower_est1_ridge, ymax = q_upper_est1_ridge, fill = "RR (ridge)"), alpha = 0.8) +
  geom_ribbon(aes(x = X_est, ymin = q_lower_est3_ridge, ymax = q_upper_est3_ridge, fill = "TSR (ridge)"), alpha = 0.8) +
  geom_boxplot(aes(x = X, y = -32), width = 7.5, fill = "grey", alpha = 0.7) +
  geom_boxplot(aes(x = X_S, y = -20), width = 7.5, fill = "grey", alpha = 0.7) +
  #geom_line(aes(x = X_est, y = apply(X_est,1,original), color = "E[Y|do(X=x)]"), linewidth = 1.3) +
  geom_line(aes(x = X_est, y = est1_ridge_mean, color = "RR (ridge)"), size = 1.5)+
  geom_line(aes(x = X_est, y = q_lower_est1_ridge, color = "RR (ridge)"),linetype="dashed", size = 1)+
  geom_line(aes(x = X_est, y = q_upper_est1_ridge, color = "RR (ridge)"),linetype="dashed", size = 1)+
  geom_line(aes(x = X_est, y = est3_ridge_mean, color = "TSR (ridge)"), size = 1.5)+
  geom_line(aes(x = X_est, y = q_lower_est3_ridge, color = "TSR (ridge)"),linetype="dashed", size = 1)+
  geom_line(aes(x = X_est, y = q_upper_est3_ridge, color = "TSR (ridge)"),linetype="dashed", size = 1)+
  geom_line(aes(x = X_est, y = apply(X_est,1,original), color = "ground truth"), linewidth = 1.3) +
  theme_minimal() +
  coord_cartesian(ylim = c(-40, 100),xlim=c(-10,10))+
  scale_fill_manual(values = c( "TSR (ridge)" = "lightgreen", "RR (ridge)" = "lightblue"), name = NULL) +
  scale_color_manual(values = c("ground truth" = "black", "RR (ridge)"="blue","TSR (ridge)"=middle_green), name = NULL) +
  theme(legend.position = "top")+
  labs(x="x",y=expression(hat(E) * group("[", Y * group("|", do(X == x), ""), "]")))+
  theme(
    axis.title.x = element_text(size = 20),  
    axis.title.y = element_text(size = 20),  
    legend.title = element_text(size = 20),  
    legend.text = element_text(size = 20),  
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.position = "right",
    legend.spacing = unit(2, "cm"),
    legend.key.height = unit(1, "cm"),  
    legend.key.width = unit(1, "cm")
  )+
  guides(fill = guide_legend(order = 2),
         color = guide_legend(order = 2))


