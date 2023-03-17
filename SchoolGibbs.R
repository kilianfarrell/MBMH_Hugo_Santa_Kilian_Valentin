#/---- Schools gibbs sampler ----/
library(mvtnorm)
library(purrr)
library(stats)
library(matlib)

# data
source("schools.data.R")
datadf = cbind(data.frame(Y,school,Gender,LRT),School_gender,School_denom,VR)
colnames(datadf) = c("Y","school","Gender","LRT","School_gender1","School_gender2","School_denom1",
                     "School_denom2","School_denom3","VR1","VR2")
View(datadf)

# init
theta <- 0

phi <- 0.01
gamma <- c(0, 0, 0)
beta <- c(0, 0, 0, 0, 0, 0, 0, 0)
Omega <- structure(c(10, -3, -3, -3, 135, -65, -3, -65, 135), .Dim = as.integer(c(3, 
                                                                         3)))
alpha <- structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), .Dim = as.integer(c(38, 
                                                                          3)))

schools_gibbs <- function(Nchain,theta_init,phi_init,gamma_init,beta_init,Omega_init,alpha_init,
                          g0,Sg0,Y,school,LRT,VR,Gender,School.gender,School.denom){ 
  n = 1978
  M = 38
  
  # To keep in chain : Beta, gamma, phi, theta
  Beta = matrix(NaN,nrow=Nchain,ncol=8,dimnames = list(c(), paste0("Beta",1:8)))
  gamma = matrix(NaN,nrow=Nchain,ncol=3,dimnames = list(c(), paste0("gamma",1:3)))
  phi = rep(NaN,Nchain)
  thetha = rep(NaN,Nchain)
  
  # For updating purposes : Omega, alpha_k
  Omega = NaN
  alpha = matrix(NaN,nrow = M,ncol=3,dimnames = list(c(), paste0("alpha",1:3)))
  
  # Initialisation : 
  Beta[1,] = beta_init
  gamma[1,] = gamma_init
  phi[1] = phi_init
  theta[1] = theta_init
  
  Omega = Omega_init
  T = inv(Omega)
  alpha = alpha_init
  
  tau0 = 0.0001
  
  
  for (t in 2:Nchain){
    # Updating  gamma #dim(Sg0) =3x3, dim(g0) = 3x1 
    Sig_star = M*T + Sg0
    Gamma_n = inv(Sig_star)%*%(M*T%*%matrix(colMeans(alpha),nrow = 3,ncol=1) + 
                                 sg0%*%matrix(g0,nrow=3,ncol=1)) 
    gamma[t,] = rmvnorm(Gamma_n,Sig_star)
    
    # Updtaing T
    somme = matrix(0,9,nrow=3,ncol=3)
    for (j in 1:M){
      a_m_g_j = alpha[j,] - gamma[t,]
      somme = somme + matrix(a_m_g_j,nrow=3,ncol=1)%*%matrix(a_m_g_j,nrow=3,ncol=1)
    }
    V = inv(somme + inv(R))
    T = rWishart(1,df = M+2,Sigma = V)[,,1]
    
    # updating theta
    prop = theta[t-1] + runif(1,min=-1,max=1) # support englobe-t-il celui de g = exp(-1/2*(IJtheta+...))
    
   
    mu = rep(NaN,n)
    ltau = rep(NaN,n)
    for (i in 1:n){
      ltau[i] <- theta[t-1] + phi[t-1] * LRT[i] # log(tau_ij)
      #sigma2[i] <- 1 / tau[i]
      LRT2[i] <- LRT[i] * LRT[i] #(LRT_ij)^2
      
      j = school[i]
      mu[i] <- alpha[j, 1] + alpha[j, 2] * LRT[i]
        + alpha[j, 3] * VR[i, 1] + Beta[t-1,1] * LRT2[i]
        + Beta[t-1,2] * VR[p, 2] + Beta[t-1,3] * Gender[i]
        + Beta[t-1,4] * School.gender[i, 1] + Beta[t-1,5] * School.gender[i, 2]
        + Beta[t-1,6] * School.denom[i, 1] + Beta[t-1,7] * School.denom[i, 2]
        + Beta[t-1,8] * School.denom[i, 3]
    }
    
    somme_t = sum(exp(phi[t-1]*LRT)*(Y-mu)^2) #somme dans la densité coditionnelle de theta
    
    top = -0.5*(- M*n*prop+prop^2*tau0+exp(prop)*somme_t)
    bottom = -0.5*(- M*n*theta[t-1]+theta[t-1]^2*tau0+exp(theta[t-1])*somme_t)
    
    acc_prob = exp(top-bottom) 
    if (runif(1,min=0,max=1) < acc_prob){
      theta[t] = prop
    }else{
      theta[t] = theta[t-1]
    }
    
    
    # Updating phi
    prop = phi[t-1] + runif(1,min=-1,max=1) # même question au niveau du suppport
    
    top = prop*sum(LRT) -0.5*(prop^2*tau0+exp(theta[t])*somme_t)
    bottom = phi[t-1]*sum(LRT) -0.5*(phi[t-1]^2*tau0+exp(theta[t])*somme_t)
    
    acc_prob = exp(top-bottom) 
    if (runif(1,min=0,max=1) < acc_prob){
      phi[t] = prop
    }else{
      phi[t] = phi[t-1]
    }
    
    #Updating alpha
    #PARTIE HUGO 
    
    
  }
  
  return(list(Beta,gamma,phi,theta))  
}


# Brouillon

matrix(NaN,nrow=5,ncol=8,dimnames = list(c(), paste0(paste0("Beta",1:8))))

Y[which(school == 1)]


datadf[which(datadf["school"] == 1),"LRT"]
alpha

School_gender


























