library(coda)
library(invgamma)
library(MASS)
library(gtools)
library(matlib)
library(Matrix)

#################################
##     Exercice 14 page 80     ##
#################################


Gibbs_School = function(N_chain = 500, init){
  
  
  chain = matrix(NA, N_chain + 1, ncol = 136)
  colnames(chain) = c('theta',
                      'phi',
                      paste0('gamma_', 1:3),
                      paste0('beta_', 1:8),
                      paste0('T_1_', 1:3), # attention, rempli par ligne
                      paste0('T_2_', 1:3),
                      paste0('T_3_', 1:3),
                      paste0('alpha_1_', 1:38),# attention, rempli par ligne
                      paste0('alpha_2_', 1:38),
                      paste0('alpha_3_', 1:38))
  
  # init de 
  chain[1,] = init
  
  for(i in 2:N_chain){
    
    # remise en forme
    theta = chain[i-1,1]
    phi = chain[i-1,2]
    gamma = chain[i-1,3:5]
    beta = chain[i-1,6:13]
    T = matrix(chain[i-1,14:22], ncol = 3, nrow = 3, byrow = TRUE)
    alpha = matrix(chain[i-1,23:136], ncol = M, nrow = 3, byrow = TRUE)
    
    
    ## maj de gamma
    Sigma_gamma_star = inv(M*T + diag(10**4, ncol = 3, nrow = 3))
    Mu_gamma_star = Sigma_gamma_star%*%(T%*%rowSums(alpha)) # moy de gamma0 c'est 0
    gamma = mvrnorm(1, Mu_gamma_star, Sigma_gamma_star)
    
    
    ## maj de T
    somme = matrix(0,9,nrow=3,ncol=3)
    for (j in 1:M){
      a_m_g_j = alpha[,j] - gamma
      somme = somme + tcrossprod(a_m_g_j)
    }
    V = inv(somme + inv(R))
    T = rWishart(1, df = M+2, Sigma = V)[,,1]
    
    
    # updating theta
    prop = theta + runif(1,min=-1,max=1) # support englobe-t-il celui de g = exp(-1/2*(IJtheta+...))
    
    mu = rep(NaN, N)
    ltau = rep(NaN, N)
    for (i in 1:N){
      ltau[i] <- theta + phi * LRT[i] # log(tau_ij)
      
      j = school[i]
      mu[i] <- alpha[1, j] + alpha[2, j] * LRT[i]
      + alpha[3, j] * VR[i, 1] + beta[1] * LRT[i]**2
      + beta[2] * VR[i, 2] + beta[3] * Gender[i]
      + beta[4] * School_gender[i, 1] + beta[5] * School_gender[i, 2]
      + beta[6] * School_denom[i, 1] + beta[7] * School_denom[i, 2]
      + beta[8] * School_denom[i, 3]
    }
    
    somme_t = sum(exp(phi*LRT) * (Y-mu)**2) #somme dans la densité coditionnelle de theta
    
    top = -0.5*(10**4*prop**2 + exp(prop)*somme_t - M*N*prop)
    bottom = -0.5*(10**4*theta**2 + exp(theta)*somme_t - M*N*theta)
    
    acc_prob = exp(top-bottom) # <- sur de ça ?
    if (runif(1,min=0,max=1) < acc_prob){
      theta = prop
    } # pas besoin du else dans mon cas
    
    
    # Updating phi
    prop = phi + runif(1,min=-1,max=1) # même question au niveau du suppport
    
    top = prop*sum(LRT) - 0.5*(prop**2*10**4 + exp(theta)*somme_t)
    bottom = phi*sum(LRT) - 0.5*(phi**2*10**4 + exp(theta)*somme_t)
    
    acc_prob = exp(top-bottom) 
    if (runif(1,min=0,max=1) < acc_prob){
      phi = prop
    }# pas besoin du else dans mon cas
    
    
    
    ## maj de alpha
    
    
    for(j in 1:M){
      nb_in_school_j = sum(datadf[,'school'] == j)
      invT = inv(T)
      
      # mais alpha et y_ij pas de la même dimension quand même
      #alpha_j_tilde = rep(alpha[,j], nb_in_school_j)
      #gamma_tilde = rep(gamma, nb_in_school_j)
      #Sigma_tilde = bdiag(lapply(1:nb_in_school_j, \(x){T}))# matrice par block (I répétitions de T)
      
      # maj de alpha 1
      mu_prior_alpha1j = gamma[1] + invT[1,-1]%*%T[-1,-1]%*%(alpha[-1,j]-gamma[-1])
      sigma2_prior_alpha1j = invT[1,1] - matrix(invT[1,-1], ncol=2, nrow = 1)%*%T[-1,-1]%*%matrix(invT[-1,1], ncol = 1, nrow = 2) #<- variance négative
      
     
      
    }
    
    
    
  }
  
  
  
  return(chain)
}


# data
source("schools.data.R")
datadf = cbind(data.frame(Y,school,Gender,LRT),School_gender,School_denom,VR)
colnames(datadf) = c("Y","school","Gender","LRT","School_gender1","School_gender2","School_denom1",
                     "School_denom2","School_denom3","VR1","VR2")
View(datadf)



R = matrix(c(1/10 , 1/200, 1/200, 
             1/200, 1/10 , 1/200, 
             1/200, 1/200, 1/10 ), 
           nrow = 3, ncol = 3, byrow = TRUE)


theta = rnorm(1, sd = 10**2)
phi = rnorm(1, sd = 10**2)
gamma = rnorm(3, sd = 10**2)
beta = rnorm(8, sd = 10**2)
T = rWishart(1, df = 3, Sigma = R)[,,1]
alpha = mvrnorm(38, mu = gamma, Sigma = T)

init = c(theta, phi, gamma, beta, T, alpha)






oups = Gibbs_Galaxies(1000, alpha = 0.1, beta = 0.1, mu0 = 20, sigma2_0 = 100, K = 6)