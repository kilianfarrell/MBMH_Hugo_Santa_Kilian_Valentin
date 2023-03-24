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
    Sigma_gamma_star = inv(M*T + diag(10**-4, ncol = 3, nrow = 3))
    Mu_gamma_star = Sigma_gamma_star%*%(T%*%rowSums(alpha)) # moy de gamma0 c'est 0
    gamma = mvrnorm(1, Mu_gamma_star, Sigma_gamma_star)
    
    
    ## maj de T
    somme = matrix(0,nrow=3,ncol=3)
    for (j in 1:M){
      a_j_moins_gamma = alpha[,j] - gamma
      somme = somme + tcrossprod(a_j_moins_gamma)
    }
    V = inv(somme + inv(R))
    #V = V + diag(10^-50, ncol = 3, nrow = 3)#print(V)
    T = rWishart(1, df = M + 3, Sigma = V)[,,1]
    
    
    ## maj de theta
    
    mu = rep(NaN, N)
    for (j in 1:M){
      eleve_in_school_j = which(datadf[,'school'] == j)
      
      mu_j = alpha[1,j]                              + 
        alpha[2,j] * LRT[eleve_in_school_j]          + 
        alpha[3,j] * VR[eleve_in_school_j,1]         +
        beta[1] * LRT[eleve_in_school_j]^2           + 
        beta[2] * VR[eleve_in_school_j,2]            +     
        beta[3] * Gender[eleve_in_school_j]          + 
        beta[4] * School_gender[eleve_in_school_j,1] +
        beta[5] * School_gender[eleve_in_school_j,2] + 
        beta[6] * School_denom[eleve_in_school_j,1]  +
        beta[7] * School_denom[eleve_in_school_j,2]  + 
        beta[8] * School_denom[eleve_in_school_j,3]
      
      mu[eleve_in_school_j] = mu_j
    }
    
    prop_sd = 0.1
    prop = theta + rnorm(1, 0, prop_sd)
    
    top =    -0.5*(-(prop*N) + (prop*sqrt(0.0001))^2 + (exp(prop)*sum(exp(phi*LRT)%*%(Y - mu)^2)))
    bottom = -0.5*(-(theta*N) + (theta*sqrt(0.0001))^2 + (exp(theta)*sum(exp(phi*LRT)%*%(Y - mu)^2)))
    
    if (runif(1) < exp(top-bottom)){
      theta = prop
    }
    
    
    ## maj de phi
    
    prop_sd = 0.1
    prop = phi + rnorm(1, 0, prop_sd)
    
    top =    -0.5*(-(prop*sum(LRT)) + (prop*sqrt(0.0001))^2 + (exp(theta)*sum(exp(prop*LRT)%*%(Y - mu)^2)))
    bottom = -0.5*(-(phi*sum(LRT)) + (phi*sqrt(0.0001))^2 + (exp(theta)*sum(exp(phi*LRT)%*%(Y - mu)^2)))
    
    if (runif(1) < exp(top-bottom)){
      phi = prop
    }
    
    
    ## maj de alpha
    
    for(j in 1:M){
      
      eleve_in_school_j = which(datadf[,'school'] == j)
      nb_in_school_j = sum(datadf[,'school'] == j)
      
      tau_j = exp(theta + phi * LRT[eleve_in_school_j])
      LRT_j = LRT[eleve_in_school_j]
      VR1_j = VR[eleve_in_school_j, 1]
      Y_j = Y[eleve_in_school_j]
      
      mu_j_sans_alphaj = beta[1] * LRT[eleve_in_school_j]^2 + beta[2] * VR[eleve_in_school_j, 2]     +     
        beta[3] * Gender[eleve_in_school_j]           + beta[4] * School_gender[eleve_in_school_j, 1] +
        beta[5] * School_gender[eleve_in_school_j, 2] + beta[6] * School_denom[eleve_in_school_j, 1]  +
        beta[7] * School_denom[eleve_in_school_j, 2]  + beta[8] * School_denom[eleve_in_school_j, 3]
      
      somme_sigma = 0
      somme_mu = 0
      for(k in 1:nb_in_school_j){
        somme_sigma = somme_sigma + tau_j[k] * tcrossprod(c(1, LRT_j[k], VR1_j[k]))
        somme_mu = somme_mu + c(1, LRT_j[k], VR1_j[k]) * (Y_j[k] - mu_j_sans_alphaj[k]) * tau_j[k]
      }
      
      
      sigma_star_alphaj = inv(T + somme_sigma)
      mu_star_alpha_j = sigma_star_alphaj %*% (T %*% gamma + matrix(somme_mu, ncol = 1, nrow = 3))
      
      alpha[,j] = mvrnorm(mu = mu_star_alpha_j, Sigma = sigma_star_alphaj)
    }
    
    
    ## maj des betas
    
    tau = exp(theta + phi * LRT)
    
    
    # maj de beta1
    mu_sans_beta1 = rep(NaN, N)
    for (j in 1:M){
      eleve_in_school_j = which(datadf[,'school'] == j)
      
      mu_j_sans_beta1 = alpha[1,j]                   + 
        alpha[2,j] * LRT[eleve_in_school_j]          + 
        alpha[3,j] * VR[eleve_in_school_j,1]         +
        #beta[1] * LRT[eleve_in_school_j]^2           + 
        beta[2] * VR[eleve_in_school_j,2]            +     
        beta[3] * Gender[eleve_in_school_j]          + 
        beta[4] * School_gender[eleve_in_school_j,1] +
        beta[5] * School_gender[eleve_in_school_j,2] + 
        beta[6] * School_denom[eleve_in_school_j,1]  +
        beta[7] * School_denom[eleve_in_school_j,2]  + 
        beta[8] * School_denom[eleve_in_school_j,3]
      
      mu_sans_beta1[eleve_in_school_j] = mu_j_sans_beta1
    }
    
    sigma2_star_beta1 = 1/(10**-4 + sum(tau * LRT^4))
    mu_star_beta1 = sigma2_star_beta1 * (sum(tau * LRT^2 * (Y - mu_sans_beta1)))
    
    beta[1] = rnorm(1, mean = mu_star_beta1, sd = sqrt(sigma2_star_beta1))
    
    
    # maj de beta2
    mu_sans_beta2 = rep(NaN, N)
    for (j in 1:M){
      eleve_in_school_j = which(datadf[,'school'] == j)
      
      mu_j_sans_beta2 = alpha[1,j]                   + 
        alpha[2,j] * LRT[eleve_in_school_j]          + 
        alpha[3,j] * VR[eleve_in_school_j,1]         +
        beta[1] * LRT[eleve_in_school_j]^2           + 
        #beta[2] * VR[eleve_in_school_j,2]            +     
        beta[3] * Gender[eleve_in_school_j]          + 
        beta[4] * School_gender[eleve_in_school_j,1] +
        beta[5] * School_gender[eleve_in_school_j,2] + 
        beta[6] * School_denom[eleve_in_school_j,1]  +
        beta[7] * School_denom[eleve_in_school_j,2]  + 
        beta[8] * School_denom[eleve_in_school_j,3]
      
      mu_sans_beta2[eleve_in_school_j] = mu_j_sans_beta2
    }
    
    sigma2_star_beta2 = 1/(10**-4 + sum(tau * VR[,2]^2))
    mu_star_beta2 = sigma2_star_beta2 * (sum(tau * VR[,2] * (Y - mu_sans_beta2)))
    
    beta[2] = rnorm(1, mean = mu_star_beta2, sd = sqrt(sigma2_star_beta2))
    
    
    # maj de beta3
    mu_sans_beta3 = rep(NaN, N)
    for (j in 1:M){
      eleve_in_school_j = which(datadf[,'school'] == j)
      
      mu_j_sans_beta3 = alpha[1,j]                   + 
        alpha[2,j] * LRT[eleve_in_school_j]          + 
        alpha[3,j] * VR[eleve_in_school_j,1]         +
        beta[1] * LRT[eleve_in_school_j]^2           + 
        beta[2] * VR[eleve_in_school_j,2]            +     
        #beta[3] * Gender[eleve_in_school_j]          + 
        beta[4] * School_gender[eleve_in_school_j,1] +
        beta[5] * School_gender[eleve_in_school_j,2] + 
        beta[6] * School_denom[eleve_in_school_j,1]  +
        beta[7] * School_denom[eleve_in_school_j,2]  + 
        beta[8] * School_denom[eleve_in_school_j,3]
      
      mu_sans_beta3[eleve_in_school_j] = mu_j_sans_beta3
    }
    
    sigma2_star_beta3 = 1/(10**-4 + sum(tau * Gender^2))
    mu_star_beta3 = sigma2_star_beta3 * (sum(tau * Gender * (Y - mu_sans_beta3)))
    
    beta[3] = rnorm(1, mean = mu_star_beta3, sd = sqrt(sigma2_star_beta3))
    
    
    # maj de beta4
    mu_sans_beta4 = rep(NaN, N)
    for (j in 1:M){
      eleve_in_school_j = which(datadf[,'school'] == j)
      
      mu_j_sans_beta4 = alpha[1,j]                   + 
        alpha[2,j] * LRT[eleve_in_school_j]          + 
        alpha[3,j] * VR[eleve_in_school_j,1]         +
        beta[1] * LRT[eleve_in_school_j]^2           + 
        beta[2] * VR[eleve_in_school_j,2]            +     
        beta[3] * Gender[eleve_in_school_j]          + 
        #beta[4] * School_gender[eleve_in_school_j,1] +
        beta[5] * School_gender[eleve_in_school_j,2] + 
        beta[6] * School_denom[eleve_in_school_j,1]  +
        beta[7] * School_denom[eleve_in_school_j,2]  + 
        beta[8] * School_denom[eleve_in_school_j,3]
      
      mu_sans_beta4[eleve_in_school_j] = mu_j_sans_beta4
    }
    
    sigma2_star_beta4 = 1/(10**-4 + sum(tau * School_gender[,1]^2))
    mu_star_beta4 = sigma2_star_beta4 * (sum(tau * School_gender[,1] * (Y - mu_sans_beta4)))
    
    beta[4] = rnorm(1, mean = mu_star_beta4, sd = sqrt(sigma2_star_beta4))
    
    
    # maj de beta5
    mu_sans_beta5 = rep(NaN, N)
    for (j in 1:M){
      eleve_in_school_j = which(datadf[,'school'] == j)
      
      mu_j_sans_beta5 = alpha[1,j]                   + 
        alpha[2,j] * LRT[eleve_in_school_j]          + 
        alpha[3,j] * VR[eleve_in_school_j,1]         +
        beta[1] * LRT[eleve_in_school_j]^2           + 
        beta[2] * VR[eleve_in_school_j,2]            +     
        beta[3] * Gender[eleve_in_school_j]          + 
        beta[4] * School_gender[eleve_in_school_j,1] +
        #beta[5] * School_gender[eleve_in_school_j,2] + 
        beta[6] * School_denom[eleve_in_school_j,1]  +
        beta[7] * School_denom[eleve_in_school_j,2]  + 
        beta[8] * School_denom[eleve_in_school_j,3]
      
      mu_sans_beta5[eleve_in_school_j] = mu_j_sans_beta5
    }
    
    sigma2_star_beta5 = 1/(10**-4 + sum(tau * School_gender[,2]^2))
    mu_star_beta5 = sigma2_star_beta5 * (sum(tau * School_gender[,2] * (Y - mu_sans_beta5)))
    
    beta[5] = rnorm(1, mean = mu_star_beta5, sd = sqrt(sigma2_star_beta5))
    
    
    # maj de beta6
    mu_sans_beta6 = rep(NaN, N)
    for (j in 1:M){
      eleve_in_school_j = which(datadf[,'school'] == j)
      
      mu_j_sans_beta6 = alpha[1,j]                   + 
        alpha[2,j] * LRT[eleve_in_school_j]          + 
        alpha[3,j] * VR[eleve_in_school_j,1]         +
        beta[1] * LRT[eleve_in_school_j]^2           + 
        beta[2] * VR[eleve_in_school_j,2]            +     
        beta[3] * Gender[eleve_in_school_j]          + 
        beta[4] * School_gender[eleve_in_school_j,1] +
        beta[5] * School_gender[eleve_in_school_j,2] + 
        #beta[6] * School_denom[eleve_in_school_j,1]  +
        beta[7] * School_denom[eleve_in_school_j,2]  + 
        beta[8] * School_denom[eleve_in_school_j,3]
      
      mu_sans_beta6[eleve_in_school_j] = mu_j_sans_beta6
    }
    
    sigma2_star_beta6 = 1/(10**-4 + sum(tau * School_denom[,1]^2))
    mu_star_beta6 = sigma2_star_beta6 * (sum(tau * School_denom[,1] * (Y - mu_sans_beta6)))
    
    beta[6] = rnorm(1, mean = mu_star_beta6, sd = sqrt(sigma2_star_beta6))
    
    
    # maj de beta7
    mu_sans_beta7 = rep(NaN, N)
    for (j in 1:M){
      eleve_in_school_j = which(datadf[,'school'] == j)
      
      mu_j_sans_beta7 = alpha[1,j]                   + 
        alpha[2,j] * LRT[eleve_in_school_j]          + 
        alpha[3,j] * VR[eleve_in_school_j,1]         +
        beta[1] * LRT[eleve_in_school_j]^2           + 
        beta[2] * VR[eleve_in_school_j,2]            +     
        beta[3] * Gender[eleve_in_school_j]          + 
        beta[4] * School_gender[eleve_in_school_j,1] +
        beta[5] * School_gender[eleve_in_school_j,2] + 
        beta[6] * School_denom[eleve_in_school_j,1]  +
        #beta[7] * School_denom[eleve_in_school_j,2]  + 
        beta[8] * School_denom[eleve_in_school_j,3]
      
      mu_sans_beta7[eleve_in_school_j] = mu_j_sans_beta7
    }
    
    sigma2_star_beta7 = 1/(10**-4 + sum(tau * School_denom[,2]^2))
    mu_star_beta7 = sigma2_star_beta7 * (sum(tau * School_denom[,2] * (Y - mu_sans_beta7)))
    
    beta[7] = rnorm(1, mean = mu_star_beta7, sd = sqrt(sigma2_star_beta7))
    
    
    # maj de beta8
    mu_sans_beta8 = rep(NaN, N)
    for (j in 1:M){
      eleve_in_school_j = which(datadf[,'school'] == j)
      
      mu_j_sans_beta8 = alpha[1,j]                   + 
        alpha[2,j] * LRT[eleve_in_school_j]          + 
        alpha[3,j] * VR[eleve_in_school_j,1]         +
        beta[1] * LRT[eleve_in_school_j]^2           + 
        beta[2] * VR[eleve_in_school_j,2]            +     
        beta[3] * Gender[eleve_in_school_j]          + 
        beta[4] * School_gender[eleve_in_school_j,1] +
        beta[5] * School_gender[eleve_in_school_j,2] + 
        beta[6] * School_denom[eleve_in_school_j,1]  +
        beta[7] * School_denom[eleve_in_school_j,2]  #+ 
        #beta[8] * School_denom[eleve_in_school_j,3]
      
      mu_sans_beta8[eleve_in_school_j] = mu_j_sans_beta8
    }
    
    sigma2_star_beta8 = 1/(10**-4 + sum(tau * School_denom[,3]^2))
    mu_star_beta8 = sigma2_star_beta8 * (sum(tau * School_denom[,3] * (Y - mu_sans_beta8)))
    
    beta[8] = rnorm(1, mean = mu_star_beta8, sd = sqrt(sigma2_star_beta8))
    
    
    
    # remise dans la chaine
    
    chain[i,] = c(theta, phi, gamma, beta, T, alpha)
    #print(chain[i,])
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


theta = 0.6#0 #rnorm(1, sd = 10**2)
phi = 0.01 #rnorm(1, sd = 10**2)
gamma = c(-0.7,0,1)#c(0, 0, 0) #rnorm(3, sd = 10**2)
beta = c(0, 0, 0, 0, 0, 0, 0, 0) #rnorm(8, sd = 10**2)
T = matrix(c(10, -3, -3, -3, 135, -65, -3, -65, 135), nrow = 3, ncol = 3, byrow = TRUE) #rWishart(1, df = 3, Sigma = R)[,,1]
alpha = structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), .Dim = as.integer(c(3, 38)))
#mvrnorm(38, mu = gamma, Sigma = T)

init = c(theta, phi, gamma, beta, T, alpha)






sample = Gibbs_School(N_chain = 200, init)

plot(sample[,'beta_1'], type = 'l', main = 'evolution de beta_1')
plot(sample[,'beta_2'], type = 'l', main = 'evolution de beta_2')
plot(sample[,'beta_3'], type = 'l', main = 'evolution de beta_3')
plot(sample[,'beta_4'], type = 'l', main = 'evolution de beta_4')
plot(sample[,'beta_5'], type = 'l', main = 'evolution de beta_5')
plot(sample[,'beta_6'], type = 'l', main = 'evolution de beta_6')
plot(sample[,'beta_7'], type = 'l', main = 'evolution de beta_7')
plot(sample[,'beta_8'], type = 'l', main = 'evolution de beta_8')

plot(sample[,'gamma_1'], type = 'l', main = 'evolution de gamma_1')
plot(sample[,'gamma_2'], type = 'l', main = 'evolution de gamma_2')
plot(sample[,'gamma_3'], type = 'l', main = 'evolution de gamma_3')

plot(sample[,'phi'], type = 'l', main = 'evolution de phi')
plot(sample[,'theta'], type = 'l', main = 'evolution de theta')

matplot(sample[,c(paste0("alpha_1_", 1:38))], type = 'l', main = 'evolution de alpha_1')
matplot(sample[,c(paste0("alpha_2_", 1:38))], type = 'l', main = 'evolution de alpha_2')
matplot(sample[,c(paste0("alpha_3_", 1:38))], type = 'l', main = 'evolution de alpha_3')
























for(j in 1:M){
  
  #nb_in_school_j = sum(datadf[,'school'] == j)
  eleve_in_school_j = which(datadf[,'school'] == j)
  #premier_eleve_in_school_j = eleve_in_school_j[1]
  
  
  tau_j = exp(theta + phi * LRT[eleve_in_school_j])
  LRT_j = LRT[eleve_in_school_j]
  VR1_j = VR[eleve_in_school_j, 1]
  Y_j = Y[eleve_in_school_j]
  
  
  ## maj de alpha_1j
  
  mu_j_sans_alpha1j = alpha[2,j] * LRT[eleve_in_school_j] + alpha[3,j] * VR[eleve_in_school_j,1]  +
    beta[1] * LRT[eleve_in_school_j]^2           + beta[2] * VR[eleve_in_school_j, 2]            +     
    beta[3] * Gender[eleve_in_school_j]           + beta[4] * School_gender[eleve_in_school_j, 1] +
    beta[5] * School_gender[eleve_in_school_j, 2] + beta[6] * School_denom[eleve_in_school_j, 1]  +
    beta[7] * School_denom[eleve_in_school_j, 2]  + beta[8] * School_denom[eleve_in_school_j, 3]
  
  
  sigma2_star_alpha1j = 1/(T[1,1] + sum(tau_j))
  mu_star_alpha1j = sigma2_star_alpha1j * (T[1,1]*gamma[1] + 
                                             sum(tau_j*(Y_j - mu_j_sans_alpha1j)))
  
  alpha[1,j] = rnorm(1, mean = mu_star_alpha1j, sd = sqrt(sigma2_star_alpha1j))
  
  
  ## maj de alpha_2j
  #  on doit recalculer mu_j avec les nouvelles valeurs de alpha_1j
  
  mu_j_sans_alpha2j = alpha[1,j] + alpha[3,j] * VR[eleve_in_school_j,1] +
    beta[1] * LRT[eleve_in_school_j]^2           + beta[2] * VR[eleve_in_school_j, 2]            +     
    beta[3] * Gender[eleve_in_school_j]           + beta[4] * School_gender[eleve_in_school_j, 1] +
    beta[5] * School_gender[eleve_in_school_j, 2] + beta[6] * School_denom[eleve_in_school_j, 1]  +
    beta[7] * School_denom[eleve_in_school_j, 2]  + beta[8] * School_denom[eleve_in_school_j, 3]
  
  
  sigma2_star_alpha2j = 1/(T[2,2] + sum(tau_j*LRT_j^2))
  mu_star_alpha2j = sigma2_star_alpha2j * (T[2,2]*gamma[2] + 
                                             sum(tau_j*LRT_j*(Y_j - mu_j_sans_alpha2j)))
  
  alpha[2,j] = rnorm(1, mean = mu_star_alpha2j, sd = sqrt(sigma2_star_alpha2j))
  
  
  ## maj de alpha_3j
  #  on doit recalculer mu_j avec les nouvelles valeurs de alpha_1j et de alpha_2j
  
  mu_j_sans_alpha3j = alpha[1,j] + alpha[2,j] * LRT[eleve_in_school_j] + 
    beta[1] * LRT[eleve_in_school_j]^2           + beta[2] * VR[eleve_in_school_j, 2]            +     
    beta[3] * Gender[eleve_in_school_j]           + beta[4] * School_gender[eleve_in_school_j, 1] +
    beta[5] * School_gender[eleve_in_school_j, 2] + beta[6] * School_denom[eleve_in_school_j, 1]  +
    beta[7] * School_denom[eleve_in_school_j, 2]  + beta[8] * School_denom[eleve_in_school_j, 3]
  
  
  sigma2_star_alpha3j = 1/(T[3,3] + sum(tau_j*VR1_j^2))
  mu_star_alpha3j = sigma2_star_alpha3j * (T[3,3]*gamma[3] + 
                                             sum(tau_j*VR1_j*(Y_j - mu_j_sans_alpha3j)))
  
  alpha[3,j] = rnorm(1, mean = mu_star_alpha3j, sd = sqrt(sigma2_star_alpha3j))
  
}





























mu_sans_beta = rep(NaN, N)
for(j in 1:M){
  
  eleve_in_school_j = which(datadf[,'school'] == j)
  nb_in_school_j = sum(datadf[,'school'] == j)
  
  LRT_j = LRT[eleve_in_school_j]
  VR1_j = VR[eleve_in_school_j, 1]
  
  mu_j_sans_beta = alpha[1,j] + alpha[2,j] * LRT_j + alpha[3,j] * VR1_j
  
  mu_sans_beta[eleve_in_school_j] = mu_j_sans_beta
}

somme_sigma_beta = 0
somme_mu_beta = 0
for(k in 1:N){
  somme_sigma_beta = somme_sigma_beta + 
    tau[k] * tcrossprod(c(LRT[k]^2, 
                          VR[k,2], 
                          Gender[k], 
                          School_gender[k,1], 
                          School_gender[k,2], 
                          School_denom[k,1], 
                          School_denom[k,2], 
                          School_denom[k,3]))
  somme_mu_beta = somme_mu_beta + 
    c(LRT[k]^2, 
      VR[k,2], 
      Gender[k], 
      School_gender[k,1], 
      School_gender[k,2], 
      School_denom[k,1], 
      School_denom[k,2], 
      School_denom[k,3]) * (Y[k] - mu_sans_beta[k]) * tau[k]
}


sigma_star_beta = inv(diag(10^-4, ncol=8, nrow=8) + somme_sigma_beta)
mu_star_beta = sigma_star_beta %*% matrix(somme_mu_beta, ncol = 1, nrow = 8)

beta = mvrnorm(mu = mu_star_beta, Sigma = sigma_star_beta)