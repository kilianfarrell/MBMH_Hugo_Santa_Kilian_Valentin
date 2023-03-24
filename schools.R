LRT = c(0.5, 5, 15,5, 7, 12)
Y = c(5, 6, 7,15, 27, 2)
mu = c(15, 6, 17,5, 27, 12)
phi = 0.01
theta = 0.58
post_theta = function(nchain, prop_sd){
  S = length(LRT)
  theta = rep(NA, nchain)
  theta[1] = 0
  for (i in 2:nchain){
    
    current = theta[i-1]
    
    prop = current + rnorm(1, 0, prop_sd)
    
    top = -0.5*(-(prop*S)+(prop*sqrt(0.0001))^2 + (exp(prop)*sum(exp(phi*LRT)%*%(Y - mu)^2)))
    bottom = -0.5*(-(current*S)+(current*sqrt(0.0001))^2 + (exp(current)*sum(exp(phi*LRT)%*%(Y - mu)^2)))
    
    if (runif(1) < exp(top-bottom)){
      theta[i] = prop
    }
    else theta[i] = current
    
  }
  return(theta)
}
th = post_theta(10000, 0.1)
plot(th, type='l')

post_phi = function(nchain, init, prop_sd){
  S = length(LRT)
  phi = rep(NA, nchain)
  phi[1] = init
  for (i in 2:nchain){
    
    current = phi[i-1]
    
    prop = current + rnorm(1, 0, prop_sd)
    
    top = -0.5*(-(prop*sum(LRT))+(prop*sqrt(0.0001))^2 + (exp(theta)*sum(exp(prop*LRT)%*%(Y - mu)^2)))
    bottom = -0.5*(-(current*sum(LRT))+(current*sqrt(0.0001))^2 + (exp(theta)*sum(exp(current*LRT)%*%(Y - mu)^2)))
    
    if (runif(1) < exp(top-bottom)){
      phi[i] = prop
    }
    else phi[i] = current
    
  }
  return(phi)
}

ph = post_phi(1000, 0.01, 0.1)
plot(ph, type = 'l')
#####################################################
init = list ("theta" = 0,
             "phi" = 0.01,
             "gamma" = c(0, 0, 0),
             "beta" = c(0, 0, 0, 0, 0, 0, 0, 0),
             "Omega" = structure(c(10, -3, -3, -3, 135, -65, -3, -65, 135), .Dim = as.integer(c(3,
                                                                                        3))),
             "alpha" = structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), .Dim = as.integer(c(38,
                                                                                       3))))

gibbs_schools = function(nchain, Y, init){
  
  #inialisation
  n = 1978
  M = 38
  
  alpha_1 = matrix(NA, nrow = nchain, ncol = M)
  alpha_2 = matrix(NA, nrow = nchain, ncol = M)
  alpha_3 = matrix(NA, nrow = nchain, ncol = M)
  beta = matrix(NA, nrow = nchain, ncol = 8)
  theta = rep(NA, nchain)
  phi = rep(NA, nchain)
  tau = rep(NA, nchain)
  mu = rep(NA, nchain)
  Omega = rep(NA, nchain)
  gamma = matrix(NA, nrow = nchain, ncol = 3)
  
  
  alpha_1[1,] = init$alpha[,1]
  alpha_2[1,] = init$alpha[,2]
  alpha_3[1,] = init$alpha[,3]
  beta[1,] = init$beta
  theta[1] = init$theta
  phi[1] = init$phi
  tau[1] = exp(theta[1]+phi[1]*LRT[1])
  mu[1] = alpha_1[1,] + alpha_2[1,]*LRT[1] + alpha_1[1,]*VR[1,1] + beta[1,1]*LRT[1]^2 + 
    beta[1,2]*VR[1,2] + beta[1,3]*Gender[1] + beta[1,4]*School_gender[1,1] + 
    beta[1,5]*School_gender[1,2] + beta[1,6]*School_denom[1,1] + beta[1,7]*School_denom[1,2] +
    beta[1,8]*School_denom[1,3]
  Omega = init$Omega
  gamma[1,] = init$gamma
  
  for (i in 2:nchain){
    
    #update gamma
    mu_gamma = 
    s_gamma = M * inv(Omega) + diag()
    gamma[i, ] = rmvnorm()
    
    
    
    
    #update alpha_1
    for (j in 1:M){
      
    }
    sd_star = 1/(inv(Omega[1,1]) + sum(tau)
    alpha_1[i,] = rnorm(1, mean = mu_star, sd=sd_star)
    
  }
  
  
  
  
}
