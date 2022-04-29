####################################################################################
################  Two treatment Bayesian adaptive randomisation   ##################
####################################################################################

BAR = function(D, D0, a, X, X0, mu, mu_prior, sigma_prior, gamma_pow = 0.5) {
  
  J = length(D)-1
  
  # a = c(a, rep(0, D[J+1]-D[1]))
  # X = c(X, rep(0, D[J+1]-D[1]))
  
  for (pn in 1:J) {
    n0 = D0[pn]
    n1 = sum(a[1:D[pn]] == 1)
    n2 = sum(a[1:D[pn]] == 2)
    n3 = sum(a[1:D[pn]] == 3)
    
    post_mean0 = (sigma_prior[1] / (1 + n0 * sigma_prior[1])) * sum(X0[1:D0[pn]]) +
      (n0 / (1 + n0 * sigma_prior[1])) * mu_prior[1]
    
    post_var0 = sigma_prior[1] / (1 + n0 * sigma_prior[1])
    
    post_mean1 = (sigma_prior[2] / (1 + n1 * sigma_prior[2])) * sum(X[a[1:D[pn]] == 1])
    + (n1 / (1 + n1 * sigma_prior[2])) * mu_prior[2]
    
    post_var1 = sigma_prior[2] / (1 + n1 * sigma_prior[2])
    
    post_mean2 = (sigma_prior[3] / (1 + n2 * sigma_prior[3])) * sum(X[a[1:D[pn]] == 2]) + 
      (n2 / (1 + n2 * sigma_prior[3])) * mu_prior[3]
    
    post_var2 = sigma_prior[3] / (1 + n2 * sigma_prior[3])
    
    post_mean3 = (sigma_prior[4] / (1 + n3 * sigma_prior[4])) * sum(X[a[1:D[pn]] == 3]) + 
      (n3 / (1 + n3 * sigma_prior[4])) * mu_prior[4]
    
    post_var3 = sigma_prior[4] / (1 + n3 * sigma_prior[4])
    
    p1 = pnorm(0, mean = post_mean0 - post_mean1,
               sd = sqrt(post_var1 + post_var0)) ^ gamma_pow
    p2 = pnorm(0, mean = post_mean0 - post_mean2,
               sd = sqrt(post_var2 + post_var0)) ^ gamma_pow
    p3 = pnorm(0, mean = post_mean0 - post_mean3,
               sd = sqrt(post_var3 + post_var0)) ^ gamma_pow
    
    pi1 = p1/(p1+p2+p3)
    pi2 = p2/(p1+p2+p3)
    # pi3 = 1 - pi1 - pi2
    
    # print(c(pi1, pi2, pi3))
    
    a[(D[pn] + 1):D[pn + 1]] = rassign2(D[pn + 1] - D[pn], pi1, pi2)
    
    X[(D[pn] + 1):D[pn + 1]] = rnorm(D[pn + 1] - D[pn], 
                                     mean = mu[a[(D[pn] + 1):D[pn + 1]]])
  }
  
  return(list(a = a, X = X))
}