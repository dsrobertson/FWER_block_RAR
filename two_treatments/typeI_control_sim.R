################################################################################
########  Two treatment response adaptive block randomisation trial   ##########
###########      with burn in and protected control allocation      ############
################################################################################

source('BAR.R')
source('wT_calc_EG.R')
source('wT_calc.R')

################################### INPUTS ###################################
################################################################################

N = 10^5  # Number of trials to simulate

# Block sizes for experimental treatments after burn-in
d = c(40, 40, 40)

# Block sizes for control after burn-in
d0 = c(20, 20, 20)


mu = 0  # True mean of control

delta1 = 0.5   # Incremental benefit of treatment 1
delta2 = 0.5   # Incremental benefit of treatment 2


alpha = 0.05   # Target type I error rate


### Burn in

r0 = 5   # Control
r1 = 5   # Treatment 1 
r2 = 5   # Treatment 2


### Allocation scheme (uncomment out)

# alloc_scheme = 'FR'; p = 0.8   # Fixed Randomization with probability p

# alloc_scheme = 'BAR'   # Bayesian Adaptive Randomization

alloc_scheme = 'inflator'   # Type I error inflator


################################################################################
################################################################################

J = length(d)
m = sum(d)
delta = c(delta1, delta2)

r = r1+r2

n = m+r
n0 = sum(d0)+r0

D = cumsum(d)
D = r+c(0, D)

D0 = cumsum(d0)
D0 = r0 + c(0, D0)


### Hypothesis testing
z_alpha = qnorm(1-alpha)
z_alphaB = qnorm(1-alpha/2)

H_12n = H_1n = H_2n = H_1m = H_2m = H_12mEG = H_1mEG = H_2mEG = rep(0,N)

H_1nb = H_2nb = H_1mb = H_2mb = H_1mbEG = H_2mbEG = rep(0,N)

T1mEG_vec = T1m_vec = T2mEG_vec = T2m_vec = rep(0,N)

### Auxiliary design allocations

b = rep(0, n)

b[1:r1] = 1
b[(r1+1):r] = 2


##### BAR #####

mu_actual = mu + delta

mu_prior = c(0, 0, 0)   # Prior means

sigma_prior = c(1, 1, 1)   # Prior variances

gamma_pow = 0.5


################ Run simulation ################ 

set.seed(7)

for (i in 1:N) {
  
  a = rep(0, n)
  a[1:r] = b[1:r]
  
  B1 = rnorm(r1, mean = mu + delta1)
  B2 = rnorm(r2, mean = mu + delta2)
  
  X0 = rnorm(n0, mean = mu)
  
  X = c(B1, B2, rep(0,m))
  
  b[(r+1):n] = rbinom(m, 1, prob = 0.5)+1   # 'Pre-plan' allocation
  
  #############################################################################
  ################### Allocation scheme  ########################
  
  
  if(alloc_scheme == 'FR'){
    
    a[(r+1):n] = rbinom(m, 1, prob = p)+1
    X[(r+1):n] = rnorm(m, mean = mu+delta[a[(r+1):n]])
    
    
  } else if (alloc_scheme == 'inflator') {
    
    pn = 1

    while(pn <= J){

      X1 = X[a == 1]
      X0j = X0[1:D0[pn]]

      if(mean(X1) - mean(X0j) > 0.5) {
        a[D[pn]+1] = 1
        X[D[pn]+1] = rnorm(1, mean = mu + delta1)
        a[(D[pn]+2):D[pn+1]] = 2
        X[(D[pn]+2):D[pn+1]] = rnorm(D[pn+1]-D[pn]-1, mean = mu + delta2)
      } else {
        a[D[pn]+1] = 2
        X[D[pn]+1] = rnorm(1, mean = mu + delta2)
        a[(D[pn]+2):D[pn+1]] = 1
        X[(D[pn]+2):D[pn+1]] = rnorm(D[pn+1]-D[pn]-1, mean = mu + delta1)
      }

      pn = pn+1
    }
  } else if (alloc_scheme == 'BAR'){


    out = BAR(D, D0, a, X, X0, mu = mu_actual, mu_prior, sigma_prior, gamma_pow)
    
    a = out$a
    X = out$X
    
  }
  
  
  ##### Calculate weights and modified z-statistic #####
  
  T1mEG = wT_calc_EG(a, b, D, D0, X, X0, r, 1)
  T2mEG = wT_calc_EG(a, b, D, D0, X, X0, r, 2)
  T12mEG = wT_calc_EG(a, b, D, D0, X, X0, r, c(1,2))
  
  b1 = c(b[1:(n-1)],1)
  b2 = c(b[1:(n-1)],2)
  
  T1m = wT_calc(a, b1, D, D0, X, X0, r, 1)
  T2m = wT_calc(a, b2, D, D0, X, X0, r, 2)
  
  
  n1 = sum(b == 1)
  n2 = sum(b == 2)
  
  T1mEG_vec[i] = T1mEG
  T1m_vec[i] = T1m/(1/n1 + 1/n0)^(1/2)
  
  T2mEG_vec[i] = T2mEG
  T2m_vec[i] = T2m/(1/n2 + 1/n0)^(1/2)
  
  
  ##### Modified Z-test for H_I #####
  
  # I = {1}
  
  H_1mEG[i] = T1mEG > z_alpha
  H_1mbEG[i] = T1mEG > z_alphaB
  
  
  H_1m[i] = T1m > z_alpha*(1/n1 + 1/n0)^(1/2)
  H_1mb[i] = T1m > z_alphaB*(1/n1 + 1/n0)^(1/2)
  
  
  # I = {2}
  
  H_2mEG[i] = T2mEG > z_alpha
  H_2mbEG[i] = T2mEG > z_alphaB
  
  H_2m[i] = T2m > z_alpha*(1/n2 + 1/n0)^(1/2)
  H_2mb[i] = T2m > z_alphaB*(1/n2 + 1/n0)^(1/2)
  
  
  # I = {1,2}
  
  H_12mEG[i] = T12mEG > z_alpha
  
  
  #### Naive Z-test for H_I #####
  
  # I = {1,2}
  
  T12 = mean(X) - mean(X0)
  H_12n[i] = T12 > z_alpha*(1/n + 1/n0)^(1/2)
  
  # I = {1}
  
  n1a = sum(a == 1)
  T1 = mean(X[a == 1]) - mean(X0)
  
  H_1n[i] = T1 > z_alpha*(1/n1a + 1/n0)^(1/2)
  
  # I = {2}
  
  n2a = sum(a == 2)
  T2 = mean(X[a == 2]) - mean(X0)
  
  H_2n[i] = T2 > z_alpha*(1/n2a + 1/n0)^(1/2)
  
  
  ### Naive Z-test with Bonferonni
  
  H_1nb[i] = T1 > z_alphaB*(1/n1a + 1/n0)^(1/2)
  H_2nb[i] = T2 > z_alphaB*(1/n2a + 1/n0)^(1/2)

}

#### 

H_1nc = c(H_1n == TRUE & H_12n == TRUE)
H_2nc = c(H_2n == TRUE & H_12n == TRUE)

H_1mc = c(H_1m == TRUE & H_12n == TRUE)
H_2mc = c(H_2m == TRUE & H_12n == TRUE)

H_1mcEG = c(H_1mEG == TRUE & H_12mEG == TRUE)
H_2mcEG = c(H_2mEG == TRUE & H_12mEG == TRUE)


H_1nbc = c(H_1n == TRUE & (H_1nb == TRUE | H_2nb == TRUE))
H_2nbc = c(H_2n == TRUE & (H_1nb == TRUE | H_2nb == TRUE))

H_1mbc = c(H_1m == TRUE & (H_1mb == TRUE | H_2mb == TRUE))
H_2mbc = c(H_2m == TRUE & (H_1mb == TRUE | H_2mb == TRUE))

H_1mbcEG = c(H_1mEG == TRUE & (H_1mbEG == TRUE | H_2mbEG == TRUE))
H_2mbcEG = c(H_2mEG == TRUE & (H_1mbEG == TRUE | H_2mbEG == TRUE))


##############################

print('Probability to reject')

if((delta1 == 0 & delta2 == 0) | (delta1 > 0 & delta2 > 0)) {

  print(sum(c(H_1mc == TRUE | H_2mc == TRUE))/N)

  print(sum(c(H_1mcEG == TRUE | H_2mcEG == TRUE), na.rm = TRUE)/N)
  
  print(sum(c(H_1mbc == TRUE | H_2mbc == TRUE))/N)
  
  print(sum(c(H_1mbcEG == TRUE | H_2mbcEG == TRUE), na.rm = TRUE)/N)

  print(sum(c(H_1nc == TRUE | H_2nc == TRUE))/N)

  print(sum(c(H_1nbc == TRUE | H_2nbc == TRUE))/N)

  print(sum(c(H_1nb == TRUE | H_2nb == TRUE))/N)

} else {

  print(c(sum(H_1mc == TRUE)/N, sum(H_2mc == TRUE)/N))
  
  print(c(sum(H_1mcEG == TRUE, na.rm = TRUE)/N, sum(H_2mcEG == TRUE, na.rm = TRUE)/N))

  print(c(sum(H_1mbc == TRUE)/N, sum(H_2mbc == TRUE)/N))
  
  print(c(sum(H_1mbcEG == TRUE, na.rm = TRUE)/N, sum(H_2mbcEG == TRUE, na.rm = TRUE)/N))

  print(c(sum(H_1nc == TRUE)/N, sum(H_2nc == TRUE)/N))

  print(c(sum(H_1nbc == TRUE)/N, sum(H_2nbc == TRUE)/N))

  print(c(sum(H_1nb == TRUE)/N, sum(H_2nb == TRUE)/N))

}


print('Non-defined test statistics')

print(c(sum(is.na(H_1mEG)), sum(is.na(H_2mEG)), sum(is.na(H_12mEG))))


print('SD of test statistics')

print(c(sd(T1m_vec), sd(T1mEG_vec, na.rm = T), sd(T2m_vec), sd(T2mEG_vec, na.rm = T)))
