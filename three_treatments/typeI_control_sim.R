################################################################################
########  Three treatment response adaptive block randomisation trial   #########
###########      with burn in and protected control allocation      ############
################################################################################

source('BAR.R')
source('rassign2.R')
source('wT_calc_EG.R')
source('WT_calc.R')

rdu = function(n,k) sample(1:k,n,replace=TRUE)

################################### INPUTS ###################################
################################################################################

N = 10^5  # Number of trials to simulate

# Block sizes for experimental treatments after burn-in
d = c(40, 40, 40)

J = length(d)
m = sum(d)


# Block sizes for control after burn-in
d0 = c(20, 20, 20)
mu = 0  # True mean of control
alpha = 0.05   # Target type I error rate


delta = c(0.5, 0.5, 0.5)  # True means of treatments


### Burn in

r0 = 5
rvec = c(5,5,5)
r = sum(rvec)

n = m+r
n0 = sum(d0)+r0

D = cumsum(d)
D = r+c(0, D)

D0 = cumsum(d0)
D0 = r0 + c(0, D0)


### Allocation scheme (uncomment out)

# alloc_scheme = 'ER'  # Equal Randomization

# alloc_scheme = 'fixed_block'; p = c(0.7, 0.5, 0.3)   # Block Randomization with probability p

# alloc_scheme = 'BAR'   # Bayesian Adaptive Randomization

alloc_scheme = 'inflator'   # Type I error inflator


################################################################################
################################################################################

### Hypothesis testing

alpha = 0.05
z_alpha = qnorm(1-alpha)

p_123 = p_1n = p_2n = p_3n = p_12n = p_13n = p_23n = rep(0,N)
p_1m = p_2m = p_3m = p_12m = p_13m = p_23m = rep(0,N)
p_1mEG = p_2mEG = p_3mEG = p_12mEG = p_13mEG = p_23mEG = rep(0,N)

### Auxiliary design allocations

b = rep(0, n)
B = rep(0, r)

rr = cumsum(rvec)
rr = c(0, rr)

for(i in 1:length(rvec)){
  b[(rr[i]+1):rr[i+1]] = i
}


##### BAR #####

mu_actual = mu + delta

mu_prior = c(0, 0, 0, 0)   # Prior means

sigma_prior = c(1, 1, 1, 1)   # Prior variances

gamma_pow = 0.5


################ Run simulation ################ 

# set.seed(7)

for (i in 1:N) {
  
  a = rep(0, n)
  a[1:r] = b[1:r]
  
  B = rep(0, r)
  
  for(j in 1:length(rvec)){
    
    B[(rr[j]+1):rr[j+1]] = rnorm(rr[j+1]-rr[j], mean = mu + delta[j])
  }
  
  
  b[(r+1):(n-1)] = rdu(m-1, 3)
  
  X0 = rnorm(n0, mean = mu)
  X = c(B, rep(0,m))
  
  #############################################################################
  ################### Allocation scheme  ########################
  
  if(alloc_scheme == 'FR'){
    
    #### Non-adaptive equal randomisation #####
    
    for(pn in 1:m) {
      
      a[r+pn] = rdu(1,3)
      X[r+pn] = rnorm(1, mean = mu + delta[a[r+pn]])
    }
    
  } else if (alloc_scheme == 'fixed_block') {
    
    for (pn in 1:J) {
      
      a[(D[pn]+1):D[pn+1]] = rbinom(D[pn+1]-D[pn], 2, prob = p[pn])+1
      
      X[(D[pn]+1):D[pn+1]] = rnorm(D[pn+1]-D[pn], mean = mu + delta[a[(D[pn]+1):D[pn+1]]])
      
    }
    
  } else if (alloc_scheme == 'inflator') {
    
    pn = 1
    
    while(pn <= J){
      
      X1 = X[a == 1]
      X0j = X0[1:D0[pn]]
      
      if(mean(X1) - mean(X0j) > 0.5) {
        a[D[pn]+1] = 1
        X[D[pn]+1] = rnorm(1, mean = mu + delta[1])
        a[(D[pn]+2):D[pn+1]] = rdu(D[pn+1]-D[pn]-1, 2)+1
        X[(D[pn]+2):D[pn+1]] = rnorm(D[pn+1]-D[pn]-1, mean = mu + delta[a[(D[pn]+2):D[pn+1]]])
      } else {
        a[D[pn]+1] = 2
        X[D[pn]+1] = rnorm(1, mean = mu + delta[2])
        a[D[pn]+2] = 3
        X[D[pn]+2] = rnorm(1, mean = mu + delta[3])
        a[(D[pn]+3):D[pn+1]] = 1
        X[(D[pn]+3):D[pn+1]] = rnorm(D[pn+1]-D[pn]-2, mean = mu + delta[1])
      }
      
      pn = pn+1
    }
    
  } else if (alloc_scheme == 'BAR'){
    
    
    out = BAR(D, D0, a, X, X0, mu = mu_actual, mu_prior, sigma_prior, gamma_pow)
    
    a = out$a
    X = out$X
    
  }
  
  
  ##### Calculate weights and modified z-statistic #####
  
  T_1mEG = wT_calc_EG(a, b, D, D0, X, X0, r, 1)
  T_2mEG = wT_calc_EG(a, b, D, D0, X, X0, r, 2)
  T_3mEG = wT_calc_EG(a, b, D, D0, X, X0, r, 3)
  
  T_12mEG = wT_calc_EG(a, b, D, D0, X, X0, r, c(1,2))
  T_13mEG = wT_calc_EG(a, b, D, D0, X, X0, r, c(1,3))
  T_23mEG = wT_calc_EG(a, b, D, D0, X, X0, r, c(2,3))
  
  b1 = c(b[1:(n-1)],1)
  b2 = c(b[1:(n-1)],2)
  b3 = c(b[1:(n-1)],3)
  
  T_1m = wT_calc(a, b1, D, D0, X, X0, r, 1)
  T_2m = wT_calc(a, b2, D, D0, X, X0, r, 2)
  T_3m = wT_calc(a, b3, D, D0, X, X0, r, 3)
  
  T_12m = wT_calc(a, b1, D, D0, X, X0, r, c(1,2))
  T_13m = wT_calc(a, b1, D, D0, X, X0, r, c(1,3))
  T_23m = wT_calc(a, b2, D, D0, X, X0, r, c(2,3))
  
  
  # I = {1}
  
  p_1mEG[i] = 1 - pnorm(T_1mEG)
  
  n1b = sum(b1 == 1)
  p_1m[i] = 1 - pnorm(T_1m*(1/n1b + 1/n0)^(-1/2))
  
  
  # I = {2}
  
  p_2mEG[i] = 1 - pnorm(T_2mEG)
  
  n2b = sum(b2 == 2)
  p_2m[i] = 1 - pnorm(T_2m*(1/n2b + 1/n0)^(-1/2))
  
  
  # I = {3}
  p_3mEG[i] = 1 - pnorm(T_3mEG)
  
  n3b = sum(b3 == 3)
  p_3m[i] = 1 - pnorm(T_3m*(1/n3b + 1/n0)^(-1/2))
  
  
  # I = {1,2}
  p_12mEG[i] = 1 - pnorm(T_12mEG)
  
  n12b = sum(b1 %in% c(1,2))
  p_12m[i] = 1 - pnorm(T_12m*(1/n12b + 1/n0)^(-1/2))
  
  
  # I = {1,3}
  p_13mEG[i] = 1 - pnorm(T_13mEG)
  
  n13b = sum(b1 %in% c(1,3))
  p_13m[i] = 1 - pnorm(T_13m*(1/n13b + 1/n0)^(-1/2))
  
  
  # I = {2,3}
  p_23mEG[i] = 1 - pnorm(T_23mEG)
  
  n23b = sum(b2 %in% c(2,3))
  p_23m[i] = 1 - pnorm(T_23m*(1/n23b + 1/n0)^(-1/2))
  
  ##### Naive Z-test for H_I #####
  
  # I = {1,2,3}
  n123 = sum(a %in% c(1,2,3))
  T_123 = mean(X[a != 0]) - mean(X0)
  
  p_123[i] = 1 - pnorm(T_123*(1/n123 + 1/n0)^(-1/2))
  
  # I = {1}
  n1a = sum(a == 1)
  T_1n = mean(X[a == 1]) - mean(X0)
  
  p_1n[i] = 1 - pnorm(T_1n*(1/n1a + 1/n0)^(-1/2))
  
  # I = {2}
  n2a = sum(a == 2)
  T_2n = mean(X[a == 2]) - mean(X0)
  
  p_2n[i] = 1 - pnorm(T_2n*(1/n2a + 1/n0)^(-1/2))
  
  # I = {3}
  n3a = sum(a == 3)
  T_3n = mean(X[a == 3]) - mean(X0)
  
  p_3n[i] = 1 - pnorm(T_3n*(1/n3a + 1/n0)^(-1/2))
  
  # I = {1,2}
  n12a = sum(a %in% c(1,2))
  T_12n = mean(X[a %in% c(1,2)]) - mean(X0)
  
  p_12n[i] = 1 - pnorm(T_12n*(1/n12a + 1/n0)^(-1/2))
  
  # I = {1,3}
  n13a = sum(a %in% c(1,3))
  T_13n = mean(X[a %in% c(1,3)]) - mean(X0)
  
  p_13n[i] = 1 - pnorm(T_13n*(1/n13a + 1/n0)^(-1/2))
  
  # I = {2,3}
  n23a = sum(a %in% c(2,3))
  T_23n = mean(X[a %in% c(2,3)]) - mean(X0)
  
  p_23n[i] = 1 - pnorm(T_23n*(1/n23a + 1/n0)^(-1/2))
  
}

#### 

H_1nc = c(p_1n < alpha & p_12n < alpha & p_13n < alpha & p_123 < alpha)
H_2nc = c(p_2n < alpha & p_12n < alpha & p_23n < alpha & p_123 < alpha)
H_3nc = c(p_3n < alpha & p_13n < alpha & p_23n < alpha & p_123 < alpha)

H_1mc = c(p_1m < alpha & p_12m < alpha & p_13m < alpha & p_123 < alpha)
H_2mc = c(p_2m < alpha & p_12m < alpha & p_23m < alpha & p_123 < alpha)
H_3mc = c(p_3m < alpha & p_13m < alpha & p_23m < alpha & p_123 < alpha)

H_1mcEG = c(p_1mEG < alpha & p_12mEG < alpha & p_13mEG < alpha & p_123 < alpha)
H_2mcEG = c(p_2mEG < alpha & p_12mEG < alpha & p_23mEG < alpha & p_123 < alpha)
H_3mcEG = c(p_3mEG < alpha & p_13mEG < alpha & p_23mEG < alpha & p_123 < alpha)

H_1nb = c(p_1n < alpha/3)
H_2nb = c(p_2n < alpha/3)
H_3nb = c(p_3n < alpha/3)

H_1nbc = c(p_1n < alpha & (p_1n < alpha/2 | p_2n < alpha/2) & (p_1n < alpha/2 | p_3n < alpha/2)
           & (p_1n < alpha/3 | p_2n < alpha/3 | p_3n < alpha/3))

H_2nbc = c(p_2n < alpha & (p_2n < alpha/2 | p_1n < alpha/2) & (p_2n < alpha/2 | p_3n < alpha/2)
           & (p_1n < alpha/3 | p_2n < alpha/3 | p_3n < alpha/3))

H_3nbc = c(p_3n < alpha & (p_3n < alpha/2 | p_1n < alpha/2) & (p_3n < alpha/2 | p_2n < alpha/2)
           & (p_1n < alpha/3 | p_2n < alpha/3 | p_3n < alpha/3))


H_1mbc = c(p_1m < alpha & (p_1m < alpha/2 | p_2m < alpha/2) & (p_1m < alpha/2 | p_3m < alpha/2)
           & (p_1m < alpha/3 | p_2m < alpha/3 | p_3m < alpha/3))

H_2mbc = c(p_2m < alpha & (p_2m < alpha/2 | p_1m < alpha/2) & (p_2m < alpha/2 | p_3m < alpha/2)
           & (p_1m < alpha/3 | p_2m < alpha/3 | p_3m < alpha/3))

H_3mbc = c(p_3m < alpha & (p_3m < alpha/2 | p_1m < alpha/2) & (p_3m < alpha/2 | p_2m < alpha/2)
           & (p_1m < alpha/3 | p_2m < alpha/3 | p_3m < alpha/3))

H_1mbcEG = c(p_1mEG < alpha & (p_1mEG < alpha/2 | p_2mEG < alpha/2) & 
               (p_1mEG < alpha/2 | p_3mEG < alpha/2)
             & (p_1mEG < alpha/3 | p_2mEG < alpha/3 | p_3mEG < alpha/3))

H_2mbcEG = c(p_2mEG < alpha & (p_2mEG < alpha/2 | p_1mEG < alpha/2) & 
               (p_2mEG < alpha/2 | p_3mEG < alpha/2)
             & (p_1mEG < alpha/3 | p_2mEG < alpha/3 | p_3mEG < alpha/3))

H_3mbcEG = c(p_3mEG < alpha & (p_3mEG < alpha/2 | p_1mEG < alpha/2) & 
               (p_3mEG < alpha/2 | p_2mEG < alpha/2)
             & (p_1mEG < alpha/3 | p_2mEG < alpha/3 | p_3mEG < alpha/3))



##############################

print('Probability to reject')

if(sum(delta) == 0 | prod(delta) > 0) {
  
  print(sum(c(H_1mc == TRUE | H_2mc == TRUE | H_3mc == TRUE))/N)
  print(sum(c(H_1mcEG == TRUE | H_2mcEG == TRUE | H_3mcEG == TRUE), na.rm = TRUE)/N)
  print(sum(c(H_1mbc == TRUE | H_2mbc == TRUE | H_3mbc == TRUE))/N)
  print(sum(c(H_1mbcEG == TRUE | H_2mbcEG == TRUE | H_3mbcEG == TRUE), na.rm = TRUE)/N)
  print(sum(c(H_1nc == TRUE | H_2nc == TRUE | H_3nc == TRUE))/N)
  print(sum(c(H_1nbc == TRUE | H_2nbc == TRUE | H_3nbc == TRUE))/N)
  print(sum(c(H_1nb == TRUE | H_2nb == TRUE | H_3nb == TRUE))/N)
  
} else if(sum(delta[1:2]) == 0) {
  
  print(c(sum(H_1mc == TRUE | H_2mc == TRUE), sum(H_3mc == TRUE))/N)
  print(c(sum(H_1mcEG == TRUE | H_2mcEG == TRUE, na.rm = TRUE), sum(H_3mcEG == TRUE, na.rm = TRUE))/N)
  print(c(sum(H_1mbc == TRUE | H_2mbc == TRUE), sum(H_3mbc == TRUE))/N)
  print(c(sum(H_1mbcEG == TRUE | H_2mbcEG == TRUE, na.rm = TRUE), sum(H_3mbcEG == TRUE, na.rm = TRUE))/N)
  print(c(sum(H_1nc == TRUE | H_2nc == TRUE), sum(H_3nc == TRUE))/N)
  print(c(sum(H_1nbc == TRUE | H_2nbc == TRUE), sum(H_3nbc == TRUE))/N)
  print(c(sum(H_1nb == TRUE | H_2nb == TRUE), sum(H_3nb == TRUE))/N)
  
} else {
  
  print(c(sum(H_1mc == TRUE), sum(H_2mc == TRUE | H_3mc == TRUE))/N)
  print(c(sum(H_1mcEG == TRUE, na.rm = TRUE), sum(H_2mcEG == TRUE | H_3mcEG == TRUE, na.rm = TRUE))/N)
  print(c(sum(H_1mbc == TRUE), sum(H_2mbc == TRUE | H_3mbc == TRUE))/N)
  print(c(sum(H_1mbcEG == TRUE, na.rm = TRUE), sum(H_2mbcEG == TRUE | H_3mbcEG == TRUE, na.rm = TRUE))/N)
  print(c(sum(H_1nc == TRUE), sum(H_2nc == TRUE | H_3nc == TRUE))/N)
  print(c(sum(H_1nbc == TRUE), sum(H_2nbc == TRUE | H_3nbc == TRUE))/N)
  print(c(sum(H_1nb == TRUE), sum(H_2nb == TRUE | H_3nb == TRUE))/N)
}


print('Non-defined test statistics')

print(c(sum(is.na(p_1mEG)), sum(is.na(p_2mEG)), sum(is.na(p_3mEG))))
print(c(sum(is.na(p_12mEG)), sum(is.na(p_13mEG)), sum(is.na(p_23mEG))))
