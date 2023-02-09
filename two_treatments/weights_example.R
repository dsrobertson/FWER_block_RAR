################################################################################
########  Two treatment response adaptive block randomisation trial   ##########
###########      with burn in and protected control allocation      ############
################################################################################

source('BAR.R')
source('weight_calc_EG.R')
source('weight_calc.R')

################################### INPUTS ###################################
################################################################################

# Block sizes for experimental treatments after burn-in
d = c(40, 40, 40)

# Block sizes for control after burn-in
d0 = c(20, 20, 20)


mu = 0  # True mean of control

delta1 = 0  # Incremental benefit of treatment 1
delta2 = 1  # Incremental benefit of treatment 2


alpha = 0.05   # Target type I error rate


### Burn in

r0 = 5   # Control
r1 = 5   # Treatment 1 
r2 = 5   # Treatment 2


### Allocation scheme (uncomment out)

alloc_scheme = 'BAR'   # Bayesian Adaptive Randomization

# alloc_scheme = 'inflator'   # Type I error inflator


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

a = rep(0, n)
a[1:r] = b[1:r]

B1 = rnorm(r1, mean = mu + delta1)
B2 = rnorm(r2, mean = mu + delta2)

X0 = rnorm(n0, mean = mu)

X = c(B1, B2, rep(0,m))

b[(r+1):n] = rbinom(m, 1, prob = 0.5)+1   # 'Pre-plan' allocation

#############################################################################
################### Allocation scheme  ########################


if(alloc_scheme == 'inflator') {
  
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

T1mEG = weight_calc_EG(a, b, D, D0, X, X0, r, 1)
T2mEG = weight_calc_EG(a, b, D, D0, X, X0, r, 2)
T12mEG = weight_calc_EG(a, b, D, D0, X, X0, r, c(1,2))

b1 = c(b[1:(n-1)],1)
b2 = c(b[1:(n-1)],2)

T1m = weight_calc(a, b1, D, D0, X, X0, r, 1)
T2m = weight_calc(a, b2, D, D0, X, X0, r, 2)



# Treatment 1

n1 = sum(a == 1)   # Number of patients on treatment 1

T1m$TI/((1/n1 + 1/n0)^(1/2))   # RW test statistic
T1mEG$Utilde   # New test statistic


T1m$wI   # RW weights for experimental treatment 1
T1mEG$wI[2:4]   # New weights for experimental treatment 1

T1m$w0   # RW weights for control
T1mEG$w0   # New weights for control



# Treatment 2

n2 = sum(a == 2)   # # Number of patients on treatment 2

T2m$TI/((1/n2 + 1/n0)^(1/2))   # RW test statistic
T2mEG$Utilde   # New test statistic


T2m$wI   # RW weights for experimental treatment 2
T2mEG$wI[2:4]


T2m$w0   # RW weights for control
T2mEG$w0   # New weights for control
