source('weight_calc.R')
source('weight_calc_EG.R')

rdu = function(n,k) sample(1:k,n,replace=T)

######################

set.seed(7)

d = c(15, 15, 15)
# Block sizes for experimental treatments after burn-in

J = length(d)
m = sum(d)

d0 = c(8, 8, 8)
# Block sizes for control after burn-in


#####

mu =  17.3/3.5  # True mean of control
delta1 = 66.2/3.5 - mu   # Incremental benefit of treatment 1
delta2 = 73.2/3.5 - mu  # Incremental benefit of treatment 2
mu_prior = c(5, 5, 5)   # Prior means

delta = c(delta1, delta2)


### Burn in

r0 = 7
r1 = 8
r2 = 8

r = r1+r2

n = m+r
n0 = sum(d0)+r0

D = cumsum(d)
D = r+c(0, D)

D0 = cumsum(d0)
D0 = r0 + c(0, D0)


### Hypothesis testing

alpha = 0.05
z_alpha = qnorm(1-alpha)


### Auxiliary design allocations

b = rep(0, n)

b[1:r1] = 1
b[(r1+1):r] = 2

b[(r+1):(n-1)] = rdu(m-1, 2)

############# Actual design allocations ################### 

a = rep(0, n)
a[1:r] = b[1:r]

X = c(rnorm(r, mean = mu + delta[a[1:r]]), rep(0,m))

X0 = rnorm(n0, mean = mu)


##### BAR #####

source('BAR.R')

# mu_prior = c(5, 5, 5)   # Prior means

sigma_prior = c(1, 1, 1)   # Prior variances

gamma_pow = 0.5

mu = mu + delta

out = BAR(D, D0, a, X, X0, mu, mu_prior, sigma_prior, gamma_pow)

a = out$a
X = out$X


##### Calculate weights and modified z-statistic #####

b1 = c(b[1:(n-1)],1)
b2 = c(b[1:(n-1)],2)

results1 = weight_calc(a, b1, D, D0, X, X0, r, 1)
results2 = weight_calc(a, b2, D, D0, X, X0, r, 2)

wI1 = results1$wI
wI2 = results2$wI

w01 = results1$w0
w02 = results2$w0

T1m = results1$TI
T2m = results2$TI

results1_EG = wT_calc_EG(a, b1, D, D0, X, X0, r, 1)
results2_EG = wT_calc_EG(a, b2, D, D0, X, X0, r, 2)

wI1_EG = results1_EG$wI
wI2_EG = results2_EG$wI

w01_EG = results1_EG$w0
w02_EG = results2_EG$w0

T1m_EG = results1_EG$TI
T2m_EG = results2_EG$TI

U1m_EG = results1_EG$Utilde
U2m_EG = results2_EG$Utilde


#####  RW test for H_I #####

# I = {1}

n1 = sum(b == 1)
p1m = 1 - pnorm(T1m*(1/n1 + 1/n0)^(-1/2))

# c1 = z_alpha*(1/n1 + 1/n0)^(1/2)
# H_1 = T1m > c1


# I = {2}

n2 = sum(b == 2)
p2m = 1 - pnorm(T2m*(1/n2 + 1/n0)^(-1/2))

# c2 = z_alpha*(1/n2 + 1/n0)^(1/2)
# H_2 = T2m > c2


#####  EG test for H_I #####

# I = {1}

p1mEG = 1 - pnorm(U1m_EG)


# I = {2}

n2 = sum(b == 2)
p2m = 1 - pnorm(U2m_EG)


#### Naive Z-test for H_I #####


# I = {1,2}

TI = mean(X[a %in% c(1,2)]) - mean(X0)
pI = 1 - pnorm(TI*(1/n + 1/n0)^(-1/2))

# I = {1}

n1a = sum(a == 1)
T1n = mean(X[a == 1]) - mean(X0)

p1n = 1- pnorm(T1n*(1/n1a + 1/n0)^(-1/2))

# I = {2}

n2a = sum(a == 2)
T2n = mean(X[a == 2]) - mean(X0)

p2n = 1 - pnorm(T2n*(1/n2a + 1/n0)^(-1/2))

