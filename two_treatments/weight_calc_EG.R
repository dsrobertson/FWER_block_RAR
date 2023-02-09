weight_calc_EG = function(a, b, D, D0, X, X0, r, Iset) {
  
  J = length(D)-1
  n0 = length(X0)
  
  indI_a = a %in% Iset
  indI_b = b %in% Iset
  
  m_Ij = cumsum(indI_b)
  
  rI1 = m_Ij[r]
  nI1 = m_Ij[n]
  
  
  T0 = (rI1/nI1)*mean(X[1:r][indI_b[1:r]])
  
  Ttilde = ntilde = rep(NA, J)
  w = rep(NA, J+1)
  
  w[1] = nI1
  
  
  for (pn in 1:J) {
    
    ntilde[pn] = sum(indI_a[(D[pn]+1):D[pn+1]])
    nI = sum(indI_b[(D[pn]+1):D[pn+1]])
    
    if(pn < J){
      nIplus = sum(indI_b[(D[pn+1]+1):D[J+1]])
    } else {
      nIplus = 0}
    
    w[pn+1] = w[pn]*sqrt((ntilde[pn] + nIplus)/(nI + nIplus))
    
    if(sum(indI_a[(D[pn]+1):D[pn+1]]) == 0){
      Ttilde[pn] = 0
    } else {
      Ttilde[pn] = (ntilde[pn]/w[pn+1])*
        mean(X[(D[pn]+1):D[pn+1]][indI_a[(D[pn]+1):D[pn+1]]])
    }
    
    
  }
  
  w0 = ((rI1/nI1) + sum(ntilde/w[-1]))
  
  Ttilde.star = (T0 + sum(Ttilde) - w0*mean(X0))
  
  std = sqrt((1/nI1) + (1/n0)*((rI1/nI1) + sum(ntilde/w[-1]))^2)
  
  Utilde = Ttilde.star/std
  
  return(list(TI = Ttilde.star, Utilde = Utilde, wI = w, w0 = n0*w0))
  
}

