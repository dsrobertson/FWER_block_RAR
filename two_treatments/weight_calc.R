weight_calc = function(a, b, D, D0, X, X0, r, Iset) {
  
  J = length(D)-1
  n0 = length(X0)
  
  indI_a = a %in% Iset
  indI_b = b %in% Iset
  
  m_Ij = cumsum(indI_b)
  
  rI = m_Ij[r]
  nI = m_Ij[n]
  mI = nI-rI
  
  m_Ij = c(0, m_Ij[1:(n-1)])
  m_Ij = nI - m_Ij
  
  m_0j = n0 - c(0, D0)
  r0 = D0[1]
  
  wI = w0 = rep(0, J+1)
  wI[1] = nI; w0[1] = n0
  
  XI = X*indI_a
  TI = sum(XI[1:r]/wI[1]) - sum(X0[1:r0]/w0[1])
  
  for (j in 1:(J-1)) {
    
    lambda1 = m_Ij[D[j]+1]/wI[j] - m_0j[j+1]/w0[j]
    lambda2 = m_Ij[D[j]+1]/(wI[j])^2 + m_0j[j+1]/(w0[j])^2
    
    mtI = m_Ij[D[j+1]+1] + sum(indI_a[(D[j]+1):D[j+1]])
    mt0 = m_0j[j+1]
    
    # mtI = m_Ij[D[j]+1] + sum(indI_a[(D[j]+1):D[j+1]] & !(indI_b[(D[j]+1):D[j+1]])) -
    #   sum(!(indI_a[(D[j]+1):D[j+1]]) & indI_b[(D[j]+1):D[j+1]])
    
    
    # if(abs(lambda1^2 - mt0*lambda2) < 1e-10){
    #   
    #   wI[j+1] = -(mt0+mtI)/(2*lambda1)
    #   
    # }
      
    wI[j+1] = (lambda1*mtI - sqrt(lambda1^2*mtI^2 - (lambda1^2 - mt0*lambda2)*mtI*(mt0+mtI)))/(lambda1^2 - mt0*lambda2)
    
    w0[j+1] = mt0*wI[j+1]/(mtI - lambda1*wI[j+1])
    
    TI = TI + sum(XI[(D[j]+1):D[j+1]]/wI[j+1]) - sum(X0[(D0[j]+1):D0[j+1]]/w0[j+1])
  }
  
  lambda1 = m_Ij[D[J]+1]/wI[J] - m_0j[J+1]/w0[J]
  lambda2 = m_Ij[D[J]+1]/(wI[J])^2 + m_0j[J+1]/(w0[J])^2
  
  mtI = sum(indI_a[(D[J]+1):D[J+1]])
  mt0 = m_0j[J+1]
  
  if (mtI > 0){
    
    # if(abs(lambda1^2 - mt0*lambda2) < 1e-10){
    #   
    #   wI[j+1] = -(mt0+mtI)/(2*lambda1)
    #   
    # } 
    
    wI[J+1] = (lambda1*mtI - sqrt(lambda1^2*mtI^2 - (lambda1^2 - mt0*lambda2)*mtI*(mt0+mtI)))/(lambda1^2 - mt0*lambda2)
    w0[J+1] = mt0*wI[J+1]/(mtI - lambda1*wI[J+1])
    
    TI = TI + sum(XI[(D[J]+1):D[J+1]]/wI[J+1]) - sum(X0[(D0[J]+1):D0[J+1]]/w0[J+1])
    
  } else {
    
    mt01 = mt0 - 1
    mt02 = mt0 - mt01
    
    w0_J1 = (-lambda1*mt01 - sqrt(lambda1^2*mt01^2 - (lambda1^2 - mt02*lambda2)*mt01*mt0))/(lambda1^2 - mt02*lambda2)
    w0_J2 =  -mt02*w0_J1/(lambda1*w0_J1 + mt01)
    
    D0_J1 = D0[J]+mt01
    
    TI = TI - sum(X0[(D0[J]+1):D0_J1]/w0_J1) - sum(X0[(D0_J1+1):D0[J+1]]/w0_J2)
    
    wI = c(wI[1:J], NA)
    w0 = c(w0[1:J], w0_J1, w0_J2)
    
  }
  
 return(list(wI = wI, w0 = w0, TI = TI))
  
}