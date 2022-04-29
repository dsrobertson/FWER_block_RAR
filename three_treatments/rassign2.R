rassign2 = function(n, p1, p2){
  
  out = rep(0,n)
  ru = runif(n)
  
  for(i in 1:n) {
    
    rui = ru[i]
    
    if(rui < p1) {
      out[i] = 1
    } else if(rui < p1+p2) {
      out[i] = 2
    } else {
      out[i] = 3
    }}
  
  return(out)
}