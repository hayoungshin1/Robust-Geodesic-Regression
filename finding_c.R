## finding c for Huber,Tukey, and L_1 efficiency

Q <- function(s,x) {
  return (Rgamma(s,x,lower=F))
}

Qinv <- function(s,x) {
  return (Rgamma.inv(s,x,lower=F))
}

P <- function(s,x) {
  return (Rgamma(s,x,lower=T))
}

Pinv <- function(s,x) {
  return (Rgamma.inv(s,x,lower=T))
}

uppergamma <- function(s,x) {
  return (Igamma(s,x,lower=F))
}

lowergamma <- function(s,x) {
  return (Igamma(s,x,lower=T))
}

efficiency <- function (c,k,m_estimator) {
  if (m_estimator == 'huber') {
    numfactor <- (k/2)*lowergamma(k/2,0.5*c^2)+c*(k-1)*(2^-1.5)*uppergamma((k-1)/2,0.5*c^2)
    denfactor <- lowergamma((k+2)/2,0.5*c^2)+(0.5*c^2)*uppergamma(k/2,0.5*c^2)
  }
  if (m_estimator == 'tukey') {
    numfactor <- (2*(k+4)/(c^4))*lowergamma((k+4)/2,0.5*c^2)-(2*(k+2)/(c^2))*lowergamma((k+2)/2,0.5*c^2)+(k/2)*lowergamma(k/2,0.5*c^2)
    denfactor <- lowergamma((k+2)/2,0.5*c^2)-(8/(c^2))*lowergamma((k+4)/2,0.5*c^2)+(24/(c^4))*lowergamma((k+6)/2,0.5*c^2)-(32/(c^6))*lowergamma((k+8)/2,0.5*c^2)+(16/(c^8))*lowergamma((k+10)/2,0.5*c^2)
  }
  result <- (numfactor^2)/(gamma((k+2)/2)*denfactor)
  return (result-0.95)
}

l1_efficiency <- function(k) {
  return (((gamma((k+1)/2))^2)/((gamma(k/2)*gamma((k+2)/2))))
}

deriv <- function(c,k,m_estimator) {
  if (m_estimator == 'huber') {
    factor1 <- (k/2)*lowergamma(k/2,0.5*c^2)+c*(k-1)*(2^-1.5)*uppergamma((k-1)/2,0.5*c^2)
    factor2 <- lowergamma((k+2)/2,0.5*c^2)+(0.5*c^2)*uppergamma(k/2,0.5*c^2)
    factor3 <- (c^(k-1))*(2^(-k/2))*exp(-0.5*c^2)+(k-1)*(2^-1.5)*uppergamma((k-1)/2,0.5*c^2)
    factor4 <- c*uppergamma(k/2,0.5*c^2)
  }
  if (m_estimator == 'tukey') {
    factor1 <- (2*(k+4)/(c^4))*lowergamma((k+4)/2,0.5*c^2)-(2*(k+2)/(c^2))*lowergamma((k+2)/2,0.5*c^2)+(k/2)*lowergamma(k/2,0.5*c^2)
    factor2 <- lowergamma((k+2)/2,0.5*c^2)-(8/(c^2))*lowergamma((k+4)/2,0.5*c^2)+(24/(c^4))*lowergamma((k+6)/2,0.5*c^2)-(32/(c^6))*lowergamma((k+8)/2,0.5*c^2)+(16/(c^8))*lowergamma((k+10)/2,0.5*c^2)
    factor3 <- -(8*(k+4)/(c^5))*lowergamma((k+4)/2,0.5*c^2)+(4*(k+2)/(c^3))*lowergamma((k+2)/2,0.5*c^2)-c^(k-1)*2^(-((k-2)/2))*exp(-0.5*c^2)
    factor4 <- (16/(c^3))*lowergamma((k+4)/2,0.5*c^2)-(96/(c^5))*lowergamma((k+6)/2,0.5*c^2)+(192/(c^7))*lowergamma((k+8)/2,0.5*c^2)-(128/(c^9))*lowergamma((k+10)/2,0.5*c^2)
  }
  numerator <- 2*factor1*factor3*factor2-(factor1^2)*factor4
  denominator <- (gamma((k+2)/2))*factor2^2
  result <- numerator/denominator
  return (result)
}

nr <- function(c,k,m_estimator) { #Newton Raphson
  old_c <- c+5
  new_c <- c
  count <- 0
  while (abs(new_c-old_c)>0.000001) {
    old_c <- new_c
    new_c <- old_c-efficiency(old_c,k,m_estimator)/deriv(old_c,k,m_estimator)
    count <- count+1
    if (count > 1000) {
      return ('fail')
    }
  }
  if (m_estimator == 'tukey') {
    new_c <- abs(new_c)
  }
  return (new_c)
}
