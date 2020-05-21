## Kendall's 2d shape space

norm <- function(v) {abs(sqrt(sum(v*Conj(v))))}

expo <- function(p,v) {
  p <- (p-mean(p))/norm(p-mean(p))
  theta <- norm(v)
  if (theta==0) {
    result <- p
  } else {
    e1 <- p
    e2 <- v/theta
    result <- cos(theta)*e1 + sin(theta)*e2
  }
  return (result)
}

loga <- function(p1,p2) {
  p1 <- (p1-mean(p1))/norm(p1-mean(p1))
  p2 <- (p2-mean(p2))/norm(p2-mean(p2))
  if (norm(p1-p2) == 0) {
    result <- integer(length(p1))
  } else {
    a <- sum(p1*Conj(p2))
    theta <- acos(max(min(abs(a),1),-1))
    tang <- (a/abs(a))*p2 - abs(a)*p1
    if (norm(tang) == 0) {
      if ((p1[1] == p2[1])) {
        result <- integer(length(p1))
      }
    } else {
      result <- theta*(tang/norm(tang))
    }
  }
  return (result)
}

dist <- function(p1,p2) {
  return (norm(loga(p1,p2)))
}

gammaprime <- function(t,p,v) {
  return (loga(expo(p,t*v),expo(p,(t+1)*v)))
}

pt <- function(p1,p2,v) {
  p1 <- (p1-mean(p1))/norm(p1-mean(p1))
  p2 <- (p2-mean(p2))/norm(p2-mean(p2))
  yi <- expo(p1,v)
  a <- sum(p1*Conj(p2))
  p2 <- (a/abs(a))*p2 #p2 is now p2star
  if (abs(a)>=1) {
    result <- loga(p2,yi)
  } else {
    b <- (1-(abs(a))^2)^0.5
    p2tilda <- (p2-abs(a)*p1)/b
    result <- v-(sum(v*Conj(p1)))*p1-(sum(v*Conj(p2tilda)))*p2tilda+((abs(a))*(sum(v*Conj(p1)))-b*(sum(v*Conj(p2tilda))))*p1+(b*(sum(v*Conj(p1)))+(abs(a))*(sum(v*Conj(p2tilda))))*p2tilda
    result <- (Conj(a/abs(a)))*result
  }
  return (result)
}

rho <- function(x,m_estimator) {
  if (m_estimator[[1]] == 'l2') {
    result <- 0.5*x^2
  } else if (m_estimator[[1]] == 'l1') {
    result <- abs(x)
  } else if (m_estimator[[1]] == 'huber') {
    if (abs(x) < m_estimator[[2]]) {
      result <- 0.5*x^2
    } else {
      result <- m_estimator[[2]]*abs(x)-0.5*m_estimator[[2]]^2
    }
  } else if (m_estimator[[1]] == 'tukey') {
    if (abs(x) < m_estimator[[2]]) {
      result <- ((m_estimator[[2]]^2)/6)*(1-(1-(x/m_estimator[[2]])^2)^3)
    } else {
      result <- (m_estimator[[2]]^2)/6
    }
  }
  return (result)
}

rho_prime <- function(x,m_estimator) {
  if (m_estimator[[1]] == 'l2') {
    result <- x
  } else if (m_estimator[[1]] == 'l1') {
    result <- sign(x)
  } else if (m_estimator[[1]] == 'huber') {
    if (abs(x) < m_estimator[[2]]) {
      result <- x
    } else {
      result <- m_estimator[[2]]*sign(x)
    }
  } else if (m_estimator[[1]] == 'tukey') {
    if (abs(x) < m_estimator[[2]]) {
      result <- x*((1-(x/m_estimator[[2]])^2)^2)*sign(x)
    } else {
      result <- 0
    }
  }
  return (result)
}

eps <- function(p,v,x,y) {
  answer <- matrix(,nrow=length(p),ncol=dim(y)[2])
  shifts <- v%*%t(x)
  for (i in 1:dim(y)[2]) {
    answer[,i] <-loga(expo(p, shifts[,i]), y[,i])
  }
  return (answer)
}

loss_sum <- function(p,v,x,y,m_estimator) {
  sum <- 0
  res <- eps(p,v,x,y)
  for (i in 1:dim(y)[2]) {
    sum <- sum + rho(norm(res[,i]),m_estimator)
  }
  return (sum)
}

j_p <- function(p,v1,v2) {
  if (norm(v1) != 0) {
    j <- (0+1i)*v1
    v2_0 <- pt(expo(p,v1),p,v2)
    w_0 <- (Re(sum(v2_0*Conj(j/norm(j)))))*(j/norm(j))
    u_0 <- v2_0-w_0
    w_tan <- (Re(sum(w_0*Conj(v1/norm(v1)))))*(v1/norm(v1))
    w_orth <- w_0-w_tan
    u_tan <- (Re(sum(u_0*Conj(v1/norm(v1)))))*(v1/norm(v1))
    u_orth <- u_0-u_tan
    L <- norm(v1)
    result <- cos(L)*u_orth + cos(2*L)*w_orth + u_tan + w_tan
  } else {
    result <- v2
  }
  return (result)
}

j_v <- function(p,v1,v2) { 
  if (norm(v1) != 0) {
    j <- (0+1i)*v1
    v2_0 <- pt(expo(p,v1),p,v2)
    w_0 <- (Re(sum(v2_0*Conj(j/norm(j)))))*(j/norm(j))
    u_0 <- v2_0-w_0
    w_tan <- (Re(sum(w_0*Conj(v1/norm(v1)))))*(v1/norm(v1))
    w_orth <- w_0-w_tan
    u_tan <- (Re(sum(u_0*Conj(v1/norm(v1)))))*(v1/norm(v1))
    u_orth <- u_0-u_tan
    L <- norm(v1)
    result <- ((sin(L))/L)*u_orth + ((sin(2*L))/(2*L))*w_orth + u_tan + w_tan
  } else {
    result <- v2
  }
  return (result)
}

grad_p <- function(p,v,x,y,m_estimator) {
  sum <- integer(length(p))
  res <- eps(p,v,x,y)
  if (dim(x)[2]==1) {
    for (i in 1:dim(y)[2]) {
      if (norm(res[,i])!=0) {
        sum <- sum - rho_prime(norm(res[,i]),m_estimator)*j_p(p,x[i]*as.vector(v),(res[,i]/norm(res[,i])))
      }
    }    
  } else {
    shifts <- v%*%t(x)
    for (i in 1:dim(y)[2]) {
      if (norm(res[,i])!=0) {
        sum <- sum - rho_prime(norm(res[,i]),m_estimator)*pt(expo(p, shifts[,i]),p,(res[,i]/norm(res[,i])))
      }
    }
  }
  return (sum)
}

grad_v <- function(p,v,x,y,m_estimator) {
  sum <- matrix(0L,nrow=length(p),ncol=dim(x)[2])
  res <- eps(p,v,x,y)
  if (dim(x)[2]==1) {
    for (i in 1:dim(y)[2]) {
      if (norm(res[,i])!=0) {
        sum <- sum - x[i]*rho_prime(norm(res[,i]),m_estimator)*j_v(p,x[i]*as.vector(v),(res[,i]/norm(res[,i])))
      }
    }    
  } else {
    shifts <- v%*%t(x)
    for (h in 1:dim(x)[2]) {
      for (i in 1:dim(y)[2]) {
        if (norm(res[,i])!=0) {
          sum[,h] <- sum[,h] - x[i,h]*rho_prime(norm(res[,i]),m_estimator)*pt(expo(p, shifts[,i]),p,(res[,i]/norm(res[,i])))
        }
      }
    }
  }
  return (sum)
}

alg <- function(p,v,x,y,m_estimator) {
  p <- (p-mean(p))/norm(p-mean(p))
  current_p <- p
  current_v <- v
  old_p <- integer(length(p))
  old_p[1] <- 0.5^0.5
  old_p[embed] <- -(0.5^0.5)
  old_v <- matrix(0L,nrow=length(p),ncol=dim(x)[2])
  count <- 0
  alt_count <- 0
  if ((m_estimator[[1]] == 'huber') | (m_estimator[[1]] == 'tukey')) { 
    xi <- (2*Pinv(dim/2,0.5))^0.5
    deviations <- vector(length=dim(y)[2])
    current_shifts <- current_v%*%t(x)
    for (i in 1:dim(y)[2]) {
      deviations[i] <- dist(expo(current_p,current_shifts[,i]),y[,i])
    }
    mad <- median(deviations)
    if (m_estimator[[1]] == 'huber') {
      c <- nr(2,(2*length(p)-4),m_estimator)
    } else if (m_estimator[[1]] == 'tukey') {
      c <- nr(15,(2*length(p)-4),m_estimator)
    }
    sigma <- mad/xi
    cutoff <- c*sigma
    rho_function <- m_estimator
    m_estimator <- vector('list',length=2)
    m_estimator[[1]] <- rho_function
    m_estimator[[2]] <- cutoff
  }
  step_p <- grad_p(current_p,current_v,x,y,m_estimator)
  step_v <- grad_v(current_p,current_v,x,y,m_estimator)
  v_diffs <- vector(length=dim(x)[2])
  for (h in 1:dim(x)[2]) {
    v_diffs[h] <- norm((pt(old_p,current_p,old_v[,h])-current_v[,h]))
  }
  lambda <- min((1/norm(step_p)),0.1)
  while ((count==0) | ((count<20000) & (alt_count<100000) & ((dist(old_p,current_p)>0.0000001) | (any(v_diffs>0.0000001))))) { ##& ((dist(copy_old_p,current_p)!=0) | (norm(current_v-copy_old_v)!=0)))) {
    new_p <- expo(current_p, -lambda*step_p)
    new_v <- matrix(,nrow=length(p),ncol=dim(x)[2])
    for (h in 1:dim(x)[2]) {
      new_v[,h] <- pt(current_p,new_p,current_v[,h]-lambda*step_v[,h])
    }
    if (loss_sum(current_p,current_v,x,y,m_estimator) >= loss_sum(new_p,new_v,x,y,m_estimator)) {
      alt_count <- 0
      old_p <- current_p
      old_v <- current_v
      current_p <- new_p
      current_v <- new_v
      if ((m_estimator[[1]] == 'huber') | (m_estimator[[1]] == 'tukey')) {
        current_shifts <- current_v%*%t(x)
        for (i in 1:dim(y)[2]) {
          deviations[i] <- dist(expo(current_p,current_shifts[,i]),y[,i])
        }
        mad <- median(deviations)
        sigma <- mad/xi
        cutoff <- c*sigma
        m_estimator[[2]] <- cutoff
      }
      step_p <- grad_p(current_p,current_v,x,y,m_estimator)
      step_v <- grad_v(current_p,current_v,x,y,m_estimator)
      for (h in 1:dim(x)[2]) {
        v_diffs[h] <- norm((pt(old_p,current_p,old_v[,h])-current_v[,h]))
      }
      lambda <- min((1/norm(step_p)),2*lambda)
      count <- count+1
    } else {
      lambda <- lambda/2
      alt_count <- alt_count+1
    }
  }
  result <- vector("list", length=5)
  result[[1]] <- current_p
  result[[2]] <- current_v
  return (result)
}
