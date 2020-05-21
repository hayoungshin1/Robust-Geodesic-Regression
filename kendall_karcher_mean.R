## Karcher mean on Kendall's 2d shape space

k_mean_loss_sum <- function(p, y) {
  sum <- 0
  for (i in 1:dim(y)[2]) {
    sum <- sum + 0.5*((dist(p, y[, i]))^2)
  }
  return (sum)
}

k_mean_grad <- function(p, y) {
  sum <- integer(length(p))
  for (i in 1:dim(y)[2]) {
    sum <- sum - loga(p, y[, i])
  }
  return (sum)
}

k_mean <- function(y) {
  current_p <- y[, 1]
  old_p <- y[, 2]
  lambda <- 0.1
  step_p <- k_mean_grad(current_p, y)
  count <- 0
  while ((count==0) | (dist(old_p, current_p) > 0.0000001)) {
    lambda <- min((1/norm(step_p)),lambda)
    new_p <- expo(current_p, -lambda*step_p)
    if (k_mean_loss_sum(current_p, y) >= k_mean_loss_sum(new_p, y)) {
      old_p <- current_p
      current_p <- new_p
      step_p <- k_mean_grad(current_p, y)
      lambda <- 2*lambda
      count <- count+1
    } else {
      lambda <- lambda/2
    }
  }
  result <- current_p
  return (result) 
}
