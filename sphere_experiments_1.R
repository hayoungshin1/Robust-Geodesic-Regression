rm(list = ls())

library(mvtnorm)
library(MASS)
library(GeodRegr)

# initializations

set.seed(23423)

manifold <- 'sphere'

#dim <- 2
dim <- 3
embed <- dim+1
L <- 1024
#n <- 1
n <- 2
max <- 8
estimator <- 'l2'
#estimator <- 'l1'
#estimator <- 'huber'
#estimator <- 'tukey'

true_p <- integer(embed)
true_p[1] <- 1
true_v <- matrix(0L,nrow=embed,ncol=n)
true_v[2,1] <- pi/4
true_v[embed,2] <- -pi/6
sigma <- pi/8
t_sigma <- pi/16
mixed_sigma_1 <- pi/6
mixed_sigma_2 <- pi/24

# data

graph_p <- integer(max-1) # used to store MSE(p)s
graph_v <- matrix(0L,nrow=(max-1),ncol=n) # used to store MSE(v^j)s

t_graph_p <- integer(max-1)
t_graph_v <- matrix(0L,nrow=(max-1),ncol=n)

mixed_graph_p <- integer(max-1)
mixed_graph_v <- matrix(0L,nrow=(max-1),ncol=n)


norm <- function(v) {sqrt(sum(v*v))}

if ((estimator == 'l2') | (estimator == 'l1')) {
  c <- NULL
} else if (estimator == 'huber') {
  c <-  are_nr('huber', dim, 2)
} else if (estimator == 'tukey') {
  c <-  are_nr('tukey', dim, 5)
}

options(warn=-1)

for (k in 2:max) {
  mse_p <- 0
  mse_v <- integer(n)
  t_mse_p <- 0
  t_mse_v <- integer(n)
  mixed_mse_p <- 0
  mixed_mse_v <- integer(n)
  for (j in 1:L) {
    train_x <- matrix(runif(n*(2^max)),ncol=n)-0.5
    train_y <- matrix(,nrow=embed,ncol=2^max)
    t_train_y <- matrix(,nrow=embed,ncol=2^max)
    mixed_train_y <- matrix(,nrow=embed,ncol=2^max)
    train_sims <- random_sphere_gaussian_tangents(n=2^max, k=dim, sigma_sq=sigma^2)
    t_train_sims <- matrix(rmvt(2^max,sigma=diag(x=t_sigma^2,nrow=dim),df=3,delta=integer(dim),type="shifted"),nrow=2^max)
    bern <- rbinom(2^max,1,0.9)
    mixed_train_sims <- (1-bern)*random_sphere_gaussian_tangents(n=2^max, k=dim, sigma_sq=mixed_sigma_1^2)+bern*random_sphere_gaussian_tangents(n=2^max, k=dim, sigma_sq=mixed_sigma_2^2)
    true_shifts <- true_v%*%t(train_x)
    for (i in 1:2^max) {
      train_error <- par_trans(manifold, true_p,exp_map(manifold, true_p,true_shifts[,i]),c(0,train_sims[i,]))
      train_y[,i] <- exp_map(manifold, exp_map(manifold, true_p,true_shifts[,i]),train_error)
      t_train_error <- par_trans(manifold, true_p,exp_map(manifold, true_p,true_shifts[,i]),c(0,t_train_sims[i,]))
      t_train_y[,i] <- exp_map(manifold, exp_map(manifold, true_p,true_shifts[,i]),t_train_error)
      mixed_train_error <- par_trans(manifold, true_p,exp_map(manifold, true_p,true_shifts[,i]),c(0,mixed_train_sims[i,]))
      mixed_train_y[,i] <- exp_map(manifold, exp_map(manifold, true_p,true_shifts[,i]),mixed_train_error)
    }
    train_x <- t(t(train_x[1:2^k,]))
    train_y <- train_y[,1:2^k]
    t_train_y <- t_train_y[,1:2^k]
    mixed_train_y <- mixed_train_y[,1:2^k]
    ans <- geo_reg(manifold, train_x,train_y,estimator, c = c, max_iter = 1000)
    t_ans <- geo_reg(manifold, train_x,t_train_y,estimator, c = c, max_iter = 1000)
    mixed_ans <- geo_reg(manifold, train_x,mixed_train_y,estimator, c = c, max_iter = 1000)
    mse_p <- mse_p + (geo_dist(manifold, ans[[1]],true_p))^2
    for (h in 1:length(train_x[1,])) {
      mse_v[h] <- mse_v[h] + (norm(par_trans(manifold, ans[[1]],true_p,ans[[2]][,h])-true_v[,h]))^2
    }
    t_mse_p <- t_mse_p + (geo_dist(manifold, t_ans[[1]],true_p))^2
    for (h in 1:length(train_x[1,])) {
      t_mse_v[h] <- t_mse_v[h] + (norm(par_trans(manifold, t_ans[[1]],true_p,t_ans[[2]][,h])-true_v[,h]))^2
    }
    mixed_mse_p <- mixed_mse_p + (geo_dist(manifold, mixed_ans[[1]],true_p))^2
    for (h in 1:length(train_x[1,])) {
      mixed_mse_v[h] <- mixed_mse_v[h] + (norm(par_trans(manifold, mixed_ans[[1]],true_p,mixed_ans[[2]][,h])-true_v[,h]))^2
    }
  }
  mse_p <- mse_p/L
  mse_v <- mse_v/L
  t_mse_p <- t_mse_p/L
  t_mse_v <- t_mse_v/L
  mixed_mse_p <- mixed_mse_p/L
  mixed_mse_v <- mixed_mse_v/L
  graph_p[k-1] <- mse_p
  graph_v[k-1,] <- mse_v
  t_graph_p[k-1] <- t_mse_p
  t_graph_v[k-1,] <- t_mse_v
  mixed_graph_p[k-1] <- mixed_mse_p
  mixed_graph_v[k-1,] <- mixed_mse_v
  print(k)
  print(mse_p)
  print(mse_v)
  print(t_mse_p)
  print(t_mse_v)
  print(mixed_mse_p)
  print(mixed_mse_v)
}

print(graph_p)
print(t(graph_v))

print(t_graph_p)
print(t(t_graph_v))

print(mixed_graph_p)
print(t(mixed_graph_v))
