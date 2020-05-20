library(mvtnorm)
library(MASS)

## initializations

set.seed(234234)

dim <- 2
#dim <- 3
embed <- dim+1
L <- 1024
n <- 1
#n <- 2
max <- 8
size <- 7
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

graph_p <- integer(max-1)
graph_v <- matrix(0L,nrow=(max-1),ncol=n)

t_graph_p <- integer(max-1)
t_graph_v <- matrix(0L,nrow=(max-1),ncol=n)

mixed_graph_p <- integer(max-1)
mixed_graph_v <- matrix(0L,nrow=(max-1),ncol=n)

p_variance <- integer(max-1)
v_variance <- matrix(0L,nrow=(max-1),ncol=n)

for (k in 2:max) {
  mse_p <- 0
  mse_v <- integer(n)
  t_mse_p <- 0
  t_mse_v <- integer(n)
  mixed_mse_p <- 0
  mixed_mse_v <- integer(n)
  all_p <- matrix(,nrow=embed,ncol=L)
  all_v <- array(0L,c(embed,n,L))
  for (j in 1:L) {
    train_x <- matrix(runif(n*(2^max)),ncol=n)-0.5
    train_y <- matrix(,nrow=embed,ncol=2^max)
    t_train_y <- matrix(,nrow=embed,ncol=2^max)
    mixed_train_y <- matrix(,nrow=embed,ncol=2^max)
    train_sims <- matrix(mvrnorm(n=2^max, mu=integer(dim), Sigma=diag(x=sigma^2, nrow=dim)),nrow=2^max)
    t_train_sims <- matrix(rmvt(2^max,sigma=diag(x=t_sigma^2,nrow=dim),df=3,delta=integer(dim),type="shifted"),nrow=2^max)
    bern <- rbinom(2^max,1,0.9)
    mixed_train_sims <- (1-bern)*matrix(mvrnorm(n=2^max, mu=integer(dim), Sigma=diag(x=mixed_sigma_1^2, nrow=dim)),nrow=2^max)+bern*matrix(mvrnorm(n=2^max, mu=integer(dim), Sigma=diag(x=mixed_sigma_2^2, nrow=dim)),nrow=2^max)
    true_shifts <- true_v%*%t(train_x)
    for (i in 1:2^max) {
      train_error <- pt(true_p,expo(true_p,true_shifts[,i]),c(0,train_sims[i,]))
      train_y[,i] <- expo(expo(true_p,true_shifts[,i]),train_error)
      t_train_error <- pt(true_p,expo(true_p,true_shifts[,i]),c(0,t_train_sims[i,]))
      t_train_y[,i] <- expo(expo(true_p,true_shifts[,i]),t_train_error)
      mixed_train_error <- pt(true_p,expo(true_p,true_shifts[,i]),c(0,mixed_train_sims[i,]))
      mixed_train_y[,i] <- expo(expo(true_p,true_shifts[,i]),mixed_train_error)
    }
    train_x <- t(t(train_x[1:2^k,]))
    train_y <- train_y[,1:2^k]
    t_train_y <- t_train_y[,1:2^k]
    mixed_train_y <- mixed_train_y[,1:2^k]
    init_p <- k_mean(train_y)
    init_v <- matrix(0L,nrow=embed,ncol=n)
    t_init_p <- k_mean(t_train_y)
    mixed_init_p <- k_mean(mixed_train_y)
    ans <- alg(init_p,init_v,train_x,train_y,estimator)
    t_ans <- alg(t_init_p,init_v,train_x,t_train_y,estimator)
    mixed_ans <- alg(mixed_init_p,init_v,train_x,mixed_train_y,estimator)
    all_p[,j] <- ans[[1]]
    all_v[,,j] <- ans[[2]]
    mse_p <- mse_p + (dist(ans[[1]],true_p))^2
    for (h in 1:length(train_x[1,])) {
      mse_v[h] <- mse_v[h] + (norm(pt(ans[[1]],true_p,ans[[2]][,h])-true_v[,h]))^2
    }
    t_mse_p <- t_mse_p + (dist(t_ans[[1]],true_p))^2
    for (h in 1:length(train_x[1,])) {
      t_mse_v[h] <- t_mse_v[h] + (norm(pt(t_ans[[1]],true_p,t_ans[[2]][,h])-true_v[,h]))^2
    }
    mixed_mse_p <- mixed_mse_p + (dist(mixed_ans[[1]],true_p))^2
    for (h in 1:length(train_x[1,])) {
      mixed_mse_v[h] <- mixed_mse_v[h] + (norm(pt(mixed_ans[[1]],true_p,mixed_ans[[2]][,h])-true_v[,h]))^2
    }
  }
  p_mean <- k_mean(all_p)
  p_variance[k-1] <- (1/L)*2*k_mean_loss_sum(p_mean,all_p)
  for (j in 1:L) {
    for (h in 1:n) {
      all_v[,h,j] <- pt(all_p[,j],p_mean,all_v[,h,j])
    }
  }
  v_mean <- rowSums(all_v,dims=2)/L
  for (j in 1:L) {
    for (h in 1:n) {
      v_variance[(k-1),h] <- v_variance[(k-1),h]+(norm(all_v[,h,j]-v_mean[,h]))^2
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
}

v_variance <- v_variance/L


p_variance
t(v_variance)

graph_p
t(graph_v)

t_graph_p
t(t_graph_v)

mixed_graph_p
t(mixed_graph_v)
