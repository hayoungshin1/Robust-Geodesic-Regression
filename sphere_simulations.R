library(mvtnorm)
library(zipfR)
library(MASS)

## initializations

dim <- 2
embed <- dim+1
L <- 1024
l <- 1
max <- 8
size <- 7
estimator <- 'huber'

true_p <- integer(embed)
true_p[1] <- 1
true_v <- matrix(0L,nrow=embed,ncol=l)
true_v[2,1] <- pi/4
true_v[embed,2] <- -pi/6
sigma <- pi/8
t_sigma <- pi/16
mixed_sigma_1 <- pi/6
mixed_sigma_2 <- pi/24

# data

set.seed(234234234)

graph_p <- integer(max-1)
graph_v <- matrix(0L,nrow=(max-1),ncol=l)

t_graph_p <- integer(max-1)
t_graph_v <- matrix(0L,nrow=(max-1),ncol=l)

mixed_graph_p <- integer(max-1)
mixed_graph_v <- matrix(0L,nrow=(max-1),ncol=l)

p_variance <- integer(max-1)
v_variance <- matrix(0L,nrow=(max-1),ncol=l)

t_p_variance <- integer(max-1)
t_v_variance <- matrix(0L,nrow=(max-1),ncol=l)

mixed_p_variance <- integer(max-1)
mixed_v_variance <- matrix(0L,nrow=(max-1),ncol=l)


for (k in 2:max) {
  mse_p <- 0
  mse_v <- integer(l)
  t_mse_p <- 0
  t_mse_v <- integer(l)
  mixed_mse_p <- 0
  mixed_mse_v <- integer(l)
  if ((estimator == 'huber') | (estimator == 'tukey')) {
    sigma_sq <- 0
    t_sigma_sq <- 0
    mixed_sigma_sq <- 0
  }
  weird <- 0
  strange <- 0
  t_weird <- 0
  t_strange <- 0
  mixed_weird <- 0
  mixed_strange <- 0
  all_p <- matrix(,nrow=embed,ncol=L)
  t_all_p <- matrix(,nrow=embed,ncol=L)
  mixed_all_p <- matrix(,nrow=embed,ncol=L)
  all_v <- array(0L,c(embed,l,L))
  t_all_v <- array(0L,c(embed,l,L))
  mixed_all_v <- array(0L,c(embed,l,L))
  for (j in 1:L) {
    train_x <- matrix(runif(l*(2^max)),ncol=l)-0.5
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
    init_v <- matrix(0L,nrow=embed,ncol=l)
    t_init_p <- k_mean(t_train_y)
    t_init_v <- matrix(0L,nrow=embed,ncol=l)
    mixed_init_p <- k_mean(mixed_train_y)
    mixed_init_v <- matrix(0L,nrow=embed,ncol=l)
    ans <- alg(init_p,init_v,train_x,train_y,estimator)
    t_ans <- alg(t_init_p,t_init_v,train_x,t_train_y,estimator)
    mixed_ans <- alg(mixed_init_p,mixed_init_v,train_x,mixed_train_y,estimator)
    all_p[,j] <- ans[[1]]
    t_all_p[,j] <- t_ans[[1]]
    mixed_all_p[,j] <- mixed_ans[[1]]
    all_v[,,j] <- ans[[2]]
    t_all_v[,,j] <- t_ans[[2]]
    mixed_all_v[,,j] <- mixed_ans[[2]]
    weird <- weird + (ans[[3]]==2000)
    strange <- strange + (ans[[4]]==100000)
    t_weird <- t_weird + (t_ans[[3]]==2000)
    t_strange <- t_strange + (t_ans[[4]]==100000)
    mixed_weird <- mixed_weird + (mixed_ans[[3]]==2000)
    mixed_strange <- mixed_strange + (mixed_ans[[4]]==100000)
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
    if ((estimator == 'huber') | (estimator == 'tukey')) {
      sigma_sq <- sigma_sq + ans[[5]]^2
      t_sigma_sq <- t_sigma_sq + t_ans[[5]]^2
      mixed_sigma_sq <- mixed_sigma_sq + mixed_ans[[5]]^2
    }
  }
  p_mean <- k_mean(all_p)
  p_variance[k-1] <- (1/L)*2*k_mean_loss_sum(p_mean,all_p)
  t_p_mean <- k_mean(t_all_p)
  t_p_variance[k-1] <- (1/L)*2*k_mean_loss_sum(t_p_mean,t_all_p)
  mixed_p_mean <- k_mean(mixed_all_p)
  mixed_p_variance[k-1] <- (1/L)*2*k_mean_loss_sum(mixed_p_mean,mixed_all_p)
  for (j in 1:L) {
    for (h in 1:l) {
      all_v[,h,j] <- pt(all_p[,j],p_mean,all_v[,h,j])
      t_all_v[,h,j] <- pt(t_all_p[,j],t_p_mean,t_all_v[,h,j])
      mixed_all_v[,h,j] <- pt(mixed_all_p[,j],mixed_p_mean,mixed_all_v[,h,j])
    }
  }
  v_mean <- rowSums(all_v,dims=2)/L
  t_v_mean <- rowSums(t_all_v,dims=2)/L
  mixed_v_mean <- rowSums(mixed_all_v,dims=2)/L
  for (j in 1:L) {
    for (h in 1:l) {
      v_variance[(k-1),h] <- v_variance[(k-1),h]+(norm(all_v[,h,j]-v_mean[,h]))^2
      t_v_variance[(k-1),h] <- t_v_variance[(k-1),h]+(norm(t_all_v[,h,j]-t_v_mean[,h]))^2
      mixed_v_variance[(k-1),h] <- mixed_v_variance[(k-1),h]+(norm(mixed_all_v[,h,j]-mixed_v_mean[,h]))^2
    }
  }
  v_variance <- v_variance/L
  t_v_variance <- t_v_variance/L
  mixed_v_variance <- mixed_v_variance/L
  mse_p <- mse_p/L
  mse_v <- mse_v/L
  t_mse_p <- t_mse_p/L
  t_mse_v <- t_mse_v/L
  mixed_mse_p <- mixed_mse_p/L
  mixed_mse_v <- mixed_mse_v/L
  if ((estimator == 'huber') | (estimator == 'tukey')) {
    sigma_sq <- sigma_sq/L
    t_sigma_sq <- t_sigma_sq/L
    mixed_sigma_sq <- mixed_sigma_sq/L
  }
  graph_p[k-1] <- mse_p
  graph_v[k-1,] <- mse_v
  t_graph_p[k-1] <- t_mse_p
  t_graph_v[k-1,] <- t_mse_v
  mixed_graph_p[k-1] <- mixed_mse_p
  mixed_graph_v[k-1,] <- mixed_mse_v
  print(k)
  print(weird)
  print(strange)
  print(mse_p)
  print(mse_v)
  print(t_weird)
  print(t_strange)
  print(t_mse_p)
  print(t_mse_v)
  print(mixed_weird)
  print(mixed_strange)
  print(mixed_mse_p)
  print(mixed_mse_v)
}

graph_p
t(graph_v)

t_graph_p
t(t_graph_v)

mixed_graph_p
t(mixed_graph_v)

p_variance
t(v_variance)

t_p_variance
t(t_v_variance)

mixed_p_variance
t(mixed_v_variance)
