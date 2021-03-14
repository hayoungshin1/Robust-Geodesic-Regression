rm(list = ls())

library(mvtnorm)
library(MASS)
library(GeodRegr)

set.seed(23423)

manifold <- 'sphere'

dim <- 2
#dim <- 3
embed <- dim+1
L <- 1024
n <- 1
#n <- 2
max <- 8
estimator <- 'l2'
#estimator <- 'l1'
#estimator <- 'huber'
#estimator <- 'tukey'

k <- max

true_p <- integer(embed)
true_p[1] <- 1
true_v <- matrix(0L,nrow=embed,ncol=n)
true_v[2,1] <- pi/4
true_v[embed,2] <- -pi/6
sigma_1 <- pi/8
sigma_2 <- pi/6
sigma_3 <- pi/4
sigma_4 <- pi/2

v_variance_1 <- vector(length=n)
v_variance_2 <- vector(length=n)
v_variance_3 <- vector(length=n)
v_variance_4 <- vector(length=n)

# data

norm <- function(v) {sqrt(sum(v*v))}

if ((estimator == 'l2') | (estimator == 'l1')) {
  c <- NULL
} else if (estimator == 'huber') {
  c <-  are_nr('huber', dim, 2)
} else if (estimator == 'tukey') {
  c <-  are_nr('tukey', dim, 5)
}

options(warn=-1)

all_p_1 <- matrix(,nrow=embed,ncol=L)
all_v_1 <- array(0L,c(embed,n,L))
all_p_2 <- matrix(,nrow=embed,ncol=L)
all_v_2 <- array(0L,c(embed,n,L))
all_p_3 <- matrix(,nrow=embed,ncol=L)
all_v_3 <- array(0L,c(embed,n,L))
all_p_4 <- matrix(,nrow=embed,ncol=L)
all_v_4 <- array(0L,c(embed,n,L))
for (j in 1:L) {
  train_x <- matrix(runif(n*(2^max)),ncol=n)-0.5
  train_y_1 <- matrix(,nrow=embed,ncol=2^max)
  train_y_2 <- matrix(,nrow=embed,ncol=2^max)
  train_y_3 <- matrix(,nrow=embed,ncol=2^max)
  train_y_4 <- matrix(,nrow=embed,ncol=2^max)
  train_sims_1 <- random_sphere_gaussian_tangents(n=2^max, k=dim, sigma_sq=sigma_1^2)
  train_sims_2 <- random_sphere_gaussian_tangents(n=2^max, k=dim, sigma_sq=sigma_2^2)
  train_sims_3 <- random_sphere_gaussian_tangents(n=2^max, k=dim, sigma_sq=sigma_3^2)
  train_sims_4 <- random_sphere_gaussian_tangents(n=2^max, k=dim, sigma_sq=sigma_4^2)
  true_shifts <- true_v%*%t(train_x)
  for (i in 1:2^max) {
    train_error_1 <- par_trans(manifold, true_p,exp_map(manifold, true_p,true_shifts[,i]),c(0,train_sims_1[i,]))
    train_y_1[,i] <- exp_map(manifold, exp_map(manifold, true_p,true_shifts[,i]),train_error_1)
    train_error_2 <- par_trans(manifold, true_p,exp_map(manifold, true_p,true_shifts[,i]),c(0,train_sims_2[i,]))
    train_y_2[,i] <- exp_map(manifold, exp_map(manifold, true_p,true_shifts[,i]),train_error_2)
    train_error_3 <- par_trans(manifold, true_p,exp_map(manifold, true_p,true_shifts[,i]),c(0,train_sims_3[i,]))
    train_y_3[,i] <- exp_map(manifold, exp_map(manifold, true_p,true_shifts[,i]),train_error_3)
    train_error_4 <- par_trans(manifold, true_p,exp_map(manifold, true_p,true_shifts[,i]),c(0,train_sims_4[i,]))
    train_y_4[,i] <- exp_map(manifold, exp_map(manifold, true_p,true_shifts[,i]),train_error_4)
  }
  train_x <- t(t(train_x[1:2^k,]))
  train_y_1 <- train_y_1[,1:2^k]
  train_y_2 <- train_y_2[,1:2^k]
  train_y_3 <- train_y_3[,1:2^k]
  train_y_4 <- train_y_4[,1:2^k]
  ans_1 <- geo_reg(manifold, train_x,train_y_1,estimator, c = c, max_iter = 1000)
  ans_2 <- geo_reg(manifold, train_x,train_y_2,estimator, c = c, max_iter = 1000)
  ans_3 <- geo_reg(manifold, train_x,train_y_3,estimator, c = c, max_iter = 1000)
  ans_4 <- geo_reg(manifold, train_x,train_y_4,estimator, c = c, max_iter = 1000)
  all_p_1[,j] <- ans_1[[1]]
  all_v_1[,,j] <- ans_1[[2]]
  all_p_2[,j] <- ans_2[[1]]
  all_v_2[,,j] <- ans_2[[2]]
  all_p_3[,j] <- ans_3[[1]]
  all_v_3[,,j] <- ans_3[[2]]
  all_p_4[,j] <- ans_4[[1]]
  all_v_4[,,j] <- ans_4[[2]]
}
p_mean_1 <- intrinsic_location(manifold, all_p_1, 'l2')
p_variance_1 <- (1/L)*2*loss(manifold, p_mean_1, integer(embed), t(t(integer(L))),all_p_1, 'l2')
for (j in 1:L) {
  for (h in 1:n) {
    all_v_1[,h,j] <- par_trans(manifold, all_p_1[,j],p_mean_1,all_v_1[,h,j])
  }
}
v_mean_1 <- rowSums(all_v_1,dims=2)/L
for (j in 1:L) {
  for (h in 1:n) {
    v_variance_1[h] <- v_variance_1[h]+(norm(all_v_1[,h,j]-v_mean_1[,h]))^2
  }
}
p_mean_2 <- intrinsic_location(manifold, all_p_2, 'l2')
p_variance_2 <- (1/L)*2*loss(manifold, p_mean_2, integer(embed), t(t(integer(L))),all_p_2, 'l2')
for (j in 1:L) {
  for (h in 1:n) {
    all_v_2[,h,j] <- par_trans(manifold, all_p_2[,j],p_mean_2,all_v_2[,h,j])
  }
}
v_mean_2 <- rowSums(all_v_2,dims=2)/L
for (j in 1:L) {
  for (h in 1:n) {
    v_variance_2[h] <- v_variance_2[h]+(norm(all_v_2[,h,j]-v_mean_2[,h]))^2
  }
}
p_mean_3 <- intrinsic_location(manifold, all_p_3, 'l2')
p_variance_3 <- (1/L)*2*loss(manifold, p_mean_3, integer(embed), t(t(integer(L))),all_p_3, 'l2')
for (j in 1:L) {
  for (h in 1:n) {
    all_v_3[,h,j] <- par_trans(manifold, all_p_3[,j],p_mean_3,all_v_3[,h,j])
  }
}
v_mean_3 <- rowSums(all_v_3,dims=2)/L
for (j in 1:L) {
  for (h in 1:n) {
    v_variance_3[h] <- v_variance_3[h]+(norm(all_v_3[,h,j]-v_mean_3[,h]))^2
  }
}
p_mean_4 <- intrinsic_location(manifold, all_p_4, 'l2')
p_variance_4 <- (1/L)*2*loss(manifold, p_mean_4, integer(embed), t(t(integer(L))),all_p_4, 'l2')
for (j in 1:L) {
  for (h in 1:n) {
    all_v_4[,h,j] <- par_trans(manifold, all_p_4[,j],p_mean_4,all_v_4[,h,j])
  }
}
v_mean_4 <- rowSums(all_v_4,dims=2)/L
for (j in 1:L) {
  for (h in 1:n) {
    v_variance_4[h] <- v_variance_4[h]+(norm(all_v_4[,h,j]-v_mean_4[,h]))^2
  }
}

v_variance_1 <- v_variance_1/(L-1)
v_variance_2 <- v_variance_2/(L-1)
v_variance_3 <- v_variance_3/(L-1)
v_variance_4 <- v_variance_4/(L-1)

print(p_variance_1)
print(v_variance_1)

print(p_variance_2)
print(v_variance_2)

print(p_variance_3)
print(v_variance_3)

print(p_variance_4)
print(v_variance_4)
