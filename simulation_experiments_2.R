library(mvtnorm)
library(MASS)
library(GeodRegr)

# initializations

set.seed(42)

manifold <- 'sphere'
# manifold <- 'hyperbolic'

# dim <- 2
dim <- 3
embed <- dim + 1
L <- 1024
# estimator <- 'l2'
# estimator <- 'l1'
# estimator <- 'huber'
estimator <- 'tukey'

n <- 8

true_p <- numeric(embed)
true_p[1] <- 1
sigmas <- c(pi / 32, pi / 16, pi / 8, pi / 4, pi / 2)

# data

norm <- function(v) {sqrt(sum(v * v))}

if ((estimator == 'l2') | (estimator == 'l1')) {
  c <- NULL
} else if (estimator == 'huber') {
  c <-  are_nr('huber', dim, 2)
} else if (estimator == 'tukey') {
  c <-  are_nr('tukey', dim, 5)
}

options(warn = -1)

p_variance <- vector(length = length(sigmas))

for (number in 1:length(sigmas)) {
  all_p <- matrix(, nrow = embed, ncol = L)
  for (j in 1:L) {
    train_y <- matrix(, nrow = embed, ncol = 2 ^ n)
    train_sims <- rnormtangents(manifold, N = 2 ^ n, n = dim, sigma_sq = sigmas[number] ^ 2)
    for (i in 1:2 ^ n) {
      train_y[, i] <- exp_map(manifold, true_p, train_sims[, i])
    }
    ans <- intrinsic_location(manifold, train_y, estimator, c = c, max_iter = 10000)
    all_p[, j] <- ans
  }
  p_mean <- intrinsic_location(manifold, all_p, 'l2')
  p_variance[number] <- (1 / L) * 2 * loss(manifold, p_mean, numeric(embed), t(t(numeric(L))), all_p, 'l2')
  print(p_variance[number])
}

print(p_variance)
