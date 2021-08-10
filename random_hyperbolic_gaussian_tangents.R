library(MASS)

# used in the Newton-Raphson process in the inflection_point function
xtanx <- function(x) {
  return(x * tan(x))
}

# used in the Newton-Raphson process in the inflection_point function
deriv_xtanx <- function(x) {
  return (tan(x) + x / (cos(x)^ 2))
}

# inflection point of G for dimension m given sigma_sq which will be in inv_func, found using the Newton-Raphson method
inflection_point <- function(m, sigma_sq) {
  old_x <- pi / 2 + 5
  new_x <- pi / 2 - 0.1
  count <- 0
  while (xtanx(new_x) < m * sigma_sq) {
    new_x <- pi - (pi - new_x) / 10
  }
  while (abs(old_x - new_x) > 0.000000001) {
    old_x <- new_x
    new_x <- old_x - (xtanx(old_x) - m * sigma_sq) / deriv_xtanx(old_x)
    count <- count + 1
    if (count > 1000) {
      return ('fail')
    }
  }
  return (new_x)
}

# the error function from the pracma package gives an error when z = 0; this erfz has been modified to remove the error
erfz <- function (z)
{
  if (is.null(z)) 
    return(NULL)
  else if (!is.numeric(z) && !is.complex(z)) 
    stop("Argument 'z' must be a numeric or complex scalar or vector.")
  a0 <- abs(z)
  c0 <- exp(-z * z)
  z1 <- ifelse(Re(z) < 0, -z, z)
  i <- a0 <= 5.8
  work.i <- i
  cer <- rep(NA, length = length(z))
  if (sum(work.i) > 0) {
    cs <- z1
    cr <- cs
    for (k in 1:120) {
      cr[work.i] <- cr[work.i] * z1[work.i] * z1[work.i]/(k + 
                                                            0.5)
      cs[work.i] <- cs[work.i] + cr[work.i]
    }
    cer[i] <- 2 * c0[i] * cs[i]/sqrt(pi)
  }
  work.i <- !i
  if (sum(work.i) > 0) {
    cl <- 1/z1
    cr <- cl
    for (k in 1:13) {
      cr[work.i] <- -cr[work.i] * (k - 0.5)/(z1[work.i] * 
                                               z1[work.i])
      cl[work.i] <- cl[work.i] + cr[work.i]
    }
    cer[!i] <- 1 - c0[!i] * cl[!i]/sqrt(pi)
  }
  cer[Re(z) < 0] <- -cer[Re(z) < 0]
  return(cer)
}

# distribution function of r = d(y, mu) for dimension m given sigma_sq
G <- function(m, sigma_sq, R) {
  sum <- 0
  for (j in 0:m) {
    sum <- sum + (factorial(m) / (factorial(j) * factorial(m - j))) * ((-1) ^ j) * exp((sigma_sq * (m - 2j)^2) / 2) * erfz(R / ((2 * sigma_sq) ^ 0.5) - ((sigma_sq / 2) ^ 0.5) * (m - 2 * j))
  }
  return(sum)
}

# derivative of G
deriv_G <- function(m, sigma_sq, R) {
  return (exp((-R^2) / (2 * sigma_sq)) * (sin(R) ^ m) * (2 ^ (m + 1)) * exp((sigma_sq * m ^ 2) / 2) / ((2 * pi * sigma_sq) ^ 0.5))
}

# inverse of G at t, using the Newton-Raphson method with G and deriv_G
inv_func <- function(m, sigma_sq, t) {
  new_x <- inflection_point(m, sigma_sq) # use the inflection point of G as the starting point for the Newton-Raphson method
  old_x <- new_x + 5
  count <- 0
  while (any(abs(old_x - new_x) > 0.0000000001)) {
    old_x <- new_x
    new_x <- old_x - (G(m, sigma_sq, old_x) - t) / deriv_G(m, sigma_sq, old_x)
    count <- count + 1
    if (count > 1000) {
      return ('fail')
    }
  }
  return (new_x)
}

# random generation of tangent vectors of n normally distributed points on S^k;
# more precisely, tangent vectors are of the form Log(mu, y) in the tangent space at mu, 
# which is equivalent to R^k, when y has a Riemannian Gaussian distribution.
random_sphere_gaussian_tangents <- function(n, k, sigma_sq) {
  u <- runif(n, 0, 1)
  t <- u * (G(k - 1, sigma_sq, pi) - G(k - 1, sigma_sq, 0)) + G(k - 1, sigma_sq, 0)
  magnitude <- inv_func(k - 1, sigma_sq, t)
  while (magnitude[1] == 'fail') {
    u <- runif(n, 0, 1)
    t <- u * (G(k - 1, sigma_sq, pi) - G(k - 1, sigma_sq, 0)) + G(k - 1, sigma_sq, 0)
    magnitude <- inv_func(k - 1, sigma_sq, t)
  }
  direction <- matrix(mvrnorm(n = n, mu = integer(k), Sigma = diag(x = 10, nrow = k)), nrow = n)
  direction <- direction / (rowSums(direction ^ 2) ^ 0.5)
  return (magnitude * direction)
}
