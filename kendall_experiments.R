library(MASS)
library(R.matlab)
library(rmatio)

#initializations

manifold <- 'kendall'

boundary_points <- 50
dim <- 2*boundary_points-4
embed <- boundary_points
n <- 1
estimator <- 'l2'
#estimator <- 'l1'
#estimator <- 'tukey'

# data

ADF_x_data <- readMat('../downloads/ADNI_RMRSS/data/ADinfoF.mat')$ADinfoF
ADF_y_data <- readMat('../downloads/ADNI_RMRSS/data/ADLdataF.mat')$ADLdataF

ages <- ADF_x_data[,9]

x_data <- t(t(ages))
x_data <- x_data-mean(x_data) ## center x

y_data <- ADF_y_data
y_data <- aperm(y_data, c(2,1,3))

for (i in 1:length(ages)) { ## remove translation
  for (j in 1:2) {
    y_data[j,,i] <- y_data[j,,i]-mean(y_data[j,,i])
  }
}

for (i in 1:length(ages)) { ## remove scaling
  y_data[,,i] <- y_data[,,i]/((sum(y_data[,,i]*y_data[,,i]))^0.5)
}

#y_data[2,,(length(ages)-19):length(ages)] <- -y_data[2,,(length(ages)-19):length(ages)] ### tampered data

y_data <- y_data[1,,]+y_data[2,,]*(0+1i)

ans <- geo_reg(manifold,x_data,y_data,estimator)
