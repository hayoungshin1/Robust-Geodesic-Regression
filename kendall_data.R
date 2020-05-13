library(zipfR)
library(MASS)
library(R.matlab)
library(rmatio)
library(abind)
library(data.table)

#data

boundary_points <- 50
dim <- 2*boundary_points-4
embed <- boundary_points
#M <- 1024
l <- 1
estimator <- 'l1'

ADF_x_data <- readMat('../downloads/ADNI_RMRSS/data/ADinfoF.mat')$ADinfoF

ADF_y_data <- readMat('../downloads/ADNI_RMRSS/data/ADLdataF.mat')$ADLdataF

ages <- ADF_x_data[,9]

x_data <- t(t(ages))
x_data <- x_data-mean(x_data) ## centering
#x_data <- c(x_data,0)

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

## y_data[2,,(length(ages)-19):length(ages)] <- -y_data[2,,(length(ages)-19):length(ages)] ### fake points

y_data <- y_data[1,,]+y_data[2,,]*(0+1i)

## solutions

init_p <- k_mean(y_data)
init_v <- t(t(integer(embed)))

ans <- alg(init_p,init_v,x_data,y_data,estimator)

total_var <- k_mean_loss_sum(k_mean(y_data),y_data)
unexp_var <- loss_sum(ans[[1]],ans[[2]],x_data,y_data,'l2')
r_sq <- 1-unexp_var/total_var

r_sq

ans[[1]]
t(ans[[2]])
sum(ans[[1]]*Conj(ans[[2]]))
norm(ans[[1]])
sum(ans[[1]])
sum(ans[[2]])
ans[[3]]
ans[[4]]
ans[[5]]
