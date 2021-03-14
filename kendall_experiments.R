library(MASS)
library(R.matlab)
library(rmatio)

# we will test the robustness of each estimator by comparing their
# performance on the original (corrupted) data set to that of the L_2
# estimator on the uncorrupted data set (with the 4 problematic data points
# removed).
data(calvaria)
manifold <- 'kendall'
contam_x_data <- calvaria$x
contam_mean_x <- mean(contam_x_data)
contam_x_data <- contam_x_data - contam_mean_x # center x data
uncontam_x_data <- calvaria$x[ -c(23, 101, 104, 160)]
uncontam_mean_x <- mean(uncontam_x_data)
uncontam_x_data <- uncontam_x_data - uncontam_mean_x # center x data
contam_y_data <- calvaria$y
uncontam_y_data <- calvaria$y[, -c(23, 101, 104, 160)] # remove corrupted
# columns
landmarks <- dim(contam_y_data)[1]
dimension <- 2 * landmarks - 4
# we ignore Huber's estimator as the L_1 estimator already has an
# (approximate) efficiency above 95% in 12 dimensions; see documentation for
# the are and are_nr functions
tol <- 1e-5
uncontam_l2 <- geo_reg(manifold, uncontam_x_data, uncontam_y_data,
'l2', p_tol = tol, V_tol = tol)
contam_l2 <- geo_reg(manifold, contam_x_data, contam_y_data,
'l2', p_tol = tol, V_tol = tol)
contam_l1 <- geo_reg(manifold, contam_x_data, contam_y_data,
'l1', p_tol = tol, V_tol = tol)
contam_tukey <- geo_reg(manifold, contam_x_data, contam_y_data,
'tukey', are_nr('tukey', dimension, 10, 0.99), p_tol = tol, V_tol = tol)
geodesics <- vector('list')
geodesics[[1]] <- uncontam_l2
geodesics[[2]] <- contam_l2
geodesics[[3]] <- contam_l1
geodesics[[4]] <- contam_tukey
loss(manifold, geodesics[[1]]$p, geodesics[[1]]$V, uncontam_x_data,
uncontam_y_data, 'l2')
loss(manifold, geodesics[[2]]$p, geodesics[[2]]$V, contam_x_data,
contam_y_data, 'l2')
loss(manifold, geodesics[[3]]$p, geodesics[[3]]$V, contam_x_data,
contam_y_data, 'l1')
loss(manifold, geodesics[[4]]$p, geodesics[[4]]$V, contam_x_data,
contam_y_data, 'tukey', are_nr('tukey', dimension, 10, 0.99))
6 exp_map
# visualization of each geodesic
oldpar <- par(mfrow = c(1, 4))
days <- c(7, 14, 21, 30, 40, 60, 90, 150)
pal <- colorRampPalette(c("blue", "red"))(length(days))
# each predicted geodesic will be represented as a sequence of the predicted
# shapes at each of the above ages, the blue contour will show the predicted
# shape on day 7 and the red contour the predicted shape on day 150
contour <- vector('list')
for (i in 1:length(days)) {
contour[[i]] <- exp_map(manifold, geodesics[[1]]$p, (days[i] -
uncontam_mean_x) * geodesics[[1]]$V)
contour[[i]] <- c(contour[[i]], contour[[i]][1])
}
plot(Re(contour[[length(days)]]), Im(contour[[length(days)]]), type = 'n',
xaxt = 'n', yaxt = 'n', ann = FALSE, asp = 1)
for (i in 1:length(days)) {
lines(Re(contour[[i]]), Im(contour[[i]]), col = pal[i])
}
for (j in 2:4) {
for (i in 1:length(days)) {
contour[[i]] <- exp_map(manifold, geodesics[[j]]$p, (days[i] -
contam_mean_x) * geodesics[[j]]$V)
contour[[i]] <- c(contour[[i]], contour[[i]][1])
}
plot(Re(contour[[length(days)]]), Im(contour[[length(days)]]), type = 'n',
xaxt = 'n', yaxt = 'n', ann = FALSE, asp = 1)
for (i in 1:length(days)) {
lines(Re(contour[[i]]), Im(contour[[i]]), col = pal[i])
}
}
# even with a mere 4 corrupted landmarks out of a total of 8 * 168 = 1344, we
# can clearly see that contam_l2, the second image, looks slightly
# different from all the others, especially near the top of the image.
par(oldpar)
