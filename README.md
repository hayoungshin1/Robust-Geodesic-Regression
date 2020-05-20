# Robust-Geodesic-Regression

Code for the numerical experiments in the paper 'Robust Geodesic Regression' by Ha-Young Shin and Hee-Seok Oh.

finding_c.R contains the code for the approximate AREs and their derivatives for the Huber and Tukey biweight estimators, used to find the tuning parameter c for 95% efficiency. Also included is the code for the ARE of the L_1 estimator.

sphere_gradient_descent.R contains the code for performing the gradient descent algorithm to solve the geodesic regression problem on the k-sphere S^k, embedded in R(k+1). It includes the code for the exponential and logarithmic maps, parallel transport and Jacobi fields.

sphere_karcher_mean.R contains the code for finding the intrinsic mean of a dataset on the k-sphere S^k, embedded in R(k+1).

sphere_simulations.R contains the code for the simulation experiments on S^k. It requires the finding_c.R, sphere_gradient_descent.R and sphere_karcher_mean.R files to be executed first. The output is the sample variance for p and v^j, used to calculate relative efficiency, and the MSEs for p and v^j in the G, T, and C scenarios. Under initializations, the dimension k, the number of independent variables n, and the type of M-type estimator can be changed.

kendall_gradient_descent.R contains the code for performing the gradient descent algorithm to solve the geodesic regression problem on the Kendall's 2-dimensional shape space. It includes the code for the exponential and logarithmic maps, parallel transport and Jacobi fields.

knendall_karcher_mean.R contains the code for finding the intrinsic mean of a dataset on Kendall's 2-dimensional shape space.

kendall_experiments.R contains the code for the experiments on Kendall's 2-dimensional shape space. It requires the finding_c.R, kendall_gradient_descent.R and kendall_karcher_mean.R files to be executed first. The type of M-type estimator can be changed under the initializations section. We have used the preprocessed data provided by Cornea and Zhu on their website http://www.bios.unc.edu/research/bias/software.html under the title 'Regression Models on Riemannian Symmetric Spaces'.
