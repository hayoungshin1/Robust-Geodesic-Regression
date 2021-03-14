# Robust-Geodesic-Regression

Code for the numerical experiments in Section 4.1 of the paper 'Robust Geodesic Regression' by Ha-Young Shin and Hee-Seok Oh. The 'GeodRegr' package is required, as are 'MASS', 'mvtnorm', 'R.matlab' and 'rmatio'.

sphere_gaussian_tangents.R contains the code for generating random vectors that are the Log of points on the sphere S^k distributed according to the Riemannian Gaussian distribution for given sigma squared. So tangent vectors are of the form Log(mu, y) in the tangent space at mu, which is equivalent to R^k, when y has a Riemannian Gaussian distribution.

sphere_experiments_1.R contains the code for the simulation experiments on S^k for type G, T and C errors. It requires the sphere_gaussian_tangents.R file to be executed first. The output is the MSEs for p and v^j. Under initializations, the dimension, the number of independent variables n and the type of M-type estimator can be changed.

sphere_experiments_2.R contains the code for the simulation experiments on S^k for various values of sigma. It requires the sphere_gaussian_tangents.R file to be executed first. The output is the sample variance for p and v^j. Under initializations, the dimension, the number of independent variables n and the type of M-type estimator can be changed.

kendall_experiments.R contains the code for the experiments on Kendall's 2-dimensional shape space. The type of M-type estimator can be changed under the initializations section. We have used the preprocessed data provided by Cornea and Zhu on their website http://www.bios.unc.edu/research/bias/software.html under the title 'Regression Models on Riemannian Symmetric Spaces'.
