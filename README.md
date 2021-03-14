# Robust-Geodesic-Regression

Code for the numerical experiments in the paper 'Robust Geodesic Regression' by Ha-Young Shin and Hee-Seok Oh.

sphere_gaussian_tangents.R

sphere_experiments_1.R contains the code for the simulation experiments on S^k. It requires the sphere_gaussian_tangents.R file to be executed first. The output is the sample variance for p and v^j, used to calculate relative efficiency, and the MSEs for p and v^j in the G, T, and C scenarios. Under initializations, the dimension, the number of independent variables n and the type of M-type estimator can be changed.

sphere_experiments_2.R contains the code for the simulation experiments on S^k. It requires the sphere_gaussian_tangents.R file to be executed first. The output is the sample variance for p and v^j, used to calculate relative efficiency, and the MSEs for p and v^j in the G, T, and C scenarios. Under initializations, the dimension, the number of independent variables n and the type of M-type estimator can be changed.

kendall_experiments.R contains the code for the experiments on Kendall's 2-dimensional shape space. The type of M-type estimator can be changed under the initializations section. We have used the preprocessed data provided by Cornea and Zhu on their website http://www.bios.unc.edu/research/bias/software.html under the title 'Regression Models on Riemannian Symmetric Spaces'.
