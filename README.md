# Robust-Geodesic-Regression

Code for the numerical experiments in Section 4.1 of the paper 'Robust Geodesic Regression' by Ha-Young Shin and Hee-Seok Oh. The '[GeodRegr](https://github.com/hayoungshin1/GeodRegr)' package is required, as are 'MASS', 'mvtnorm', 'R.matlab' and 'rmatio'.

simulation_experiments_1.R contains the code for the simulation experiments on the sphere S^n and hyperbolic space \mathbb{H}^n for type N, T and C errors. The output is the MSEs for p and v^j. Under initializations, the dimension, the number of independent variables k and the type of M-type estimator can be changed.

simulation_experiments_2.R contains the code for the simulation experiments on the sphere S^n and hyperbolic space \mathbb{H}^n for various values of sigma. The output is the sample variance for p. Under initializations, the dimension, the number of independent variables k and the type of M-type estimator can be changed.

kendall_experiments.R contains the code for the experiments on Kendall's 2-dimensional shape space. The type of M-type estimator can be changed under the initializations section. We have used the preprocessed data provided by Cornea and Zhu on their website http://www.bios.unc.edu/research/bias/software.html under the title 'Regression Models on Riemannian Symmetric Spaces'.
