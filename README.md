# USAGE:

`make -B`

`./main path_to_data MAX_ITERS`

Note that the data should be in a .npy file with format n_features x n_samples (as opposed to the standard n_samples x n_features).

This is implementation of a hill climbing structure learning algorithm for directed gaussian graphical models. The project serves more as a way to brush up on skills
and isn't as fast nor as general as it could be (for instance I don't use any BLAS routines for matrix inversion, and I use LU decomposition instead of cholesky 
decomposition for linear regression). Nonetheless, I plan incorporate these things, along with a better CLI, and generalise it with better search algorithms. (see TODO.md). 
