# TODO
Things that require somewhat immediate attention:

	(1) linalg.c Fix/merge together operations in variance_of_residuals.
	(2) Serpate VOR into seperate residual computation and linear regression functions. 
	(2) io.c: Support float 32
	(3) io.c Support verbose mode.
	(3) io.c: write adjency or python list to file
	(4) score_functions.c: Support simulated annealing or other
	(5) main.c: Terminate if cycle in moves found. 
	(6) Support better binary file than .npy.
	
Plans for the future:
	
	(2) Use cholesky decomp for linear regession than LU decomp.
	(1) Implement multithreading for linear regression, and graph search.
	(2) Implement different score functions. 
	(3) implement probabilitic choice. 
	(4) Better CMD interface. 
	(5) Better modularity to support (a) score function variety (b) prediction function variety for
	different data distributions, and (c) search algorthim variety.
	(6) Test functions.
		
	
	
	
