#include <stdio.h>
#include <stdlib.h>
#include "linalg.h"


double *** LU_decomposition(matrix * A)
{
    /* LU decomposition with pivoting.
       Function is used as way to invert matrix, so (1) matrix `A` will be altered but we
       won't need it after function call anyways, and (2) `L` and `U` aren't filled with zeros
       but are a sequence of pointers to arrays of shrinking/growing size.

       TODO:
       Implement pivots
       Check if matrix is singular and check if it is square.
      
     */
    
    double ** L;
    double ** U;
    int D = A->dims[0];

    L = malloc(D*sizeof(double *));   // initialise lower triangular matrix.
    U = malloc(D*sizeof(double *));   // initialise upper triangular matrix.
    for (int i = 0; i < D; i++)
	{
	    L[i] = malloc(D*sizeof(double *));
	    L[i][i] = 1;

	    U[i] = malloc((D-i)*sizeof(double *));
	}

    for (int n = 0; n < D-1; n++)
	{
	    U[n][n] = A->data[n][n];
	    
	    for (int i = n + 1; i < D; i++)
		{
		    L[i][n] = (1/A->data[n][n])*A->data[i][n];   // L[n+1:, n] = l
		    U[n][i] = A->data[n][i];  // U[n, n+1:] = u
		}

	    // A = A - outer_product(l, u)
	    for (int i = n + 1; i < D; i++)
		{
		for (int j = n + 1; j < D; j++)
		    {
    		        A->data[i][j] = A->data[i][j] - L[i][n]*U[n][j];
		    }
		}
	}
    U[D-1][D-1] = A->data[D-1][D-1];
    
    double *** return_array = malloc(2*sizeof(double**));
    return_array[0] = L;
    return_array[1] = U;
    return return_array;
 }
