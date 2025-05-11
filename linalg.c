#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#include "linalg.h"

extern matrix *SCRATCH;

matrix *invert_matrix(matrix * A)
{
    /*
      Constructs L and U arrays from LU-decomposition and solves for the inverse
      using backward and forward substitution.

      Arguments: square matrix A dimension N to be inverted.
      
      Creates: (1) double ** arrays L, U of same dimension N.
               (2) Inverted square matrix B of dimension N.
	       (3) int array pivots of length N.
	       (4) copy of matrix A

      Frees:   (1) arrays L, U.
               (2) array pivots
      
      Returns: Inverted matrix B 

      
     */
    
    double **L, **U;
    int *pivots;
    
    matrix *B = malloc(sizeof(*B));
    
    B->data = malloc(A->dims[0] * sizeof(double *));
    B->dims = malloc(2 * sizeof(int));
    memcpy(B->dims, A->dims, 2 * sizeof(int));
    
    for (int i = 0; i < A->dims[0]; i++)
	{
	    B->data[i] = calloc(A->dims[0], sizeof(double));
	}

    int idx;
    pivots = LU_decomposition(A, &L, &U);
    for (int i = 0; i < A->dims[0]; i++)                                      // Transpose/inverse of permutation matrix
	{
	    idx = pivots[i];
	    B->data[idx][i] = 1.;
	}

    for (int i = 0; i < A->dims[0]; i++)                                      // Solving with each row of P.T matrix
	{
	    solve_Ax_b(L, &(B->data[i]), A->dims[0], 0);
	    solve_Ax_b(U, &(B->data[i]),  A->dims[0], 1);
	}
    
    freeArray(L, A->dims[0]);
    freeArray(U, A->dims[0]);
    
    free(pivots);
    return B;  // !! B is transposed !!
}



int * LU_decomposition(matrix *copy_A, double *** ptrL, double *** ptrU)
{
    /* LU decomposition with pivoting. (see https://en.wikipedia.org/wiki/LU_decomposition)
       
       
       Function is used as way to invert a matrix. 

       TODO:
       - Check if matrix is singular and check if it is square, make return value specify this if true.

       Creates: (1) matrix A, a copy of matrix Ad
                (2) 
    */
    int D = copy_A->dims[0];
    *ptrL = malloc(D*sizeof(double *)), *ptrU = malloc(D*sizeof(double *));
    int *pivots;
    
    for (int i = 0; i < D; i++)
	{
	    (*ptrL)[i] = calloc(D, sizeof(double));                           // Learned the hard way that its safer to make them square.
	    (*ptrL)[i][i] = 1.0;	                                      // initialise lower triangular matrix with identity on diagonal.

	    (*ptrU)[i] = calloc(D, sizeof(double));
	}

    
    /* ======== PIVOTING ======= */
    pivots = malloc(D * sizeof(int));
    for (int i = 0; i < D; i++)
	{
	    pivots[i] = i;
	}

    int j, max_val_index;
    for (int n = 0; n < D-1; n++)
	{
	    max_val_index = n;
	    double max_value = fabs(copy_A->data[n][n]);

	    for (int i = n; i < D; i ++) 	                              // Finding max value of column n, starting from row n.
		{
		    if (fabs(copy_A->data[i][n]) > max_value)
			{
			    max_value = fabs(copy_A->data[i][n]);
			    max_val_index = i;
			}
		}
	    // Swapping rows
	    double *ptrRow = copy_A->data[max_val_index];
	    copy_A->data[max_val_index] = copy_A->data[n];
	    copy_A->data[n] = ptrRow;

	    
	    j = pivots[n];
	    pivots[n] = pivots[max_val_index];
	    pivots[max_val_index] = j;
	}
   

    /* ====== LU decomposition - Doolittle's algorithm. ======*/
    for (int n = 0; n < D; n++)
	{
	    (*ptrU)[n][n] = copy_A->data[n][n];
	    for (int i = n + 1; i < D; i++)
		{
		    (*ptrL)[i][n] = (1/copy_A->data[n][n])*copy_A->data[i][n];          // L[n+1:, n] = l
		    (*ptrU)[n][i] = copy_A->data[n][i];                            // U[n, n+1:] = u
		}
	    // A = A - outer_product(l, u)
	    for (int i = n + 1; i < D; i++)
		{
		    for (int j = n + 1; j < D; j++)
			{
			    copy_A->data[i][j] -= (*ptrL)[i][n]*(*ptrU)[n][j];
			}
		}
	}

    return pivots;
}



void solve_Ax_b(double **A, double **b, int D, int upper)
// Solves a system of linear equations; modifies b
{
    double sum;
    switch (upper)
	{
	case 0:  // Solve with A as a lower triangular matrix with unit diagonal.
	    for (int row = 0; row < D; row++)
		{
		    sum = 0;
		    for (int col = 0; col < row; col++)
			{
			    sum += (*b)[col]*A[row][col];
			}
		    (*b)[row] -= sum; // L has unit diag so no need to divide.

		}
	    break;

	case 1:  // Solve with A as an upper triangular matrix.
	    for (int row = D-1; row >= 0; row--)
		{
		    sum = 0;
		    for (int col = row + 1; col < D; col++) 
			{
			    sum += (*b)[col]*A[row][col];
			}
		    (*b)[row] = ((*b)[row] - sum)/(A[row][row]);
		}
	    break;
	}
}


		
double variance_of_residuals(matrix *X, double *Y)
{
    // First calculate regression coefficients using [(X.T X)^-1] X.T y
    // X X.T is symmetric so transpose of inverse(X.T X) doesnt matter
    // Samples are row wise not column-wise.
    // This is basically the MSE-error

    matrix *Inv, *xtx, *C;

    printf(" ASDASDAS " );
    xtx = malloc(sizeof(*xtx));
    xtx->dims = malloc(2 * sizeof(int));
    xtx->dims[0] = X->dims[0]; xtx->dims[1] = X->dims[0];
	
    alloc_array(xtx);

    
    double dot;
    for (int i = 0; i < X->dims[0]; i++)  // Gram matrix; dotting every row with eachother.
	{
	    for (int j = i; j < X->dims[0]; j++)
		{
		    dot = 0;
		    for (int k = 0; k < X->dims[1]; k++)
			{
			    xtx->data[i][j] += X->data[i][k]*X->data[j][k];
			    xtx->data[j][i] += X->data[i][k]*X->data[j][k];
			}
		}
	}

    printf(" \n BEFORE INVERT MAT \n");
    
    Inv = invert_matrix(xtx);
    printf("( %d, %d )\n", Inv->dims[0], Inv->dims[1]);
    freeMatrix(xtx);

    C = matmul(Inv, X);
    
    freeMatrix(Inv);

    
    double betas[X->dims[0]];  // Slopes
    for (int i = 0; i < X->dims[0]; i++)
	{
	    dot = 0;
	    for (int j = 0; j < X->dims[1]; j++)
		{
		    dot += Y[j]*C->data[i][j];
		}
	    betas[i] = dot;
	}
    
    // Calculating intercept. / mean of Y  - sample mean of each feature times its slope
    double mu_Y = 0;
    for (int i = 0; i < X->dims[1]; i++) {mu_Y += Y[i]; }
    mu_Y /= X->dims[1];

    double mu_XB = 0;
    for (int i = 0; i < X->dims[0]; i++)
	{
	    for (int j = 0; j < X->dims[1]; j++)
		{
		    mu_XB += X->data[i][j];  // Sum of samples of feature j.
		}
	    mu_XB *= betas[i]/X->dims[1];
	}
    double intercept = mu_Y - mu_XB;

    // Calculate the variance of residuals.
    double VOR;
    VOR = 0;
    for (int i = 0; i < X->dims[1]; i++)
	{
	    dot = 0;
	    for (int j = 0; j < X->dims[0]; j++)
		{
		    dot += (X->data[j][i] * betas[j]) + intercept;
		}
	    VOR += (Y[i] - dot)*(Y[i] - dot);
	}
    
    VOR /= (X->dims[1] + 1);
	
    return VOR;

}
    
