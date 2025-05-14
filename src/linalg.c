#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "linalg.h"
#include "graph.h"

matrix *invert_matrix(matrix * A)
{
    /*
      Constructs L and U arrays from LU-decomposition and solves for the inverse
      using backward and forward substitution.

      `Inv` will be transposed because it's more efficient to just copy the results to
      rows than copy it to the columns. This can ammended later to obtain
      the correct coefficients for linear regression.

      Arguments: square matrix A dimension N to be inverted.
      
      Creates: (1) double ** arrays `L`, `U` of same dimension N.
               (2) Inverted square matrix `B` of dimension N.
	       (3) int array `pivots` of length N.

      Frees:   (1) arrays `L`, `U`.
               (2) array `pivots`.
      
      Returns: (3) Inverted matrix B.
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
    /* LU decomposition with partial pivoting. (see https://en.wikipedia.org/wiki/LU_decomposition)
       Computes lower and upper triangular matrices, L and U, so that mathmul(L, U) = A

       Functions allocations memory for L, and U, and Pivots. Caller should free these when finished.  
    */

    int D = copy_A->dims[0];
    *ptrL = malloc(D*sizeof(double *)), *ptrU = malloc(D*sizeof(double *));
    int *pivots;
    
    for (int i = 0; i < D; i++)
	{
	    (*ptrL)[i] = calloc(D, sizeof(double));         
	    (*ptrL)[i][i] = 1.0;	        

	    (*ptrU)[i] = calloc(D, sizeof(double));
	}

    
    /* ====== PIVOTING ====== */
    // We swap rows so that the maximum absolute values lie on the diagonal.
    // The rows of A are permuted and the pivots are logged in `pivots` which is just a permutation map.
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

	    // logging swapping the corresponding indices in `pivots`.
	    j = pivots[n];
	    pivots[n] = pivots[max_val_index];
	    pivots[max_val_index] = j;
	}

    
    /* ====== LU decomposition - Doolittle's algorithm. ======*/
    for (int n = 0; n < D; n++)
	{
	    if (copy_A->data[n][n] < 1e-12)                               // matrix is singular
        {
            pivots[0] = -1;
            return pivots;
        }
	    (*ptrU)[n][n] = copy_A->data[n][n];
	    for (int i = n + 1; i < D; i++)
		{
		    (*ptrL)[i][n] = (1/copy_A->data[n][n])*copy_A->data[i][n];     // L[n+1:, n] = l
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
// Solves a system of linear equations; modifies b in place.
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

void linear_regression(matrix *X, double *Y, double *betas, double *intercept)
{
    // We first calculate the regression coefficients by solving B = [(X.T X)^-1] X.T y
    // We then get the variance of the residuals.

    // The proceeding code looks different than the expression above. This is mainly
    // because the inverse (X.T X^-1) is transposed and X is tranposed for speed.
    // We first obtain Z = `X`@`Y`, and then compute Z.T @ Inv. 

    matrix *Inv, *XtX;

    XtX = malloc(sizeof(*XtX));
    XtX->dims = malloc(2 * sizeof(int));
    XtX->dims[0] = X->dims[0];
    XtX->dims[1] = X->dims[0];
    
    alloc_array(XtX);

    // The expression X^tX is the dot product each column of X (a feature with N samples.)
    // with every other column giving an P x P matrix with index (i, j) being the dotproduct 
    ///of column i with column j. We have passed X in the form of P x N instead of N x P,
    // so we instead dot every row with everyother row.
    
    double dot;
    for (int i = 0; i < X->dims[0]; i++) 
	{
	    for (int j = i; j < X->dims[0]; j++)
		{
		    for (int k = 0; k < X->dims[1]; k++)
			{
			    XtX->data[i][j] += X->data[i][k]*X->data[j][k];
			    XtX->data[j][i] += X->data[i][k]*X->data[j][k];
			}
		}
	}
    
    Inv = invert_matrix(XtX);
    
    double XtY[X->dims[0]];
    memset(XtY, 0., X->dims[0] * sizeof(double));
    for (int i = 0; i < X->dims[0]; i++)
	{
	    for (int j = 0; j < X->dims[1]; j++)
		{
		    XtY[i] += X->data[i][j]*Y[j];
		}
	}

    // Our inverse is also tranposed, so we need to multiply XtY by the columns of Inv.
    // Its a bit more cache efficient iteratively build each element of Betas
    // so that we stay on the same row in Inv.
    
    /* double betas[X->dims[0]] */;  // our regression coefficients.
    memset(betas, 0., X->dims[0] * sizeof(double));
    for (int i = 0; i < Inv->dims[0]; i++)
	{
	    for (int j = 0; j < Inv->dims[1]; j++)
		{
		    betas[j] += XtY[i]*Inv->data[i][j];
		}
	}

    // Calculating intercept: mean(Y) - mean(Beta*X)
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

    *intercept = mu_Y - mu_XB;

    
    freeMatrix(Inv);
    freeMatrix(XtX);
}
		
double variance_of_residuals(matrix *X, double *Y)
{
    double betas[X->dims[0]]; double intercept;
    
    linear_regression(X, Y, betas, &intercept);
    // Calculate the variance of residuals.
    double VOR;
    VOR = 0;

    double dot;
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
    
