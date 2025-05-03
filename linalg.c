#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "linalg.h"


int * LU_decomposition(matrix * A, double *** ptrL, double *** ptrU)
{
    /* LU decomposition with pivoting. (see https://en.wikipedia.org/wiki/LU_decomposition)
       
       
       Function is used as way to invert matrix, so (1) matrix `A` will be altered but we
       won't need it after function call anyways, and (2) `L` and `U` aren't filled with zeros
       but are a sequence of pointers to arrays of shrinking/growing size.

       TODO:
       - Check if matrix is singular and check if it is square, make return value specify this if true.
    */

    int D = A->dims[0];

    *ptrL = malloc(D*sizeof(double *));   
    *ptrU = malloc(D*sizeof(double *));   // initialise upper triangular matrix.
    
    for (int i = 0; i < D; i++)
	{
	    (*ptrL)[i] = calloc(D, sizeof(double));                           // Learned the hard way that its safer to make them square.
	    (*ptrL)[i][i] = 1.0;	// initialise lower triangular matrix with identity on diagonal.

	    (*ptrU)[i] = calloc(D, sizeof(double)); 
	}

    
    /* ======== PIVOTING ======= */
    int *pivots = malloc(D * sizeof(int));
    for (int i = 0; i < D; i++)
	{
	    pivots[i] = i;
	}

    int j, max_val_index;
    for (int n = 0; n < D-1; n++)
	{
	    max_val_index = n;
	    double max_value = fabs(A->data[n][n]);

	    for (int i = n; i < D; i ++) 	                              // Finding max value of column n, starting from row n.
		{
		    if (fabs(A->data[i][n]) > max_value)
			{
			    max_value = fabs(A->data[i][n]);
			    max_val_index = i;
			}
		}
	    // Swapping rows
	    double *ptrRow = A->data[max_val_index];
	    A->data[max_val_index] = A->data[n];
	    A->data[n] = ptrRow;

	    
	    j = pivots[n];
	    pivots[n] = pivots[max_val_index]; /// <======= CAUSING ISSUES HRMMMMM 
	    pivots[max_val_index] = j;
	}
   

    /* ====== LU decomposition - Doolittle's algorithm. ======*/
    for (int n = 0; n < D; n++)
	{
	    (*ptrU)[n][n] = A->data[n][n];
	    for (int i = n + 1; i < D; i++)
		{
		    (*ptrL)[i][n] = (1/A->data[n][n])*A->data[i][n];          // L[n+1:, n] = l
		    (*ptrU)[n][i] = A->data[n][i];                            // U[n, n+1:] = u
		}
	    // A = A - outer_product(l, u)
	    for (int i = n + 1; i < D; i++)
		{
		    for (int j = n + 1; j < D; j++)
			{
			    A->data[i][j] -= (*ptrL)[i][n]*(*ptrU)[n][j];
			}
		}
	}
    
    return pivots;
}


matrix * invert_matrix(matrix * A)
{
    /*
      TAKES L and U arrays from LU decomposition and solves for the inverse. Not that the rows of L are permuted so it's
      not a lower triangular matrix in standard form.
      Let L, U be the decomp of matrix X and let Z be its inverse. Let I be the indentity matrix and P be the permutation matrix.
     */

    double **L, **U;
    int idx;
    int *pivots = LU_decomposition(A, &L, &U);

    
    matrix * B = malloc(sizeof(*B));
    
    B->data = malloc((A->dims[0])*sizeof(double *));
    B->dims = malloc(2*sizeof(int));
    memcpy((B->dims), (A->dims), 2*sizeof(int));

    
    for (int i = 0; i < A->dims[0]; i++)
	{
	    B->data[i] = calloc(A->dims[0], sizeof(double));
	}
	        
    for (int i = 0; i < A->dims[0]; i++)                                      // Transpose permutation matrix
	{
	    idx = pivots[i];
	    B->data[idx][i] = 1.;
	}

    for (int i = 0; i < A->dims[0]; i++)                                      // Solving with each row of P.T matrix
	{
	    solve_Ax_b(L, &(B->data[i]), A->dims[0], 0);
	    solve_Ax_b(U, &(B->data[i]),  A->dims[0], 1);
	}
    
    return B;  // !! B is transposed !!
}



void solve_Ax_b(double **A, double **b, int D, int upper)
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


double ** matmul(matrix * A, matrix * B)
{
    assert(A->dims[1] == B->dims[0] && "Matrix dimensions are incompatabilble.");
    
    double ** C = malloc(A->dims[0]*sizeof(double *));
    for (int row = 0; row < A->dims[0]; row++)
	{
	    C[row] = calloc(B->dims[1], sizeof(double));
	    for (int i = 0; i < A->dims[1]; i++)
		{
		    for (int j = 0; j < B->dims[0]; j++)
			{
			    C[row][j] += A->data[row][i]*B->data[i][j];
			}
		}
	}
	    
    return C;
}
