#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "linalg.h"


matrix * invert_matrix(matrix * A)
{
    /*
      Constructs L and U arrays from LU-decomposition and solves for the inverse
      using backward and forward substitution.

      Note, that LU_decomposition modifies A in place.

      
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



int * LU_decomposition(matrix * A, double *** ptrL, double *** ptrU)
{
    /* LU decomposition with pivoting. (see https://en.wikipedia.org/wiki/LU_decomposition)
       
       
       Function is used as way to invert a matrix. 

       TODO:
       - Check if matrix is singular and check if it is square, make return value specify this if true.
    */

    int D = A->dims[0];

    *ptrL = malloc(D*sizeof(double *));   
    *ptrU = malloc(D*sizeof(double *));  
    
    for (int i = 0; i < D; i++)
	{
	    (*ptrL)[i] = calloc(D, sizeof(double));                           // Learned the hard way that its safer to make them square.
	    (*ptrL)[i][i] = 1.0;	                                      // initialise lower triangular matrix with identity on diagonal.

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
	    pivots[n] = pivots[max_val_index];
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



void solve_Ax_b(double **A, double **b, int D, int upper)
// Solves a system of linear equations, modidies b
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



void freeArray(double ** array, int dim0) // Frees each row first and then frees initial pointer array.
{
    for (int i = 0; i < dim0; i++)
	{
	    free(array[i]);
	}
    free(array);
}

void alloc_array(matrix *A)
{
    A->data = malloc(A->dims[0] * sizeof(double *));
    for (int i = 0; i < dims[0]; i++)
	{
	    A->data[i] = malloc(A->dims[1] * sizeof(double));
	}
}

void freeMatrix(matrix *M)
{
    freeArray(M->data, M->dims[0]);
    free(M->dims);
    free(M);
}

void printArray(double **Array, int dim0, int dim1)
{
    printf("\n====================\n");
    for (int i = 0; i < dim0; i++)
	{
	    for (int j = 0; j < dim1; j++)
		{
		    printf(" %f ", Array[i][j]);
		}
	    printf(" \n ");
	}
    printf("\n====================\n");

}


double variance_of_residuals(matrix *X, double *Y)
{
    // First calculate regression coefficients using [(X X.T)^-1]X.Ty
    // X X.T is symmetric so transpose of inverse(X X.T) doesnt matter

    matrix * xtx - malloc(sizeof(*B));
    xtx->dims = malloc(2 * sizeof(int));
    alloc_array(xtx);
    
    xtx->dims[0] = X->dims[0];
    xtx->dims[1] = X->dims[0];

    double dot;
    for (int i = 0; i < (X->dims[0]/2)+1; i++)  // Gram matrix; dotting every row with eachother.
	{
	    for (int j = 0; j < X->dims[0]; j++)
		dot = 0;
		{
		    for (int k = 0; k < X->dims[1]; k++)
			{
			    dot += X->data[i][k]*X->data[j][k];
			}
		    xtx->data[i][j] = dot;
		    xtx->data[j][i] = dot;
		}
	}
    matrix *NEXT = inverse(xtx);
    freeMatrix(xtx);

    // FINISH
}
    
	    
		
