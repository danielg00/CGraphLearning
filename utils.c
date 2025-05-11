#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#include "linalg.h"
#include "io.h"
#include "graph.h"
#include "score_functions.h"


matrix *matmul(matrix * A, matrix * B)
{
    assert(A->dims[1] == B->dims[0] && "Matrix dimensions are incompatabilble.");

    matrix *C = malloc(sizeof(*C));
    C->dims[0] = A->dims[0]; C->dims[1] = A->dims[1];
				 
    C->data = malloc(C->dims[0]*sizeof(double *));
    for (int row = 0; row < A->dims[0]; row++)
	{
	    C->data[row] = calloc(C->dims[1], sizeof(double));
	    for (int i = 0; i < A->dims[1]; i++)
		{
		    for (int j = 0; j < B->dims[0]; j++)
			{
			    C->data[row][j] += A->data[row][i]*B->data[i][j];
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
    printf("\n ptr = %d = \n", A->dims[0]);
    A->data = malloc(A->dims[0] * sizeof(double *));
    for (int i = 0; i < A->dims[0]; i++)
	{
	    A->data[i] = calloc(A->dims[1], sizeof(double));
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

void arraycpy(double **Copy, double **A, int *dims)
{
    for (int i = 0; i < dims[0]; i++)
	{
	    memcpy(Copy[i], A[i], dims[1] * sizeof(double));
	}
}


	    

void free_graph(DAG *G)
{
    int D = G->num_nodes;
    for (int i = 0; i < D; i++)
	{
	    free(G->nodes[i].children);
	    free(G->nodes[i].parents);

	}
    free(G->nodes);
    free(G);
}


void print_graph(DAG *G)
{
    printf("G: \n"); 

    for (int i = 0; i < G->num_nodes; i++)
	{
	    printf(" ( %d ) ====> (", i);
	    for (int j = 0; j < G->nodes[i].num_children; j++)
		{
		    printf(" %d ", G->nodes[i].children[j]->id);
		}
	    printf(")\n"); 
	}
}
