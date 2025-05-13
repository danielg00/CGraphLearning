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
    C->dims = malloc(2 * sizeof(int));
    C->dims[0] = A->dims[0]; C->dims[1] = B->dims[1];

    C->data = malloc(C->dims[0] * sizeof(double *));
    for (int i = 0; i < A->dims[0]; i++)
	{
	    C->data[i] = calloc(C->dims[1], sizeof(double));
	    
	    for (int j = 0; j < B->dims[0]; j++)
		{
		    for (int k = 0; k < B->dims[1]; k++)
			{
			    C->data[i][k] += A->data[i][j]*B->data[j][k];
			}
		}
	}
    return C;
}

    
void alloc_array(matrix *A)
{
    A->data = malloc(A->dims[0] * sizeof(double *));
    for (int i = 0; i < A->dims[0]; i++)
	{
	    A->data[i] = calloc(A->dims[1], sizeof(double));
	}
}


void freeArray(double ** array, int dim0) // Frees each row first and then frees initial pointer array.
{
    for (int i = 0; i < dim0; i++)
	{
	    free(array[i]);
	}
    free(array);
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
