#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "linalg.h"
#include "graph.h"
#include "score_functions.h"


double gaussian_BIC_score(DAG *G, matrix *A)
{
    /* Computes the BIC score of a graph
       Iterates over each node and computes the local bic score of that node given its immediate parents
       WE assume that the rows of A give obeservations 
     */
    double score = 0;
    int D = G->num_nodes;
    for (int i = 0; i < D; i++)
	{
	    score += node_bic_score(G->nodes[i], matrix *A);
	}
}


double node_bic_score(DAG *G, matrix *A)
{
    

