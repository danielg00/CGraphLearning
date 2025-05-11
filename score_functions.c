#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#include "linalg.h"
#include "graph.h"
#include "score_functions.h"

double BIC_score(vertex *v)
{
    // Let N be the number of samples and s be the variance of the residuals of vertex v given its P parents.
    // Then BIC(v) = -(N/2)*log(s) - (P+2)/2  *log(N)
    int P = v->num_parents;
    int N = v->num_samples;

    /* if (!(v->score_changed)) */
    /* 	{ */
    /* 	    return v->score; */
    /* 	}     */
    
    if (P == 0)
	{
	    return log(N);
	}
    
    matrix *X = malloc(sizeof(*X));
    X->dims = malloc(2  * sizeof(int));
    X->dims[0] = P;
    X->dims[1] = N;
    
    X->data = calloc(P, sizeof(double*));

    printf("\n HERE2 %d\n ");
	
    for (int i = 0; i < P; i++)
	{
	    X->data[i] = v->parents[i]->data;
	}
    
    double var = variance_of_residuals(X, v->data);
    
    free(X->dims);
    free(X->data);
    X->data = NULL;
    free(X);
    
    double score = (-N/2)*log(var) - log(N)*(P+1)/2;
    
    printf(" \n HERE %d\n");
    
    return score;
}



