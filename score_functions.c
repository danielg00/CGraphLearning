#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "linalg.h"
#include "graph.h"
#include "score_functions.h"

double BIC_score(vertex *v)
{
    // Let N be the number of samples and s be the variance of the residuals of vertex v given its P parents.
    // Then BIC(v) = -(N/2)*log(s) - (P+2)/2  *log(N)
    /* int P = v->num_parents; */
    /* int N = v->num_samples; */
    
    /* matrix *X = malloc(sizeof(*X)); */
    /* X->data = malloc(P * sizeof(double*));  // We don't have to transpose anything in var_of_res if we do P features x N samples */
    /* X->dims = malloc(2  * sizeof(int)); */
    
    /* X->dims[0] = P; */
    /* X->dims[1] = N; */
    /* printf("score functions.c P: \nHERE!\n = %d =", P); */
    /* for (int i = 0; i < P; i++)  // HAVE TO COPY ARRAY BECAUSE X GETS MODIFIED, FIX LATER. ALSO HAVE TO PUT ROWS ON COLUMNS */
    /* 	{ */
    /* 	    X->data[i] = malloc(N * sizeof(double)); */
    /* 	    for (int j = 0; j < P; j++) */
    /* 		{ */
    /* 		    X->data[i] = v->parents[i]->data; */
    /* 		} */
    /* 	} */

    /* double var = variance_of_residuals(X, v->data); */
    /* freeMatrix(X); */
    
    /* double score = (-N/2)*log(var) - log(N)*(P+1)/2; */
    /* printf("%f", score); */
    /* printf("\nHERE!\n"); */
    double r = ((double)rand()/(double)RAND_MAX);
    /* printf("%f", r); */
    return r ;
}
