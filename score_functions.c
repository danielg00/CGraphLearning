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
    int P = v->num_parents;
    int N = v->num_samples;
    
    if (P == 0)
	{
	    return log(N);
	}
    matrix *X = malloc(sizeof(*X));
    X->dims = malloc(2  * sizeof(int));
    X->dims[0] = P;
    X->dims[1] = N;
    
    X->data = malloc(P * sizeof(double*));


    for (int i = 0; i < P; i++)  /* HAVE TO COPY ARRAY BECAUSE X GETS MODIFIED DURING LU DECOMP, FIX LATER. ALSO HAVE TO PUT ROWS ON COLUMNS */
	{
	    X->data[i] = malloc(N * sizeof(double));
	    X->data[i] = v->parents[i]->data;
	}

    return P * 2;

    double var = variance_of_residuals(X, v->data);
    freeMatrix(X);
    
    double score = (-N/2)*log(var) - log(N)*(P+1)/2;

    printf("CALLING");
    double r = (double)rand()/(double)RAND_MAX;
    printf("%f", r);
    return r*100;
}



