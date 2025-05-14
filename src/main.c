#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "linalg.h"
#include "io.h"
#include "graph.h"
#include "score_functions.h"
#include "utils.h"


/* int main() */
/* { */
/*     matrix *A = load_matrix_from_file("test_data/test_structured.npy");  // features x samples */
/*     DAG *G = init_graph(A); */
/*     add_child(&(G->nodes[1]), &(G->nodes[0])); */
/*     add_child(&(G->nodes[2]), &(G->nodes[0])); */
/*     add_child(&(G->nodes[3]), &(G->nodes[0])); */
/*     vertex *v = &(G->nodes[0]); */
				      
/*     int P = v->num_parents; */
/*     int N = v->num_samples; */

/*     matrix *X = malloc(sizeof(*X)); */
/*     X->dims = malloc(2  * sizeof(int)); */
/*     X->dims[0] = P; */
/*     X->dims[1] = N; */
    
/*     X->data = malloc(P * sizeof(double*)); */
	
/*     for (int i = 0; i < P; i++) */
/* 	{ */
/* 	    X->data[i] = v->parents[i]->data; // This simply points to the original data, we don't want to free/modify what it points to.  */
/* 	} */

/*     double betas[X->dims[0]]; double intercept; */
/*     linear_regression(X, v->data, betas, &intercept); */

/*     for (int i = 0; i < P; i++) { printf(" b%f ", betas[i]);} */
/*     return 0; */
/* } */

int main(int arvc, char *argv[])
{   
    char *p;
    char * fname = argv[1];
    int MAX_ITERS = strtol(argv[2], &p, 10);
    
    matrix *A = load_matrix_from_file(fname);  // features x samples
    
    DAG *G = init_graph(A);
    
    int opt[3];  // opt = {v1, v2, mod_type}
    int mod = 1;
    int i = 0;
    
    while (i < MAX_ITERS || mod == 0)
	{
	    printf("\nITERATION %d", i);
	    
	    mod = find_best_mod(G, opt, &BIC_score);
	    printf("\n   Modications found: %d", mod);

	    if (mod == 0)
		{
		    printf("\n   FINAL SCORE: %f \n", graph_score(G, *BIC_score));
		    printf("\n FINAL GRAPH: "); print_graph(G);
			
		    printf("\nCONVERGED \n");
		    return 0;
		}
	    
	    // Applying the best option
	    if (opt[2] == 0) // Apply edge deletion
		{
		    delete_edge(&(G->nodes[opt[0]]), &(G->nodes[opt[1]]));
		}

	    else  // else add or reverse and edge.
		{
		    if(opt[2] == -1) // reverse the edge
			{
			    delete_edge(&(G->nodes[opt[0]]), &(G->nodes[opt[1]]));
			    add_child(&(G->nodes[opt[1]]), &(G->nodes[opt[0]]));
			}
		    
		    else  // Add and edge
			{
			    add_child(&(G->nodes[opt[0]]), &(G->nodes[opt[1]]));
			    /* printf("Total mods: %d, Best is add edge ( %d ) --> ( %d ) ", mod, opt[0], opt[1]); */
			    /* printf("With score of %f\n", BIC_score( &(G->nodes[opt[0]]) ) + BIC_score( &(G->nodes[opt[1]]) ) ); */
			}
		}
	    
	    printf("\n   CURRENT SCORE: %f\n", graph_score(G, *BIC_score));

	    i++;
	}
    
    print_graph(G);

    free_graph(G);
    freeMatrix(A);	
    return 0;
}
