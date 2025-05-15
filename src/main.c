#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "linalg.h"
#include "io.h"
#include "graph.h"
#include "score_functions.h"
#include "utils.h"

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
	    printf("\n (%d, %d, %d)", opt[0],  opt[1], opt[2]);

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
	    print_graph(G);

	    i++;
	}
    
    print_graph(G);

    free_graph(G);
    freeMatrix(A);
    return 0;
}
