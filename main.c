#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "linalg.h"
#include "io.h"
#include "graph.h"



/* typedef struct { */
/*     int * dims; */
/*     double ** data;     */
/* } matrix; */
#define ITERS = 100

int main() {
    char * fname = "test_data/test_d.npy";
    
    matrix *A = load_matrix_from_file(fname);  // features x samples

    DAG *G = init_graph(A);

    double c_score = 0;
    double s1, s2;
    vertex *contender_nodes[2];
    int c;
    for (int j = 0, i < ITERS; i++)
	{
	    i = c % G->num_nodes;
	    s1, s2 = 0, 0;
	    c_score = score_graph(G, A);
	    contender_nodes = choose_empty_nodes(G);
	    n1, n2 = contender_nodes[0], contender_nodes[0];
	    
	    add_child(n1, n2);
	    s1 = score_graph(G, A);
	    
	    if (check_cyclic(G) == 1)
		{
		    delete_child(n1, n2);
		    s1 = -1;
		}
	    
	    delete_child(n1, n2);

	    add_child(n2, n1);
	    s1 = score_graph(G, A);
	    delete_child(n2, n1);


	    // score_if_add_child { add_child  score  delete_child } 
	    
	    
	    c++;
	}	    
	    
    freeMatrix(A);	
    return 0;
}
