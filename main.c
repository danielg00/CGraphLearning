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

#define ITERS 10

DAG *init_graph(matrix *data);
int find_best_mod(DAG *G, int *option, double (*scoreFunc)(vertex *));
int best_mod_if_connected(vertex *v1, vertex *v2, int *option, double (*scoreFunc)(vertex *), double *max_score_diff);

// SEEMS MOST IS OK, CURRENT BUG IS THAT I CANT CALL BIC_SCORE MORE THAN ONCE.
// THINK IT HAS SOMETHING TO DO WITH COPYING MATRIX DURING INVERSION, SINCE INVERSION USED TO WORK FIND,
// ALTHOUGH WHEN TESTING I ONLY CALLED THE FUNCTION ONCE.

// ANOTHER ONE IS THAT DELETE EDGE STILL GETTING CALLED FOR NON EXISTENT EDGES

// TODO FIGURE OUT COPY ARRAY PART, MOST LIKELY VAR_OF_RES MOST LIKELY NEEDS complete rework.

int main()
{
    char * fname = "test_data/test_d.npy";
    
    matrix *A = load_matrix_from_file(fname);  // features x samples
    
    DAG *G = init_graph(A);

    add_child(&(G->nodes[0]), &(G->nodes[1]));
    add_child(&(G->nodes[2]), &(G->nodes[1]));
    double var1 = BIC_score(&(G->nodes[1]));
    
    add_child(&(G->nodes[3]), &(G->nodes[1]));
    
    double var2 = BIC_score(&(G->nodes[1]));
    
    return 0;
}



int main2()
{
    char * fname = "test_data/test_d.npy";
    
    matrix *A = load_matrix_from_file(fname);  // features x samples
    
    DAG *G = init_graph(A);
    
    int opt[3];  // opt = {v1, v2, mod_type}
    int mod;
    int i = 0;
    
    srand(time(NULL));
    
    while (i < 10)
	{
	    printf("\n ############## EPOCH %d ##############\n", i);
	    
	    mod = find_best_mod(G, opt, &BIC_score);
	    printf(" \n HERE \n");
	    printf(" (%d, %d, %d) ", opt[0], opt[1], opt[2]);
	    if (mod == 0)
		{
		    printf("CONVERGED");
		    return 0;
		}
	    
	    // Applying the best option
	    if (opt[2] == 0) // Apply edge deletion
		{
		    delete_edge(&(G->nodes[opt[0]]), &(G->nodes[opt[1]]));
		    printf("With score of %f\n", BIC_score( &(G->nodes[opt[0]]) ) + BIC_score( &(G->nodes[opt[1]]) ) );
		}
	    
	    else
		{
		    add_child(&(G->nodes[opt[0]]), &(G->nodes[opt[1]]));		    
		    printf("Total mods: %d, Best is add edge ( %d ) --> ( %d ) ", mod, opt[0], opt[1]);
		    printf("With score of %f\n", BIC_score( &(G->nodes[opt[0]]) ) + BIC_score( &(G->nodes[opt[1]]) ) );
		}
	    
	    print_graph(G);

	    i++;
	}
    
    print_graph(G);
    freeMatrix(A);	
    return 0;
}



int find_best_mod(DAG *G, int *option, double (*scoreFunc)(vertex *))  // Returns 0 if no best mod was found
{
    

    // We're basically looking through all possible connections of G,
    // making all possible modifications of an edge, computing the best score, and then
    // undoing it.
    //
    // Since scoreFunc is decomposable along the vertices, we just need to compute
    // the score of unmodified vertices, and see if the modification increase it.
    // The largest increase in local score between any two vertices are the vertices
    // chosen with the corresponding modification.
    // We want to avoid checking for cycles, calling is_child, and calling scoreFunc
    // too much because these are relatively expensive.

    //
    // Modifies options to specify information for the best modification type:
    //    options = {index1, index2, mod_type}
    //    Where index1 index2 are the indices specifying vertices v1, v2, whose edge is to be modified.
    //    And where mod_type == 1: if v1 --> v2 is best option.
    //              mod_type == 0: if v1 v2 is to have no connecting edge.

    
    double max_score_diff = 0.;  // The modification with the largest value score_diff is the 'chosen one'.
    int num_mods = 0.; // Number of possible modifications. If this is zero, then no better mod was found.

    for (int i = 0; i < G->num_nodes-1; i++)
	{
	    for (int j = i+1; j < G->num_nodes; j++)
		{
		    vertex *v1 = &(G->nodes[i]); vertex *v2 = &(G->nodes[j]);
		    // Two vertices can either be connected or not-connected,
		    // Calling is_child is O(n) so we want to avoid constantly doing that.
		    // We also want to skip over cases that cause cycles.
		    
		    if (is_child(v1, v2))
			{
			    num_mods += best_mod_if_connected(v1, v2, option, scoreFunc, &max_score_diff); 
			}

		    
		    else if(is_child(v2, v1))
			{
			    num_mods +=  best_mod_if_connected(v2, v1, option, scoreFunc, &max_score_diff);
			}


		    // No edge between v1 and v2;
		    else  
			{
			    double old_score;
			    double diff_for_12 = max_score_diff;
			    double diff_for_21 = max_score_diff;
			    
			    old_score = (*scoreFunc)(v1) + (*scoreFunc)(v2);

			    if(!check_if_path(v1, v1, v2))
				{
				    add_child(v2, v1);
				    diff_for_12 = (*scoreFunc)(v1) + (*scoreFunc)(v2) - old_score;
				    delete_edge(v2, v1);    
				}

			    if(!check_if_path(v2, v2, v1))
				{
				    add_child(v1, v2);
				    diff_for_21 = (*scoreFunc)(v1) + (*scoreFunc)(v2) - old_score;
				    delete_edge(v1, v2);
				}
			    
			    else  { continue; }  // Adding an edge either way causes a cycle so continue.
			    

			    if (diff_for_12 <= max_score_diff && diff_for_21 <= max_score_diff)  // Neither changes are better.
				{
				    continue;
				}
			    else if (diff_for_12 > max_score_diff)  {  // v2 --> v1 is best choice.
				num_mods += 1; 
				option[0] = i; option[1] = j; option[2] = 1;
				max_score_diff = diff_for_12;
			    }

			    else {  // v1 --> v2 is best choice.
				num_mods += 1;
				option[0] = j; option[1] = i; option[2] = 1;
				max_score_diff = diff_for_21;
			    }
			    
			}
		}
	}

    return num_mods;
}


int best_mod_if_connected(vertex *v1, vertex *v2, int *option, double (*scoreFunc)(vertex *), double *max_score_diff)
{
    // Takes two vertices, and tries to find a config that is better than current max_score_diff.
    // The local score is the change that increases the score the most.
    // If the local change increase the score than max_score_diff, then tentatively that change is adopted.
    // Otherwise, we continue onto next vertex pair.
    
    // Returns 0 if no better config was found.

    double diff_for_delete, diff_for_reverse;

    double old_score = (*scoreFunc)(v1) + (*scoreFunc)(v2);
	
    delete_edge(v1, v2);
    diff_for_delete = (*scoreFunc)(v2) + (*scoreFunc)(v1) - old_score;
    
    // If a path between v1 and v2 exists after deleting, then reversing the edge causes a  cycle.
    
    if (!check_if_path(v1 ,v1, v2)) // If there's no path, we can proceed to check score for reversed.
	{	    
	    add_child(v2, v1);
	    diff_for_reverse = (*scoreFunc)(v2) + (*scoreFunc)(v1) - old_score;
	    
	    if (diff_for_reverse < *max_score_diff && diff_for_delete < *max_score_diff)  // Doing nothing was best option.
		{
		    delete_edge(v2, v1); add_child(v1, v2);
		    return 0;
		}
	    
	    else if (diff_for_delete <= diff_for_reverse )  // Reversing the edge was best option.
		{
		    *max_score_diff =diff_for_reverse;
		    option[0] = v2->id; option[1] = v1->id; option[2] = 1;
		    
		    delete_edge(v2, v1); add_child(v1, v2);
		    return 1;
		}
	    
	    else // deleting the edge was best option.
		{
		    *max_score_diff = diff_for_delete;
		    option[0] = v1->id; option[1] = v2->id; option[2] = 0;
		    add_child(v1, v2);
		    return 1;
		}
	}
    
    
    else  // Reversing the edge isn't an option, so we just compare with delete edge.
	{
	    add_child(v1, v2);
		
		if (diff_for_delete > *max_score_diff)
		    {
			*max_score_diff = diff_for_delete;
			option[0] = v1->id; option[1] = v2->id; option[2] = 0;
			return 1;
		    }
	    
		else  
		    {
			return 0;
		    }
	}
    
    return 0;
}




DAG *init_graph(matrix *data)  // Each row of data are observations for a feature. 
{
    int num_nodes = data->dims[0];
    DAG *G = malloc(sizeof(*G));
    G->nodes = malloc(num_nodes*sizeof(vertex));
    G->num_nodes = num_nodes;
    for (int i = 0; i < num_nodes; i++)
	{
	    G->nodes[i].id = i;

	    G->nodes[i].num_parents = 0;
	    G->nodes[i].parents = calloc(num_nodes, sizeof(vertex*));
	    
	    G->nodes[i].num_children = 0;
	    G->nodes[i].children = calloc(num_nodes, sizeof(vertex*)); // 2X Maximum size of graph; REALLY MEMORY INEFFICIENT.

	    G->nodes[i].data = data->data[i];
	    G->nodes[i].num_samples = data->dims[1];

	    G->nodes[i].calc_score_needed = 1;
		
	}
    return G;
}


