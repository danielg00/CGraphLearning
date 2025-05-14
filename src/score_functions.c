#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "linalg.h"
#include "graph.h"
#include "score_functions.h"

double BIC_score(vertex *v)
{
    // Let N be the number of samples and var be the variance of the residuals of vertex v given its P parents.
    // Then BIC(v) = -(N/2)*log(var) - log(N)*(P+2)/2
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
	
    for (int i = 0; i < P; i++)
	{
	    X->data[i] = v->parents[i]->data; // This simply points to the original data, we don't want to free/modify what it points to. 
	}
    
    double var = variance_of_residuals(X, v->data);
    
    free(X->data);
    free(X->dims);
    free(X);
    
    double score = -(N/2)*log(var) - log(N)*(P+2)/2;
    return score;
}


double graph_score(DAG *G, double (*ScoreFunc)(vertex*))
{
    double score = 0;
    for (int i = 0; i < G->num_nodes; i++)
	{
	    score += (*ScoreFunc)(&(G->nodes[i]));
	}
    return score;
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
    //
    // Modifies options to specify information for the best modification type:
    //    options = {index1, index2, mod_type}
    //    Where index1 index2 are the indices specifying vertices v1, v2, whose edge is to be modified.
    //    And where mod_type == 1: if v1 --> v2 is best option.
    //              mod_type == 0: if v1 v2 is to have no connecting edge.
    //              mod_type == -1:if v2 --> v1 is the best option.

    
    double max_delta = 0.;  // The modification with the largest value score_diff is the 'chosen one'.
    int num_mods = 0.; // Number of possible modifications. If this is zero, then no better mod was found.

    for (int i = 0; i < G->num_nodes-1; i++)
	{
	    for (int j = i+1; j < G->num_nodes; j++)
		{
		    vertex *v1 = &(G->nodes[i]); vertex *v2 = &(G->nodes[j]);
		    
		    // Two vertices can either be connected or not-connected,
		    if (is_child(v1, v2))
			{
			    num_mods += best_mod_if_connected(v1, v2, option, scoreFunc, &max_delta); 
			}

		    else if(is_child(v2, v1))
			{
			    num_mods +=  best_mod_if_connected(v2, v1, option, scoreFunc, &max_delta);
			}


		    // No edge between v1 and v2;
		    else  
			{
			    double old_score;
			    double delta_for_12 = 0;
			    double delta_for_21 = 0;
			    
			    old_score = (*scoreFunc)(v1) + (*scoreFunc)(v2);

			    if(!check_if_path(v1, v1, v2))
				{
				    add_child(v2, v1);
				    delta_for_12 = (*scoreFunc)(v1) + (*scoreFunc)(v2) - old_score;
				    delete_edge(v2, v1);    
				}

			    if(!check_if_path(v2, v2, v1))
				{
				    add_child(v1, v2);
				    delta_for_21 = (*scoreFunc)(v1) + (*scoreFunc)(v2) - old_score;
				    delete_edge(v1, v2);
				}
			    
			    else  { continue; }  // Adding an edge either way causes a cycle so continue.
			    

			    if (delta_for_12 < max_delta && delta_for_21 < max_delta)  // No change is best.
				{
				    continue;
				}
			    else if (delta_for_12 > delta_for_21)  {  // v1 --> v2 is best choice.
				num_mods += 1; 
				option[0] = i; option[1] = j; option[2] = 1;
 				max_delta = delta_for_12;
			    }

			    else {  // v2 --> v1 is best choice.
				num_mods += 1;
				option[0] = j; option[1] = i; option[2] = 1;
				max_delta = delta_for_21;
			    }
			    
			}
		}
	}

    return num_mods;
}


int best_mod_if_connected(vertex *v1, vertex *v2, int *option, double (*scoreFunc)(vertex *), double *max_delta)
{
    // Takes two vertices, and tries to find a config that is better than current max_delta.
    // The local score is the change that increases the score the most.
    // If the local change increase the score than max_delta, then tentatively that change is adopted.
    // Otherwise, we continue onto next vertex pair.
    
    // Returns 0 if no better config was found.

    double delta_for_delete, delta_for_reverse;

    double old_score = (*scoreFunc)(v1) + (*scoreFunc)(v2);
	
    delete_edge(v1, v2);
    delta_for_delete = (*scoreFunc)(v2) + (*scoreFunc)(v1) - old_score;
    
    // If a path between v1 and v2 exists after deleting, then reversing the edge causes a  cycle.
    
    if (!check_if_path(v1 ,v1, v2)) // If there's no path, we can proceed to check score for reversed.
	{	    
	    add_child(v2, v1);
	    delta_for_reverse = (*scoreFunc)(v2) + (*scoreFunc)(v1) - old_score;
	    
	    if (delta_for_reverse < *max_delta && delta_for_delete < *max_delta)  // Doing nothing was best option.
		{
		    delete_edge(v2, v1); add_child(v1, v2);
		    return 0;
		}
	    
	    else if (delta_for_delete <= delta_for_reverse )  // Reversing the edge was best option.
		{
		    *max_delta = delta_for_reverse;
		    option[0] = v1->id; option[1] = v2->id; option[2] = -1;
		    
		    delete_edge(v2, v1); add_child(v1, v2);
		    return 1;
		}
	    
	    else // deleting the edge was best option.
		{
		    *max_delta = delta_for_delete;
		    option[0] = v1->id; option[1] = v2->id; option[2] = 0;
		    add_child(v1, v2);
		    return 1;
		}
	}
    
    
    else  // Reversing the edge isn't an option, so we just compare with delete edge.
	{
	    add_child(v1, v2);
		
		if (delta_for_delete > *max_delta)
		    {
			*max_delta = delta_for_delete;
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
