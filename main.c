#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "linalg.h"
#include "io.h"
#include "graph.h"
#include "score_functions.h"

#define ITERS 10
DAG *init_graph(matrix *A);


int main()
{
    char * fname = "test_data/test_d.npy";
    
    matrix *A = load_matrix_from_file(fname);  // features x samples

    DAG *G = init_graph(A);
    printf("%d", G->num_nodes);
    
    // HILL-CLIMB
    int iter = 0; double current_score;
    vertex *v1, *v2;
    while(iter < ITERS)
	{
	    printf("============== ITERATION: %d =============", iter);
	    for (int i = 0; i < G->num_nodes; i++)
		{
		    for (int j = i+1; j < G->num_nodes-1; j++)
			{
			    printf("\n %d, %d \n", i, j);
			    v1 = &(G->nodes[i]); v2 = &(G->nodes[j]);
			    double score_v1 = BIC_score(v1); double score_v2 = BIC_score(v2);
			    current_score = score_v1 + score_v2;
			    double score_delete, score_reverse;
			    if (is_child(v1, v2))  // i.e. v2 --> v1
				{
				    printf("GERE  ");
				    delete_edge(v1, v2);
				    score_v2 = BIC_score(v2);
				    score_delete = score_v1 + score_v2;  // only changes to parents affect score.

				    // If there exists a path from v1 to v2 after deleting the edge,
				    // then reversing the edge causes a cycle, so we skip over that. 
				    
				    add_child(v2, v1);  // reverse edge
				    score_reverse = BIC_score(v1) + score_v2;
					
					


				    if (score_delete > current_score  && score_delete > score_reverse) { delete_edge(v1, v2); continue; }
				    else if (score_reverse > current_score && !check_if_path(v1, v1, v2)) { continue; }
				    else { delete_edge(v2, v1); add_child(v1, v2); continue; }
				}
				    

			    if (is_child(v2, v1)) // i.e. v1 --> v2 we do the same as above but reverse variable name.
				{
				    delete_edge(v2, v1);
				    score_v1 = BIC_score(v1);
				    score_delete = score_v2 + score_v1;  // only changes to parents affect score.
				    
				    if (check_if_path(v2, v2, v1))
					{
					    score_reverse = current_score;
					    add_child(v2, v1);
					}
				    
				    else
					{
					    add_child(v1, v2); 
					    score_reverse = BIC_score(v2) + score_v1;
					}


				    if (score_reverse >= current_score)  { continue; }
				    else if (score_delete >= current_score) { delete_edge(v2, v1); continue; }
				    else { delete_edge(v1, v2); add_child(v2, v1); continue; }
				}

			    else  // No edge between v1 and v2
				{
				    double score1, score2;
				    add_child(v1, v2);
				    score1 = BIC_score(v1) + BIC_score(v2);

				    delete_edge(v1, v2);
				    add_child(v2, v1);
				    score2 = BIC_score(v1) + BIC_score(v2);
				    if (score1 > score2)
					{
					    delete_edge(v2, v1);
					    add_child(v1, v2);
					}
				}

			}
		}
	    
	    iter++;
	}

    printf(" %d ", G->nodes[0].num_children);
    freeMatrix(A);	
    return 0;
}


int main2()
{
    char * fname = "test_data/test_d.npy";
    
    matrix *A = load_matrix_from_file(fname);  // features x samples

    DAG *G = init_graph(A);

    add_child(&G->nodes[0], &G->nodes[1]);
    add_child(&G->nodes[1], &G->nodes[2]);
    add_child(&G->nodes[2], &G->nodes[3]);

    int C = check_if_path(&G->nodes[0], &G->nodes[0], &G->nodes[0]);
    printf("\n asdasdasdasdasdasd %d \n", C);
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
		
	}
    return G;
}
