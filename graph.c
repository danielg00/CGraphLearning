#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "graph.h"
#include "score_functions.h"


// 2 Philosophies:
// (1) A graph is a list of vertices, each vertex contains properties in relation to the others in the list.
// So graph.vertex.parents makes sense.
// Graph is a single object and easy to manage.

// (2) A graph is a structure we impose on a list of vertices.
// So, vertex.graph.parents makes sense. We can easily union and manipulate graphs with this since vertices are ultimately just
// locations in memory and not able to be copied, only pointed to.

DAG *unicycle_graph();
DAG *nocycle_graph();
DAG *disjoint_1cycle_graph();
int graph_test();


void delete_edge(vertex *v, vertex *child)
{
    /* printf("DELETING %d --> %d\n", v->id, child->id); */
    
    
    int index = v->num_children;  // If child not in children then it just skips
    for (int i = 0; i < v->num_children; i++)
	{
	    if (v->children[i] == child) {index = i; break;}
	}
    
    if (index == v->num_children) { printf(" %d NOT FOUND in %d.children, exiting \n", child->id, v->id); return; }

    else{
	
	for (; index < v->num_children; index++)
	    {
		v->children[index] = v->children[index+1];
	    }
    }
    v->num_children--;

    if (v->num_children < 0) {fprintf(stderr, "Num children less than zero"); exit(1);}
    
    index = child->num_parents; 
    for (int i = 0; i < child->num_parents; i++)
	{
	    if (child->parents[i] == v) {index = i; break;}
	}

    
    for (; index < child->num_parents; index++)
	{
	    child->parents[index] = child->parents[index+1];
	}
	
    child->num_parents--;
}



int is_child(vertex *v, vertex *child)  // Returns 1 if has edge v --> child
{
    for (int i = 0; i < v->num_children; i++)
	{
	    if (v->children[i] == child)  { return 1;}
	}
    return 0;
}


int add_child(vertex *v, vertex *child)  // Add checks,
{
    /* printf("CALLED ADD_CHILD WITH %d --> %d\n", v->id, child->id); */

    if (check_if_path(child, child, v))
	{
	    /* printf("ADD_CHILD WITH %d --> %d CREATES CYCLE, IGNORING \n", v->id, child->id); */
	    
	    return 0;  // JUST DO NOTHING AND RETURN 0
	}
    v->children[v->num_children++] = child;
    child->parents[child->num_parents++] = v;
    return 1;
}


int check_if_path(vertex *start, vertex *current, vertex *end)  // Returns 1 if there is a path from v1 to v2
{
    if (end == current)
	{
	    return 1;
	}
    
    else
	{
	    for (int i = 0; i < current->num_children; i++)
		{
		    if (check_if_path(start, current->children[i], end))
			{
			    return 1;
			}
		}
	}
    
    return 0;
}
	    

void free_graph(DAG *G)
{
    int D = G->num_nodes;
    for (int i = 0; i < D; i++)
	{
	    free(G->nodes[i].children);
	    free(G->nodes[i].parents);

	}
    free(G->nodes);
    free(G);
}


void print_graph(DAG *G)
{
    for (int i = 0; i < G->num_nodes; i++)
	{
	    printf(" ( %d ) ====> (", i);
	    for (int j = 0; j < G->nodes[i].num_children; j++)
		{
		    printf(" %d ", G->nodes[i].children[j]->id);
		}
	    printf(")\n"); 
	}
}
