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
    printf("DELETING %d --> %d\n", v->id, child->id);

    int index = v->num_children;  // If child not in children then it just skips
    for (int i = 0; i < v->num_children; i++)
	{
	    if (v->children[i] == child) {index = i; break;}
	}

    
    for (; index < v->num_children; index++)
	{
	    v->children[index] = v->children[index+1];
	}
	
    v->num_children--;
    if (v->num_children < 0) {printf("!!!!!CANT DELETE %d --> %d !!!!!!\n", v->id, child->id);}
    
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


void add_child(vertex *v, vertex *child)  // Add checks,
{
    printf("CALLED ADD_CHILD WITH %d --> %d\n", v->id, child->id);
    v->children[v->num_children++] = child;
    child->parents[child->num_parents++] = v;
}


int check_list(vertex current_node, int *checked_vertices) // returns 1 if current node is in node list.
{
    if (checked_vertices[current_node.id] == 0)
	{
	    return 0;
	}
    else
	{
	    return 1;
	}
}


int check_node_in_cycle(vertex current_node, int *checked_vertices) // Want to try and do everything here on the stack. 
{
    int C = check_list(current_node, checked_vertices); // Checks if current node in check_list
    int C2;

    switch(C)
	{
	case 1:
	    
	    return 1; // Current node is in the checked_list -> its a cycle.

	case 0:

	    checked_vertices[current_node.id] += 1; // Add current node to checked vertices.
	    
	    for (int i = 0; i < current_node.num_children; i++)
		{
		    C2 = check_node_in_cycle(*(current_node.children[i]), checked_vertices);
		    
		    if (C2 == 1) { return 1; }
		}
	}

    return 0;  // No cycles found if it reaches this part.
}

    
int check_cyclic(DAG *G)  // This is slighly broken; it assumes every node has only one child. need at another loop over children.
{
    // We start at node 0 = current_node;
    // If current_node is checked, we continue onto next iter
    // else if current node is unchecked, we check it, if it is part of a cycle
    // we return 1
    // else if its not apart of a cycle we continue.

    int N = G->num_nodes;
    int checked_vertices[N];
    memset(checked_vertices, 0, N*sizeof(int));
    
    
    int C;  // 1 if node apart of cycle
    for (int i = 0; i < N; i++)
	{
	    if (checked_vertices[i] == 0)
		{
		    C = check_node_in_cycle(G->nodes[i], checked_vertices);
		    if (C == 1) { return 1; }
		}
	    
	    else
		{
		    continue;
		}
	}
    return C;  // Returns 1 if G is cyclic.
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


// ==================== TOY GRAPH FUNCTIONS =================

/* DAG *unicycle_graph()  // 1 -> 2 -> ... -> 4 -> 1 ... */
/* { */
/*     DAG *G = init_graph(5); */
/*     add_child(&(G->nodes[4]), &(G->nodes[0])); */
/*     for (int i = 1; i < 5; i++) */
/* 	{ */
/* 	    add_child(&(G->nodes[i-1]), &(G->nodes[i])); */
/* 	} */

/*     return G; */
/* } */



/* DAG *nocycle_graph()  // 1 -> 2 -> ... -> 4 ... */
/* { */
/*     DAG *G = init_graph(5); */
/*     for (int i = 1; i < 5; i++) */
/* 	{ */
/* 	    add_child(&(G->nodes[i-1]), &(G->nodes[i])); */
/* 	} */

/*     return G; */
/* } */

/* DAG *disjoint_1cycle_graph()  // 1 -> 2 -> 3; 4 -> 5 -> 6 -> 4 */
/* { */
/*     DAG *G = init_graph(6); */
    
/*     add_child(&(G->nodes[0]), &(G->nodes[1])); */
/*     add_child(&(G->nodes[1]), &(G->nodes[2])); */
/*     add_child(&(G->nodes[0]), &(G->nodes[2])); */
    
/*     add_child(&(G->nodes[3]), &(G->nodes[4])); */
/*     add_child(&(G->nodes[4]), &(G->nodes[5])); */
/*     add_child(&(G->nodes[3]), &(G->nodes[5])); */

/*     return G; */
/* } */


/* int graph_test() */
/* { */
/*     DAG *G = disjoint_1cycle_graph(); */
/*     int C = check_cyclic(G); */
/*     printf("\n cycle: %d", C); */

/*     free(G); */
/*     return 0; */
/* } */
