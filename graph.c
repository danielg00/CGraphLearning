#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "graph.h"


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

void add_child(vertex *v, vertex *child)  // Add checks,
{
    v->children[v->num_children] = child;
    v->num_children += 1;
    
    child->children[child->num_parents] = v;
    v->num_parents += 1;
}


DAG *init_graph(int num_nodes)
{
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
	}
    return G;
}


DAG *unicycle_graph()  // 1 -> 2 -> ... -> 4 -> 1 ...
{
    DAG *G = init_graph(5);
    add_child(&(G->nodes[4]), &(G->nodes[0]));
    for (int i = 1; i < 5; i++)
	{
	    add_child(&(G->nodes[i-1]), &(G->nodes[i]));
	}

    return G;
}



DAG *nocycle_graph()  // 1 -> 2 -> ... -> 4 ...
{
    DAG *G = init_graph(5);
    for (int i = 1; i < 5; i++)
	{
	    add_child(&(G->nodes[i-1]), &(G->nodes[i]));
	}

    return G;
}

DAG *disjoint_1cycle_graph()  // 1 -> 2 -> 3; 4 -> 5 -> 6 -> 4
{
    DAG *G = init_graph(6);
    
    add_child(&(G->nodes[0]), &(G->nodes[1]));
    add_child(&(G->nodes[1]), &(G->nodes[2]));
    add_child(&(G->nodes[0]), &(G->nodes[2]));
    
    add_child(&(G->nodes[3]), &(G->nodes[4]));
    add_child(&(G->nodes[4]), &(G->nodes[5]));
    add_child(&(G->nodes[3]), &(G->nodes[5]));

    return G;
}


int graph_test()
{
    DAG *G = disjoint_1cycle_graph();
    int C = check_cyclic(G);
    printf("\n cycle: %d", C);

    free(G);
    return 0;
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
