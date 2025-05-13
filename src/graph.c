#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "graph.h"
#include "linalg.h"
#include "score_functions.h"


void delete_edge(vertex *v, vertex *child)  // Deletes the edge assuming v --> child
{    
    
    int index = v->num_children;  // If child not in children then it just skips
    for (int i = 0; i < v->num_children; i++)
	{
	    if (v->children[i] == child) {index = i; break;}
	}

    // Ideally this shouldn't be getting used but it helps debugging.
    if (index == v->num_children) { printf(" %d NOT FOUND in %d.children, exiting \n", child->id, v->id); return; }

    else
	{
	for (; index < v->num_children; index++)
	    {
		v->children[index] = v->children[index+1];
	    }
	}
    v->num_children--;
    
    
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


int is_child(vertex *v, vertex *child)  // Returns 1 if there's an edge v --> child
{
    for (int i = 0; i < v->num_children; i++)
	{
	    if (v->children[i] == child)  { return 1;}
	}
    return 0;
}


void add_child(vertex *v, vertex *child)
{
    if (check_if_path(child, child, v))
	{
	    printf("Cannot add edge %d --> %d", v->id, child->id);
	    return;  // JUST DO NOTHING, IDEALLY, THIS SHOULD NOT BE GETTING CALLED ANYWAYS.
	}
    
    v->children[v->num_children++] = child;
    child->parents[child->num_parents++] = v;
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



