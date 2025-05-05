#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "linalg.h"
#include "io.h"
#include "graph.h"


int main() {
    char * fname = "test_data/test_d.npy";
    
    matrix *A = load_matrix_from_file(fname);

    graph_test();

    
    freeMatrix(A);	
    return 0;
}
