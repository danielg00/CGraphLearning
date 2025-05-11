void print_graph(DAG *G);
void free_graph(DAG *G);
void arraycpy(double **Copy, double **A, int *dims);
void printArray(double **Array, int dim0, int dim1);
void freeMatrix(matrix *M);
void alloc_array(matrix *A);
void freeArray(double ** array, int dim0);
matrix *matmul(matrix * A, matrix * B);
