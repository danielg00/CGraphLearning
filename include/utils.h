void arraycpy(double **Copy, double **A, int *dims);
void printArray(double **Array, int dim0, int dim1);

void alloc_array(matrix *A);
void freeMatrix(matrix *M);
void freeArray(double ** array, int dim0);

matrix *matmul(matrix * A, matrix * B);

DAG *init_graph(matrix *data);
void print_graph(DAG *G);
void free_graph(DAG *G);
