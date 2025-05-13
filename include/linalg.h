typedef struct {
    int * dims;
    double ** data;    
} matrix;


matrix *matmul(matrix * A, matrix * B);
int * LU_decomposition(matrix * A, double *** L, double *** U);
matrix * invert_matrix(matrix * A);
void solve_Ax_b(double **A, double **b, int D, int upper);

double variance_of_residuals(matrix *X, double *Y);

void printArray(double **Array, int dim0, int dim1);
void alloc_array(matrix *A);

void freeArray(double ** array, int N);
void freeMatrix(matrix *M);
void arraycpy(double **A, double **Copy, int *dims);
