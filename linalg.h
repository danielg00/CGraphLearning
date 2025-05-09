typedef struct {
    int * dims;
    double ** data;    
} matrix;


double ** matmul(matrix * A, matrix * B);
int * LU_decomposition(matrix * A, double *** L, double *** U);
matrix * invert_matrix(matrix * A);
void solve_Ax_b(double **A, double **b, int D, int upper);

void printArray(double **Array, int dim0, int dim1);
void freeArray(double ** array, int N);
void freeMatrix(matrix *M);
