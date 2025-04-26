/* float ** matmul(float ** A, float **B); */
typedef struct Matrix {
    int * dims;
    double ** data;    
} matrix;

double *** LU_decomposition(matrix * A);
