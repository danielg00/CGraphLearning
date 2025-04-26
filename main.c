#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "linalg.h"

int * get_array_size(FILE * f) {
    int * dimension = malloc(2 * sizeof(int));
    char buffer[7];
    /* buffer[7] = '\0'; */
    fread(buffer, 1, 6, f);
    char file_name[] = {0x93, 'N', 'l', 'M', 'P', 'Y', '\0'};  // FIX
    assert(memcmp(buffer, file_name, 6) && "dasdas");
	
    char dummy[2]; //FILE VERSIONS - DONT CARE
    fread(&dummy[0], 1, 1, f);
    fread(&dummy[1], 1, 1, f);

    unsigned short int len; // HEADER LENGTH
    fread(&len, 2, 1, f);

    
    char buffer2[200];
    for (int i = 0; i < 100; i++) 
	{
	fscanf(f, "%s", buffer2);
	if (strncmp(buffer2, "'shape':", 8) == 0)
	    {
		break;
	    }
	}

    fgetc(f); /* HACK */
    fscanf(f, "(%d,%d)", &dimension[0], &dimension[1]);

    
    int offset = (floor(len/64) + 1) * 64; // MOVE POINTER TO START OF FLOAT VALUES
    fseek(f, offset, SEEK_SET);
    return dimension;
}

double ** load_array_from_file(FILE * f, int * dim) {
    double ** array = malloc(dim[0]*sizeof(float *));

    for (int i = 0; i < dim[0]; i++)
	{
	    array[i] = malloc(dim[1]*sizeof(double));
	    fread(array[i], sizeof(double), dim[1], f);
	}
    
    assert(fgetc(f) == EOF && "EOF NOT REACHED AFTER READING MATRIX VALUES");
    return array;
}

int main() {
    FILE * f = fopen("test_data/test_d.npy", "rb");
    
    int * dim = get_array_size(f);
    
    printf("\n dim -> (%d %d) \n", dim[0], dim[1]);
    
    double ** array = load_array_from_file(f, dim);
    matrix * A = malloc(sizeof(*A));
    A->dims = malloc(2*sizeof(int));
    A->dims[0] = 6; A->dims[1] = 6;
    A->data = array;
    printf(" asd %f sdddsa", A->data[0][0]);
    double *** ret_v;
    ret_v = LU_decomposition(A);
    printf("\n %f", ret_v[0][3][2]);
    printf("\n %f", ret_v[1][2][3]);
	
    return 0;
}




