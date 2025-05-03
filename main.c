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
    double ** array = malloc(dim[0]*sizeof(double *));

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
    
    matrix * A = malloc(sizeof(*A));
    A->dims = get_array_size(f);
    A->data = load_array_from_file(f, A->dims);
       
    matrix * inv;
    inv = invert_matrix(A);
    
    /* printf("FINAL: %f", inv->data[0][0]); */
	
    return 0;
}




