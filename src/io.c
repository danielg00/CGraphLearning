#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "linalg.h"


matrix *load_matrix_from_file(char *fname) {
    FILE *f = fopen(fname, "rb");
    
    matrix *A = malloc(sizeof(matrix));
    A->dims = malloc(2 * sizeof(int));

    
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
    fscanf(f, "(%d,%d)", &(A->dims[0]), &(A->dims[1]));

    
    int offset = (floor(len/64) + 1) * 64; // MOVE POINTER TO START OF FLOAT VALUES
    fseek(f, offset, SEEK_SET);

    A->data = malloc(A->dims[0]*sizeof(double *));

    for (int i = 0; i < A->dims[0]; i++)
	{
	    A->data[i] = malloc(A->dims[1]*sizeof(double));
	    fread(A->data[i], sizeof(double), A->dims[1], f);
	}

    assert(fgetc(f) == EOF && "EOF NOT REACHED AFTER READING MATRIX VALUES");
    
    fclose(f);
    return A;
}
