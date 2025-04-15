#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/* #include "graph.h" */


float index_val(int i, int j, float ** array) {
    return array[i][j];
}

float dot_product(void) {
    return 0;
}

int * get_array_size(FILE * f) {
    int * dimension = malloc(2 * sizeof(int));
    char buffer[6];
    fread(buffer, 1, 6, f);
    char file_name[] = {0x93, 'N', 'U', 'M', 'P', 'Y'};
    printf("%s \n", file_name);
    printf("%s", buffer);
    /* assert(memcmp(buffer, file_name, sizeof(buffer)) && "dasdas"); */
	
    char dummy[2]; /* FILE VERSIONS - DONT CARE */
    fread(&dummy[0], 1, 1, f);
    fread(&dummy[1], 1, 1, f);

    unsigned short int len; /* HEADER LENGTH */
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

    
    int offset = (floor(len/64) + 1) * 64; /* MOVE POINTER TO START OF FLOAT VALUES */
    fseek(f, offset, SEEK_SET);
    return dimension;
}

float ** load_array_from_file(FILE * f, int * dim) {
    float ** array = malloc(dim[0]*sizeof(float *));

    for (int i = 0; i < dim[0]; i++)
	{
	    array[i] = malloc(dim[1]*sizeof(float));
	    fread(array[i], sizeof(float), dim[1], f);
	}
    
    assert(fgetc(f) == EOF && "EOF NOT REACHED AFTER READING MATRIX VALUES");
    return array;
}

int main() {
    FILE * f = fopen("test_d.npy", "rb");
    
    int * dim = get_array_size(f);
    
    printf("\n dim -> (%d %d) \n", dim[0], dim[1]);
    
    float ** array = load_array_from_file(f, dim);
    for (int i = 0; i < dim[0]; i++)
	{
	    for (int j = 0; j < dim[1]; j++)
		{
		    printf("%f ", array[i][j]);
		}
	    printf("\n");
	}
    
    
    return 0;
}




