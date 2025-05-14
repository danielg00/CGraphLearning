struct
{
    char * FNAME;
    char * OUT_FNAME;
    int MAX_ITERS;
    int samples_on_rows;
    int verbose; 
} arguments ; 

matrix *load_matrix_from_file(char * fname);
    
