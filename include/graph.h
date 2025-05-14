typedef struct vertex
{
    int id;
    int num_children;
    int num_parents;
    int num_samples;
    
    struct vertex **children;
    struct vertex **parents;
    double * data;  // todo: change to `samples`

} vertex;


typedef struct {
    struct vertex *nodes;
    int num_nodes;
} DAG;


void add_child(vertex *v, vertex *child);
int is_child(vertex *v1, vertex *child);
void delete_edge(vertex *v, vertex *child);


int check_if_path(vertex *start, vertex *current, vertex *end);


void free_graph(DAG *G);
void print_graph(DAG *G);
