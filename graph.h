typedef struct vertex
{
    int id;
    int num_children;
    int num_parents;
    struct vertex **children;
    struct vertex **parents;
    double * data;
    int num_samples;
    
    int calc_score_needed; // score needs to calculated if 1.
    double score;
} vertex;


typedef struct {
    struct vertex *nodes;
    int num_nodes;
} DAG;

int add_child(vertex *v, vertex *child);
int is_child(vertex *v1, vertex *child);
void delete_edge(vertex *v, vertex *child);

int check_cyclic(DAG *G);
int check_node_in_cycle(vertex current_node, int *checked_vertices);
int check_if_path(vertex *start, vertex *current, vertex *end);

void free_graph(DAG *G);
void print_graph(DAG *G);
int graph_test();
