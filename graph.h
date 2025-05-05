typedef struct vertex
{
    int id;
    int num_children;
    int num_parents;
    struct vertex **children;
    struct vertex **parents;
} vertex;


typedef struct {
    struct vertex *nodes;
    int num_nodes;
} DAG;

DAG *init_graph(int num_nodes);
void add_child(vertex *v, vertex *child);

int check_cyclic(DAG *G);
int check_node_in_cycle(vertex current_node, int *checked_vertices);

void free_graph(DAG *G);
int graph_test();
