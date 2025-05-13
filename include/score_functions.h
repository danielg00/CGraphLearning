double BIC_score(vertex *v);
double graph_score(DAG *G, double (*ScoreFunc)(vertex*));

int find_best_mod(DAG *G, int *option, double (*scoreFunc)(vertex *));
int best_mod_if_connected(vertex *v1, vertex *v2, int *option, double (*scoreFunc)(vertex *), double *max_delta);
