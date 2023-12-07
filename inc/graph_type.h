#ifndef GRAPH_TYPE
#define GRAPH_TYPE

typedef struct {
  int64_t* ptr;
  int64_t* nn;
  uint8_t* len_nb;
  int64_t* edge_node;
  uint8_t* edge_idx;
  int64_t* fusion_partner;
  bool* fusion_success;
  bool* fusion_node;
  bool* static_node;
  bool* start;
  bool* stop;
  int64_t nnode, num_edges;
  uint8_t num_nb_max;
  bool get_size;
  bool failed;
} Graph;

Graph new_graph(int64_t nnode, uint8_t num_nb_max, bool get_size, bool edge_list, bool fusion_list);
void free_graph(Graph* g);

#endif
