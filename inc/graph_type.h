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

typedef struct {
  uint8_t** nb; // adjacency list of node neighbors
  uint8_t** i_vec; // same structure as nb, storing connection vector position of edge in blk_vec
  uint8_t* num_nb; // number of neighbors stored in adjacency list nb (<= valency in p-graph)
  int8_t* blk_vec; // connection vectors in block representation
  uint8_t nnode, dim, num_nb_max; // number of nodes, dimension of p-graph, maximum number of neighbors of a node in p-graph
} UnitGraph;

Graph new_graph(int64_t nnode, uint8_t num_nb_max, bool get_size, bool edge_list, bool fusion_list);
UnitGraph new_unit_graph(uint8_t* blk_edges, int8_t* blk_vec, uint8_t nedge, uint8_t dim);
void free_graph(Graph* g);
void free_unit_graph(UnitGraph* u);

#endif
