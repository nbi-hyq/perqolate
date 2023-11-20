#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include "../inc/graph_type.h"

/* new_graph: number of nodes/qubits
   num_nb_max: maximum number of neighbors for a node
   get_size: check for size or cluster spanning
   edge_list: keep track of edges e.g. for bond percolation
   fusion_list: keep track of qubits used for fusions (true only needed for reference code that does not use any Newman-Ziff algorithm) */
Graph new_graph(int64_t nnode, uint8_t num_nb_max, bool get_size, bool edge_list, bool fusion_list){
  Graph g;
  g.ptr = malloc(nnode * sizeof(int64_t));  // several meanings: (if ptr[i]>0: parent index ("pointer"), elif ptr[i]==LOST: detached or lost node, ptr[i]==APPLYZ: to be measured in Z ,elif ptr[i]<0: - size of component)
  g.nn = malloc(nnode * (size_t)num_nb_max * sizeof(int64_t)); // neighbors of a node
  g.len_nb = malloc(nnode); // until which index there are neighbors (255 neighbors max)
  g.edge_node = edge_list ? malloc((nnode * (size_t)num_nb_max / 2 + 1) * sizeof(int64_t)) : NULL; // 1st node in edge (only for looping over edges)
  g.edge_idx = edge_list ? malloc(nnode * (size_t)num_nb_max / 2 + 1) : NULL; // index of 2nd node of edge in nn (only for looping over edges)
  g.fusion_partner = fusion_list ? malloc(nnode * sizeof(int64_t)) : NULL; // partner node for fusion, -1 means no fusion node
  g.fusion_success = fusion_list ? malloc(nnode * sizeof(bool)) : NULL; // successful fusion or not
  g.fusion_node = fusion_list ? malloc(nnode * sizeof(bool)) : NULL; // node that is used for fusion
  g.static_node = malloc(nnode * sizeof(bool)); // node that cannot be lost (or removed in classical case)
  g.start = get_size ? NULL : malloc(nnode * sizeof(bool)); // nodes where percolation starts
  g.stop = get_size ? NULL : malloc(nnode * sizeof(bool)); // nodes where percolation stops
  g.num_edges = 0; // count number of edges
  g.num_nb_max = num_nb_max; // maximum number of neighbors per node
  g.nnode = nnode; // number of nodes or central qubits
  g.get_size = get_size; // true: get largest connected component size, false: cluster spanning
  g.failed = 0; // flag that can be set to indicate that something is wrong with the graph construction
  return g;
}

void free_graph(Graph* g){
  free(g->ptr);
  free(g->nn);
  free(g->edge_node);
  free(g->edge_idx);
  free(g->len_nb);
  free(g->fusion_node);
  free(g->fusion_partner);
  free(g->fusion_success);
  free(g->static_node);
  free(g->start);
  free(g->stop);
}
