#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
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

/* create new unit graph that defines a periodic graph
input: unidirectional unit graph in block-representation with node numbering 0..N-1
output: UnitGraph struct (undirected unit graph in adjacency list representation)
blk_edges: one edge is two node numbers next to each other
blk_vec: one block of dim numbers is a connection vector
nedge: number of edges in unit graph
dim: dimension of p-graph */
UnitGraph new_unit_graph(uint8_t* blk_edges, int8_t* blk_vec, uint8_t nedge, uint8_t dim){
  UnitGraph u;
  u.dim = dim;
  u.nnode = 0;
  for(uint8_t i=0; i<nedge; i++){
    if(blk_edges[2*i]+1 > u.nnode) u.nnode = blk_edges[2*i]+1; // assume that numbering is node 0..N-1
    if(blk_edges[2*i+1]+1 > u.nnode) u.nnode = blk_edges[2*i+1]+1; // assume that numbering is node 0..N-1
  }
  u.num_nb = malloc(sizeof(uint8_t) * u.nnode);
  memset(u.num_nb, 0, sizeof(uint8_t) * u.nnode);
  for(uint8_t i=0; i<nedge; i++){
    u.num_nb[blk_edges[2*i]]++;
    u.num_nb[blk_edges[2*i+1]]++;
  }

  u.nb = malloc(sizeof(uint8_t*) * u.nnode);
  u.i_vec = malloc(sizeof(uint8_t*) * u.nnode);
  for(uint8_t i=0; i<u.nnode; i++){
    u.nb[i] = malloc(sizeof(uint8_t) * u.num_nb[i]);
    u.i_vec[i] = malloc(sizeof(uint8_t) * u.num_nb[i]);
  }

  /* blk_vec: unidirection representation, u.blk_vec: two-directional representation of unit-graph */
  u.blk_vec = malloc(2 * sizeof(uint8_t) * nedge * dim);
  int blk_pos = 0;
  uint8_t* nb_pos = malloc(sizeof(uint8_t) * u.nnode);
  memset(nb_pos, 0, sizeof(uint8_t) * u.nnode);
  for(uint8_t i=0; i<nedge; i++){
    uint8_t nd0 = blk_edges[2*i];
    uint8_t nd1 = blk_edges[2*i+1];
    u.nb[nd0][nb_pos[nd0]] = nd1;
    u.i_vec[nd0][nb_pos[nd0]++] = blk_pos;
    for(uint8_t d=0; d<dim; d++) u.blk_vec[blk_pos * dim + d] = blk_vec[i * dim + d];
    blk_pos++;
    u.nb[nd1][nb_pos[nd1]] = nd0;
    u.i_vec[nd1][nb_pos[nd1]++] = blk_pos;
    for(uint8_t d=0; d<dim; d++) u.blk_vec[blk_pos * dim + d] = - blk_vec[i * dim + d]; // also store vector in opposite direction
    blk_pos++;
  }

  u.num_nb_max = 0;
  for(uint8_t i=0; i<u.nnode; i++){
    if(u.num_nb[i] > u.num_nb_max) u.num_nb_max = u.num_nb[i];
  }
  free(nb_pos);

  return u;
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

void free_unit_graph(UnitGraph* u){
  for(uint8_t i=0; i<u->nnode; i++){
    free(u->nb[i]);
    free(u->i_vec[i]);
  }
  free(u->nb);
  free(u->i_vec);
  free(u->num_nb);
  free(u->blk_vec);
}
