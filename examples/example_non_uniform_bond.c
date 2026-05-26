#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "inc/graph_type.h"
#include "inc/graph_construction.h"
#include "inc/percolation_main.h"

/* bond-percolation simulation with non-uniform edge probabilities (no Newman-Ziff method) */
int main(){
  srand(time(NULL));

  /* define unit graph of simple cubic lattice */
  uint8_t nedge = 3; // number of different edges in unit graph representation
  uint8_t n_dim = 3; // number of dimensions

  uint8_t* blk_edge_label = malloc((size_t)nedge);
  for(uint8_t i=0; i<nedge; i++) blk_edge_label[i] = i;
  uint8_t* blk_edges = malloc((size_t)nedge * n_dim);
  memset(blk_edges, 0, (size_t)nedge * n_dim);
  int8_t* blk_vec = malloc((size_t)nedge * n_dim);
  memset(blk_vec, 0, (size_t)nedge * n_dim);
  blk_vec[0] = 1; blk_vec[4] = 1; blk_vec[8] = 1;

  UnitGraph u = new_unit_graph(blk_edges, blk_vec, blk_edge_label, nedge, n_dim);
  free(blk_edges);
  free(blk_vec);
  free(blk_edge_label);

  /* define the actual graph */
  int64_t lsize = 20; // size in each dimension
  bool periodic = false; // no periodic boundaries for start-stop percolation
  bool get_size = false; // do define start/stop set
  bool edge_list = false; // obsolete here since edge labels are specified
  Graph g = get_lattice_from_unit_graph(&u, lsize, n_dim, false, periodic, get_size, edge_list);

  /* run start-stop bond percolation */
  uint64_t num_rep = 100000; // repetitions for averaging
  float* pBond = malloc((size_t)nedge * sizeof(float));

  for (float p_bnd = 0.18; p_bnd < 0.30; p_bnd += 0.02){
    for(uint8_t i=0; i<nedge; i++) pBond[i] = p_bnd; // bond probabilities can be different for different edge labels
    uint64_t cnt = non_uniform_bond_percol(&g, pBond, num_rep);
    printf("%.15f\n", (double)cnt / num_rep);
  }

  free_unit_graph(&u);
  free_graph(&g);
  free(pBond);
  return 0;
}
