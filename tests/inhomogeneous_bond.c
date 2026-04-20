#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"
#include "../inc/percolation_main.h"

/* bond-percolation simulation with inhomogeneous edge probabilities (no Newman-Ziff method) */
int main(){
  srand(456789);

  /* define unit graph of triangular lattice */
  uint8_t nedge = 3; // number of different edges in unit graph representation
  uint8_t n_dim = 2; // number of dimensions

  uint8_t* blk_edge_label = malloc((size_t)nedge);
  for(uint8_t i=0; i<nedge; i++) blk_edge_label[i] = i;
  uint8_t* blk_edges = malloc((size_t)nedge * n_dim);
  memset(blk_edges, 0, (size_t)nedge * n_dim);
  int8_t* blk_vec = malloc((size_t)nedge * n_dim);
  memset(blk_vec, 0, (size_t)nedge * n_dim);
  blk_vec[0] = 1; blk_vec[3] = 1; blk_vec[4] = 1; blk_vec[5] = -1;

  UnitGraph u = new_unit_graph(blk_edges, blk_vec, blk_edge_label, nedge, n_dim);
  free(blk_edges);
  free(blk_vec);
  free(blk_edge_label);

  /* define the actual graph */
  int64_t lsize = 50; // size in each dimension
  bool periodic = false; // no periodic boundaries for start-stop percolation
  bool get_size = false; // do define start/stop set
  bool edge_list = false; // obsolete here since edge labels are specified
  Graph g = get_lattice_from_unit_graph(&u, lsize, n_dim, false, periodic, get_size, edge_list);

  /* run start-stop bond percolation */
  uint64_t num_rep = 1000; // repetitions for averaging
  float* pBond = malloc((size_t)nedge * sizeof(float));
  for(uint8_t i=0; i<nedge; i++) pBond[i] = 0.3; // some numbers such that by adapting pBond[0] one can reach criticality
  pBond[0] = 0.4; // start close to criticality for test speed

  int num_step = 10;
  for (int i=0; i<num_step; i++){
    uint64_t cnt = non_uniform_bond_percol(&g, pBond, num_rep);
    double p_ratio = (double)cnt / num_rep;
    pBond[0] += (0.5 - p_ratio) / pow(lsize, 0.75); // use slope at criticality and fact that percolation probabilility if 0.5 in 2d at criticality
    printf("%.15f %f %f\n", p_ratio, pBond[0], 1 - pBond[0] - pBond[1] - pBond[2] + pBond[0]*pBond[1]*pBond[2]);
  }

  int r = (1 - pBond[0] - pBond[1] - pBond[2] + pBond[0]*pBond[1]*pBond[2] > 0.015); // check criterion from https://doi.org/10.1063/1.1704215
  r |= (1 - pBond[0] - pBond[1] - pBond[2] + pBond[0]*pBond[1]*pBond[2] < -0.015); // check criterion from https://doi.org/10.1063/1.1704215
  free_unit_graph(&u);
  free_graph(&g);
  free(pBond);
  return r;
}
