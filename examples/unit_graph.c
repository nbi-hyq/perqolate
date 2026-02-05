#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"
#include "../inc/percolation_main.h"

/* construct periodic graph from unit-graph */
int main(){
  srand(6735829);
  int64_t lsize = 200; // size in each dimension
  int n_avg = 1000; // number of repetitions for averaging
  bool periodic = false; // no periodic boundaries
  bool get_size = false; // check percolation, not determine size of largest connected component
  bool static_center = false;
  uint8_t dimension = 2;
  bool list_edges = true; // must be true for bond-percolation

  /* create the unit graph (for honey-comb lattice) */
  uint8_t nedge = 3;
  uint8_t* blk_edges = malloc((size_t)nedge * dimension);
  int8_t* blk_vec = malloc((size_t)nedge * dimension);
  blk_edges[0] = 0; blk_edges[1] = 1; blk_edges[2] = 0; blk_edges[3] = 1; blk_edges[4] = 0; blk_edges[5] = 1;
  blk_vec[0] = 0; blk_vec[1] = 0; blk_vec[2] = -1; blk_vec[3] = 0; blk_vec[4] = 0; blk_vec[5] = -1;
  UnitGraph unt = new_unit_graph(blk_edges, blk_vec, nedge, dimension);
  free(blk_edges);
  free(blk_vec);

  /* create periodic graph */
  Graph g = get_lattice_from_unit_graph(&unt, lsize, dimension, static_center, periodic, get_size, list_edges);

  /* run simulation */
  int64_t idxLambda;
  double avg = 0;
  for(int i=0; i < n_avg; i++) {
    int64_t* percolated = percolate_site(&g, 1.0, &idxLambda); // site-percolation given fixed bond probability (Newman-Ziff method)
    free(percolated);
    avg += (double)idxLambda/g.nnode;
  }
  printf("%f\n", avg/n_avg);

  avg = 0;
  for(int i=0; i < n_avg; i++) {
    int64_t* percolated = percolate_bond(&g, 1.0, &idxLambda); // bond-percolation given fixed site probability (Newman-Ziff method)
    free(percolated);
    avg += (double)idxLambda/g.num_edges;
  }
  printf("%f\n", avg/n_avg);

  free_graph(&g);
  free_unit_graph(&unt);
  return 0;
}
