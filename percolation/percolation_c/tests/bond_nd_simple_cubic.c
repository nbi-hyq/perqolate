#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"
#include "../inc/percolation_main.h"

#define SWEEP_LEN 2

int main(){
  srand(2347624);
  int64_t lsize = 100; // size in each dimension
  uint8_t n_dim = 3; // dimension
  bool periodic = false;
  bool static_center = false;
  bool get_size = false;
  bool edge_list = true; // keep track of edges in list (needed for bond percolation)
  int r; // return
  Graph g = get_nd_simple_cubic(lsize, n_dim, static_center, periodic, get_size, edge_list);

  float sweep_prob[SWEEP_LEN] = {0.243, 0.255}; // note: expected bond-percolation threshold at 0.2488
  int64_t idxLambda;
  int64_t* largest_component = percolate_bond(&g, 1.0, &idxLambda); // specify (fixed) site probability
  float* expectation_value = get_expectation_value(sweep_prob, SWEEP_LEN, largest_component, g.num_edges, 0, 1, true, false);  // compute a bond-percolation curve (from bond sweep)
  r = expectation_value[0] < 0.1 ? 0 : -1;
  r |= expectation_value[1] > 0.5 ? 0 : -1;
  printf("%f\n", (double)idxLambda/g.num_edges);
  free(expectation_value);
  free(largest_component);
  free_graph(&g);
  return r;
}

