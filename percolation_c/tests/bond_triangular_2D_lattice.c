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
  srand(23635);
  int64_t lsize = 1000; // size in each dimension
  bool edge_list = true; // keep track of edges in list (needed for bond percolation)
  int r;
  Graph g = get_2d_triangular(lsize, false, true, edge_list);

  float sweep_prob[SWEEP_LEN] = {0.30, 0.40}; // note: expected bond-percolation threshold at 0.347 in this example
  int64_t idxLambda;
  int64_t* largest_component = percolate_bond(&g, 1.0, &idxLambda); // specify (fixed) site probability
  float* expectation_value = get_expectation_value(sweep_prob, SWEEP_LEN, largest_component, g.num_edges, 0, g.nnode, true, true);  // compute a bond-percolation curve (from bond sweep)
  r = expectation_value[0] < 0.1 ? 0 : -1;
  r |= expectation_value[1] > 0.5 ? 0 : -1;
  free(expectation_value);
  free(largest_component);

  free_graph(&g);
  return r;
}

