#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"
#include "../inc/percolation_main.h"

#define SWEEP_LEN 2

/* test bond-percolation with get_nd_simple_cubic_block (using Newman-Ziff method) */
int main(){
  srand(25798723);
  int r = 0; // return
  int64_t idxLambda;

  int lblock = 1;
  uint8_t n_dim = 2;
  bool static_center = false;
  int64_t lsize = 100; // size in each dimension
  bool periodic = true; // periodic boundaries
  bool get_size = true; // determine size of largest connected component
  bool edge_list = true; // keep track of edges in list (needed for bond percolation)
  float sweep_prob[SWEEP_LEN] = {0.66, 1.0}; // bond probability (tuned)

  Graph g = get_nd_simple_cubic_block(lsize, lblock, n_dim, static_center, periodic, get_size, edge_list);
  int64_t* largest_component = percolate_bond(&g, 1.0, &idxLambda); // bond-percolation given fixed site probability (Newman-Ziff method)
  float* expectation_value = get_expectation_value(sweep_prob, SWEEP_LEN, largest_component, g.num_edges, 0, g.nnode, true, false); // compute a bond-percolation curve
  r |= (expectation_value[0] > 0.08);
  r |= (expectation_value[1] < 0.73 || expectation_value[1] > 0.77);
  free(expectation_value);
  free(largest_component);
  free_graph(&g);

  sweep_prob[0] = 0.68;
  sweep_prob[1] = 0.74;
  periodic = false;
  get_size = false;
  g = get_nd_simple_cubic_block(lsize, lblock, n_dim, static_center, periodic, get_size, edge_list);
  largest_component = percolate_bond(&g, 1.0, &idxLambda); // bond-percolation given fixed site probability (Newman-Ziff method)
  expectation_value = get_expectation_value(sweep_prob, SWEEP_LEN, largest_component, g.num_edges, 0, 1, true, false); // compute a bond-percolation curve
  r |= (expectation_value[0] > 0.0001);
  r |= (expectation_value[1] < 0.9999);
  free(expectation_value);
  free(largest_component);
  free_graph(&g);

  return r;
}
