#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "inc/graph_type.h"
#include "inc/graph_construction.h"
#include "inc/percolation_main.h"

/* bond-percolation simulation given fixed site probability (using Newman-Ziff method) */
int main(){
  srand(time(NULL));
  int64_t lsize = 1000; // size in each dimension
  bool periodic = true; // periodic boundaries
  bool get_size = true; // determine size of largest connected component
  bool edge_list = true; // keep track of edges in list (needed for bond percolation)
  Graph g = get_2d_square(lsize, periodic, get_size, edge_list);

  int sweep_len = 100; // number of different probabilities
  float* sweep_prob; // bond probability (tuned)
  sweep_prob = malloc(sweep_len * sizeof(float));
  for(int i=0;i<sweep_len;i++) sweep_prob[i] = 0.5/sweep_len + 1.0*i/sweep_len;
  float* expectation_value; // relative size of largest connected component

  int64_t idxLambda;
  int64_t* largest_component = percolate_bond(&g, 0.75, &idxLambda); // bond-percolation given fixed site probability (Newman-Ziff method)
  expectation_value = get_expectation_value(sweep_prob, sweep_len, largest_component, g.num_edges, 0, g.nnode, true, false); // compute a bond-percolation curve

  free(expectation_value);
  free(largest_component);
  free(sweep_prob);
  free_graph(&g);
  return 0;
}
