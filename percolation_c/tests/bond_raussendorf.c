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

/* bond-percolation simulation of Raussendorf lattice (using Newman-Ziff method) */
int main(){
  srand(573467);
  int r = 0;
  int64_t lsize = 40; // size in each dimension
  bool periodic = false; // periodic boundaries
  bool static_center = false; // only meaningful if nodes/qubits can be lost
  bool get_size = false; // not determine size of largest connected component but check for cluster spanning
  bool edge_list = true; // keep track of edges in list (needed for bond percolation)
  Graph g = get_raussendorf_lattice(lsize, static_center, periodic, get_size, edge_list);

  float sweep_prob[SWEEP_LEN] = {0.380, 0.394};
  float* expectation_value; // relative size of largest connected component
  int64_t idxLambda;
  int64_t* largest_component = percolate_bond(&g, 1.0, &idxLambda); // bond-percolation given fixed site probability of 1.0 (Newman-Ziff method)
  expectation_value = get_expectation_value(sweep_prob, SWEEP_LEN, largest_component, g.num_edges, 0, 1, true, false); // compute a bond-percolation curve

  r |= (expectation_value[0] > 0.0001);
  r |= (expectation_value[1] < 0.9999);

  free(expectation_value);
  free(largest_component);
  free_graph(&g);
  return r;
}
