#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "inc/graph_type.h"
#include "inc/graph_construction.h"
#include "inc/percolation_main.h"

/* construct graph by fusing corresponding GHZ-states with photon loss (using Newman-Ziff method) */
int main(){
  srand(time(NULL));
  int64_t lsize = 30; // size in each dimension
  uint8_t dimension = 3; // lattice dimension
  bool static_center = false; // nodes can be lost
  bool periodic = true; // periodic boundaries
  bool get_size = true; // determine size of largest connected component
  Graph g = get_nd_simple_cubic(lsize, dimension, static_center, periodic, get_size, false);

  int sweep_len = 50; // number of different probabilities
  float* sweep_prob; // 1 - p_loss (tuned)
  sweep_prob = malloc(sweep_len * sizeof(float));
  for(int i=0;i<sweep_len;i++) sweep_prob[i] = 0.90 + 0.1*(0.5/sweep_len + 1.0*i/sweep_len);

  int64_t numQbts; // 2 fusion qubits per edge -- overall number of qubits is more than number of nodes
  int64_t numStatic; // counts static nodes (is set 0 when static_center==false)
  int64_t idxLambda;
  int64_t* largest_component = fusion_and_loss_nz(&g, 0.5, static_center, &numQbts, &numStatic, &idxLambda);
  float* expectation_value = get_expectation_value(sweep_prob, sweep_len, largest_component, numQbts, numStatic, g.nnode, true, false); // compute fusion/loss-percolation curve

  free(expectation_value);
  free(largest_component);
  free(sweep_prob);
  free_graph(&g);
  return 0;
}
