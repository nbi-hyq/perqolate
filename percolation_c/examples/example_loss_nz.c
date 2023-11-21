#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "inc/graph_type.h"
#include "inc/graph_construction.h"
#include "inc/percolation_main.h"

/* simulate size of largest component when applying photon loss to existing graph (using modified Newman-Ziff method) */
int main(){
  srand(time(NULL));
  int64_t lsize = 500; // size in each dimension
  bool get_size = true; // determine size of largest connected component
  bool periodic = false; // no periodic boundaries
  Graph g = get_2d_triangular(lsize, periodic, get_size, false);

  int sweep_len = 200; // number of different probabilities
  float* sweep_prob; // probability that photon is not LOST (tuned)
  sweep_prob = malloc(sweep_len * sizeof(float));
  for(int i=0;i<sweep_len;i++) sweep_prob[i] = 0.6 + 0.4*(0.5/sweep_len + 1.0*i/sweep_len);

  float pBond = 1.0; // bond missing on graph before loss happens
  bool considerStaticNode = false; // cosider all qubits as lossy (ignoring g->static_node)
  int64_t numStatic; // counter of static nodes (not used here: could be set 0 and use NULL rather than &numStatic)
  int64_t idxLambda; // not used here
  int64_t* largest_component = apply_loss_nz(&g, pBond, considerStaticNode, &numStatic, &idxLambda); // size of largest component (number of nodes is position of array)
  float* expectation_value = expectation_value = get_expectation_value(sweep_prob, sweep_len, largest_component, g.nnode, numStatic, g.nnode, true, false); // compute a loss-percolation curve

  free(expectation_value);
  free(largest_component);
  free(sweep_prob);
  free_graph(&g);
  return 0;
}

