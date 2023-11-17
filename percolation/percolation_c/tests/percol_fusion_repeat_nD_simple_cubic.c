#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"
#include "../inc/percolation_main.h"

int main(){
  srand(37381);
  int64_t lsize = 300; // size in each dimension
  uint8_t dimension = 2; // lattice dimension
  uint8_t fusion_attempts; // number of times a fusion attempt is made
  bool static_center = true; // central node (QD) cannot be lost
  bool periodic = false; // no periodic boundaries
  bool get_size = false; // check for percolation
  int64_t numQbts; // 2 or more fusion qubits per edge (overall number of qubits is more than number of nodes)
  int64_t numStatic; // counts static nodes
  int sweep_len = 2; // number of different probabilities

  Graph g = get_nd_simple_cubic(lsize, dimension, static_center, periodic, get_size, false);

  // 1 fusion attempt
  fusion_attempts = 1;
  float sweep_prob[2] = {0.97, 1.0}; // 1 - p_loss (tuned)
  int64_t idxLambda;
  int64_t* largest_component = fusion_repeat_and_loss_nz(&g, 0.5, fusion_attempts, &numQbts, &numStatic, &idxLambda);
  float* expectation_value = get_expectation_value(sweep_prob, sweep_len, largest_component, numQbts, numStatic, 1, true, false); // compute fusion/loss-percolation curve
  // below percolation threshold
  if(expectation_value[0] > 0.1) return -1;
  free(expectation_value);
  free(largest_component);
 
  // 2 fusion attempts
  fusion_attempts = 2;
  sweep_prob[0] = 0.94; sweep_prob[1] = 0.99; // 1 - p_loss (tuned)
  largest_component = fusion_repeat_and_loss_nz(&g, 0.5, fusion_attempts, &numQbts, &numStatic, &idxLambda);
  expectation_value = get_expectation_value(sweep_prob, sweep_len, largest_component, numQbts, numStatic, 1, true, false); // compute fusion/loss-percolation curve
  // below percolation threshold
  if(expectation_value[0] > 0.1) return -1;
  // above percolation threshold
  if(expectation_value[1] < 0.5) return -1;
  free(expectation_value);
  free(largest_component);

  // 3 fusion attempts
  fusion_attempts = 3;
  sweep_prob[0] = 0.94; sweep_prob[1] = 0.99; // 1 - p_loss (tuned)
  largest_component = fusion_repeat_and_loss_nz(&g, 0.5, fusion_attempts, &numQbts, &numStatic, &idxLambda);
  expectation_value = get_expectation_value(sweep_prob, sweep_len, largest_component, numQbts, numStatic, 1, true, false); // compute fusion/loss-percolation curve
  // below percolation threshold
  if(expectation_value[0] > 0.1) return -1;
  // above percolation threshold
  if(expectation_value[1] < 0.5) return -1;
  free(expectation_value);
  free(largest_component);
  free_graph(&g);

  return 0;
}

