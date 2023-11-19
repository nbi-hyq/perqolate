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
  bool static_center = false; // central node (photon) can be lost (not static)
  bool periodic = false; // no periodic boundaries
  bool get_size = false; // do not determine size of largest connected component (but check for percolation)

  int n_rep = 10; // number of repetitions
  int sweep_len = 50; // number of different probabilities
  float* sweep_prob; // 1 - p_loss (tuned)
  sweep_prob = malloc(sweep_len * sizeof(float));
  for(int i=0;i<sweep_len;i++) sweep_prob[i] = 0.90 + 0.1*(0.5/sweep_len + 1.0*i/sweep_len);

  int64_t* percolated_sum = NULL;
  int64_t numQbts; // 2 fusion qubits per edge -- overall number of qubits is more than number of nodes
  int64_t numStatic; // counter of static nodes (not used here)
  Graph g = get_nd_simple_cubic(lsize, dimension, static_center, periodic, get_size, false);
  for(int r=0;r<n_rep;r++){
    int64_t idxLambda; // number of existing photons where percolation appears
    int64_t* percolated = fusion_and_loss_nz(&g, 0.5, static_center, &numQbts, &numStatic, &idxLambda);
    printf("%f\n", ((double)idxLambda -(double)numStatic - 0.5) / ((double)numQbts -(double)numStatic));
    if(r==0){
      percolated_sum = percolated;
    } else{
      for(int64_t i=0; i<numQbts+1;i++) percolated_sum[i] += percolated[i];
      free(percolated);
    }
  }

  float* expectation_value = get_expectation_value(sweep_prob, sweep_len, percolated_sum, numQbts, numStatic, n_rep, true, false); // compute fusion/loss-percolation curve
  free(expectation_value);
  free(percolated_sum);
  free(sweep_prob);
  free_graph(&g);
  return 0;
}
