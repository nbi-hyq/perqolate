#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "inc/graph_type.h"
#include "inc/graph_construction.h"
#include "inc/percolation_main.h"

/* construct graph by fusing corresponding GHZ-states with photon loss only on qubits used for the fusion
   Do the same thing twice, with/without Newman-Ziff method (consistency check) */
int main(){
  srand(time(NULL));
  
  // set parameters
  int64_t lsize = 50; // size in each dimension
  uint8_t dimension = 3; // lattice dimension
  bool static_center = true; // central node (photon) can be lost (not static)
  bool periodic = false; // no periodic boundaries
  bool get_size = false; // do not determine size of largest connected component (but check for percolation)
  int sweep_len = 50; // number of different probabilities
  int n_rep = 5; // number of repetitions
  float percol_probability; // probability of percolation
  float* sweep_prob; // 1 - p_loss (tuned)
  sweep_prob = malloc(sweep_len * sizeof(float));
  for(int i=0;i<sweep_len;i++) sweep_prob[i] = 0.90 + 0.1*(0.5/sweep_len + 1.0*i/sweep_len);

  // simulation without Newman-Ziff algorithm
  for(int i=0;i<sweep_len;i++){
    percol_probability = 0;
    for(int r=0;r<n_rep;r++){
      Graph g = get_GHZ_stars_for_nd_simple_cubic(lsize, dimension, static_center, periodic, get_size);
      percol_probability += fusion_and_loss_no_newman_ziff(&g, 0.5, sweep_prob[i]);
      free_graph(&g);
    }
    printf("%f %f\n", sweep_prob[i], percol_probability / n_rep); // average gives percolation probability
  }

  // same simulation with Newman-Ziff algorithm
  int64_t* percolated_sum = NULL;
  int64_t numQbts; // 2 fusion qubits per edge -- overall number of qubits is more than number of nodes
  int64_t numStatic; // counter of static nodes (not used here)
  Graph g = get_nd_simple_cubic(lsize, dimension, static_center, periodic, get_size, false);
  for(int r=0;r<n_rep;r++){
    int64_t idxLambda; // not used here
    int64_t* percolated = fusion_and_loss_nz(&g, 0.5, static_center, &numQbts, &numStatic, &idxLambda);
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
