#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "inc/graph_type.h"
#include "inc/graph_construction.h"
#include "inc/percolation_main.h"

/* construct square lattice graph with lossy fusion photons, repeat fusion if not successful (requires static_center == true), with Newman-Ziff method */
int main(){
  srand(time(NULL));
  int64_t lsize = 400; // size in each dimension
  uint8_t dimension = 2; // lattice dimension
  uint8_t fusion_attempts; // number of times a fusion attempt is made
  bool static_center = true; // central node (QD) cannot be lost
  bool periodic = false; // periodic boundaries
  bool get_size = false; // determine size of largest connected component
  int n_rep = 10; // number of repetitions for averaging
  int n_fusion_max = 5; // largest number for maximum number of fusion attempts simulated

  int sweep_len = 100; // number of different probabilities
  float* sweep_prob = malloc(sweep_len * sizeof(float)); // probability that photon survives (tuned)
  for(int i=0;i<sweep_len;i++) sweep_prob[i] = 0.9 + 0.1*(0.5/sweep_len + 1.0*i/sweep_len);

  float* expectation_value_sum = NULL;
  for(fusion_attempts=1;fusion_attempts<n_fusion_max;fusion_attempts++){
    Graph g = get_nd_simple_cubic(lsize, dimension, static_center, periodic, get_size, false);
    for(int rep=0; rep<n_rep; rep++){ // possibility to repeat everything for averaging
      int64_t numQbts; // 2 or more fusion qubits per edge (overall number of qubits is more than number of nodes)
      int64_t numStatic; // counts static nodes
      int64_t idxLambda;
      int64_t* largest_component = fusion_repeat_and_loss_nz(&g, 0.5, fusion_attempts, &numQbts, &numStatic, &idxLambda);
      printf("%f\n", ((double)idxLambda -(double)numStatic - 0.5) / ((double)numQbts -(double)numStatic));
      float* expectation_value = get_expectation_value(sweep_prob, sweep_len, largest_component, numQbts, numStatic, 1, true, true); // compute fusion/loss-percolation curve
      free(largest_component);
      if(rep==0){
        expectation_value_sum = expectation_value;
      } else{
        for(int i=0; i<sweep_len;i++) expectation_value_sum[i] += expectation_value[i]; // neccessary to average like this since numQbts is not constant here
        free(expectation_value);
      }
    }
    free_graph(&g);
    for(int k=0; k<sweep_len; k++) printf("%f %f\n", sweep_prob[k], expectation_value_sum[k] / n_rep);
    free(expectation_value_sum);
  }
  free(sweep_prob);
  return 0;
}

