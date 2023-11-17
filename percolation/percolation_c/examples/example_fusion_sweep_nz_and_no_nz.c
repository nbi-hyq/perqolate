#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "inc/graph_type.h"
#include "inc/graph_construction.h"
#include "inc/percolation_main.h"

/* construct simple cubic graph with lossy fusion photons with and without modified Newman-Ziff method */
int main(){
  srand(time(NULL));
  int64_t lsize = 50; // size in each dimension
  uint8_t dimension = 3; // lattice dimension
  bool static_center = true; // central node (QD) cannot be lost
  bool periodic = true; // periodic boundaries
  bool get_size = true; // determine size of largest connected component
  int n_rep = 10; // number of repetitions for averaging
  int sweep_len = 50; // number of different probabilities
  float* sweep_prob; // probability that photon survives (tuned)
  bool edge_list = true; // keep track of edges in list (needed for analyse_fusion_net)
  sweep_prob = malloc(sweep_len * sizeof(float));
  for(int i=0;i<sweep_len;i++) sweep_prob[i] = 0.90 + 0.1*(0.5/sweep_len + 1.0*i/sweep_len);

  Graph g = get_nd_simple_cubic(lsize, dimension, static_center, periodic, get_size, edge_list);

  /* simulation with modified Newman-Ziff method */
  int64_t* largest_component_sum = NULL;
  int64_t numQbts; // 2 or more fusion qubits per edge (overall number of qubits is more than number of nodes)
  int64_t numStatic; // counts static nodes
  for(int rep=0; rep<n_rep; rep++){
    int64_t idxLambda; // not used here
    int64_t* largest_component = fusion_and_loss_nz(&g, 0.5, static_center, &numQbts, &numStatic, &idxLambda);
    if(rep==0){
      largest_component_sum = largest_component;
    } else{
      for(int64_t i=0; i<numQbts+1;i++) largest_component_sum[i] += largest_component[i];
      free(largest_component);
    }
  }
  float* expectation_value = get_expectation_value(sweep_prob, sweep_len, largest_component_sum, numQbts, numStatic, g.nnode * n_rep, true, false); // compute fusion/loss-percolation curve
  free(expectation_value);
  free(largest_component_sum);

  /* simulation without modified Newman-Ziff method */
  uint8_t nFusionPhoton = 2;
  float pFail = 0.5;
  for(int i=0;i<sweep_len;i++){
    int64_t largest_component_size = 0;
    for(int rep=0; rep<n_rep; rep++){
      largest_component_size += analyse_fusion_net(&g, sweep_prob[i], nFusionPhoton, pFail);
    }
    printf("%f %f\n", sweep_prob[i], (float)largest_component_size / n_rep / g.nnode);
  }

  free_graph(&g);
  free(sweep_prob);
  return 0;
}
