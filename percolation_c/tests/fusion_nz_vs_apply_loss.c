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

/* construct simple cubic graph with lossy fusion photons with and without modified Newman-Ziff method */
int main(){
  srand(235679);
  int r = 0; // return

  int64_t lsize = 30; // size in each dimension
  uint8_t dimension = 3; // lattice dimension
  bool static_center = true; // central node (QD) cannot be lost
  bool periodic = true; // periodic boundaries
  bool get_size = true; // determine size of largest connected component
  int n_rep = 5; // number of repetitions for averaging
  float sweep_prob[SWEEP_LEN] = {0.93, 0.95}; // probability that photon survives (tuned)
  bool edge_list = true; // keep track of edges in list (needed for analyse_fusion_net)

  Graph g = get_nd_simple_cubic(lsize, dimension, static_center, periodic, get_size, edge_list);

  /* largest connected component: simulation with modified Newman-Ziff method */
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
  float* expectation_value = get_expectation_value(sweep_prob, SWEEP_LEN, largest_component_sum, numQbts, numStatic, g.nnode * n_rep, true, false); // compute fusion/loss-percolation curve
  r |= (expectation_value[0] > 0.02);
  r |= (expectation_value[1] < 0.2 || expectation_value[1] > 0.5);
  free(expectation_value);
  free(largest_component_sum);

  /* largest connected component: simulation without modified Newman-Ziff method */
  uint8_t nFusionPhoton = 2;
  float pFail = 0.5;
  for(int i=0;i<SWEEP_LEN;i++){
    int64_t largest_component_size = 0;
    for(int rep=0; rep<n_rep; rep++){
      largest_component_size += analyse_fusion_net(&g, sweep_prob[i], nFusionPhoton, pFail);
    }
    float gcc = (float)largest_component_size / n_rep / g.nnode;
    if(i==0) r |= (gcc > 0.02);
    else r |= (gcc < 0.2 || gcc > 0.5);
  }

  /* rebuild graph for cluster spaning */
  free_graph(&g);
  periodic = false; // no periodic boundaries
  get_size = false; // cluster spanning
  g = get_nd_simple_cubic(lsize, dimension, static_center, periodic, get_size, edge_list);

  /* cluster spanning: simulation with modified Newman-Ziff method */
  int64_t* percol_sum = NULL;
  for(int rep=0; rep<n_rep; rep++){
    int64_t idxLambda; // not used here
    int64_t* largest_component = fusion_and_loss_nz(&g, 0.5, static_center, &numQbts, &numStatic, &idxLambda);
    if(rep==0){
      percol_sum = largest_component;
    } else{
      for(int64_t i=0; i<numQbts+1;i++) percol_sum[i] += largest_component[i];
      free(largest_component);
    }
  }
  expectation_value = get_expectation_value(sweep_prob, SWEEP_LEN, percol_sum, numQbts, numStatic, n_rep, true, false); // compute fusion/loss-percolation curve
  r |= (expectation_value[0] > 0.2);
  r |= (expectation_value[1] < 0.3);
  free(expectation_value);
  free(percol_sum);

  /* cluster spaning: simulation without modified Newman-Ziff method */
  nFusionPhoton = 2;
  pFail = 0.5;
  for(int i=0;i<SWEEP_LEN;i++){
    int64_t percol = 0;
    for(int rep=0; rep<n_rep; rep++){
      percol += analyse_fusion_net(&g, sweep_prob[i], nFusionPhoton, pFail);
    }
    if(i==0) r |= ((float)percol / n_rep > 0.2);
    else r |= ((float)percol / n_rep < 0.3);
  }

  free_graph(&g);
  return r;
}
