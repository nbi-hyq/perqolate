#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "inc/graph_type.h"
#include "inc/graph_construction.h"
#include "inc/percolation_main.h"

/* percolation simulation when applying photon loss to existing graph (no Newman-Ziff method) */
int main(){
  srand(time(NULL));
  int64_t lsize = 200; // size in each dimension
  bool get_size = false; // do not determine size of largest connected component (but check for percolation)
  bool periodic = false; // no periodic boundaries
  int sweep_len = 50; // number of different probabilities
  int n_rep = 10; // number of repetitions
  float* sweep_prob = malloc(sweep_len * sizeof(float)); // 1 - loss probability (tuned)
  Graph g = get_2d_triangular(lsize, periodic, get_size, false);
  for(int i=0;i<sweep_len;i++) sweep_prob[i] = 0.6 + 0.4*(0.5/sweep_len + 1.0*i/sweep_len);
  
  float pBond = 1.0; // bond missing on graph before loss happens
  bool considerStaticNode = false; // cosider all qubits as lossy (ignoring g->static_node)
  int64_t* percolated_sum = NULL;
  int64_t numStatic; // counter of static nodes (always 0 since considerStaticNode==false)
  for(int r=0;r<n_rep;r++){
    int64_t idxLambda; // number of existing photons where percolation appears
    int64_t* percolated = apply_loss_nz(&g, pBond, considerStaticNode, &numStatic, &idxLambda); // percolation or not (number of nodes is position of array)
    printf("%f\n", ((double)idxLambda -(double)numStatic - 0.5)/ ((double)g.nnode -(double)numStatic));
    if(r==0){
      percolated_sum = percolated;
    } else{
      for(int64_t i=0; i<g.nnode+1;i++) percolated_sum[i] += percolated[i];
      free(percolated);
    }
  }

  float* expectation_value = get_expectation_value(sweep_prob, sweep_len, percolated_sum, g.nnode, numStatic, n_rep, true, false); // compute a loss-percolation curve
  free(expectation_value);
  free(percolated_sum);
  free(sweep_prob);
  free_graph(&g);
  return 0;
}
