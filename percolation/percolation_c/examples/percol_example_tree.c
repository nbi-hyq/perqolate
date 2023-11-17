#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"
#include "../inc/percolation_main.h"

/* site percolation simulation of tree lattices
   here: Bethe lattice with vertex degree 3, expected site-percolation threshold: 1/(3-1) */
int main(){
  srand(time(NULL));
  float pBond = 1.0; // all bonds are there: pure site-percolation
  bool static_center = false;
  bool edge_list = false;
  bool get_size = false;
  int64_t idxLambda; // not used here
  int sweep_len = 100;
  float* sweep_prob = malloc(sweep_len * sizeof(float));
  float* analytic = malloc(sweep_len * sizeof(float));
  for(int i=0; i < sweep_len;i++){
    sweep_prob[i] = 0.0 + 1.0*(0.5/sweep_len + 1.0*i/sweep_len);
    analytic[i] = sweep_prob[i] * (1 - (1 - sweep_prob[i])/sweep_prob[i] * (1 - sweep_prob[i])/sweep_prob[i] * (1 - sweep_prob[i])/sweep_prob[i]); // Eq. 16 in Stuffer/Aharony
    printf("%f %f\n", sweep_prob[i], analytic[i]);
  }

  uint8_t depth_max = 20;
  uint8_t* num_child = malloc(depth_max+1);
  uint8_t vertex_degree = 3;
  num_child[0] = vertex_degree; // root node
  num_child[depth_max] = 0; // leave nodes have no children
  for(int i=1; i < depth_max;i++) num_child[i] = vertex_degree-1;

  Graph g = get_tree_lattice(depth_max, num_child, static_center, get_size, edge_list);
  int64_t* percolated_sum = NULL;
  int64_t* percolated;
  int n_rep = 200;
  for(int rep=0; rep<n_rep; rep++){ //repeat everything for averaging
    percolated = percolate_site(&g, pBond, &idxLambda);
    if(rep==0){
     percolated_sum = percolated;
    } else {
      for(int64_t i=0; i<g.nnode+1;i++) percolated_sum[i] += percolated[i];
      free(percolated);
    }
  }
  float* expectation_value = get_expectation_value(sweep_prob, sweep_len, percolated_sum, g.nnode, 0, n_rep, true, false);  // compute a site-percolation curve (from site sweep)

  free(expectation_value);
  free(percolated_sum);
  free(sweep_prob);
  free(analytic);
  free(num_child);
  free_graph(&g);

  return 0;
}

