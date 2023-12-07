#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "inc/graph_type.h"
#include "inc/graph_construction.h"
#include "inc/percolation_main.h"

/* construct graph by fusion of micro-cluster states and do percolation simulation (no Newman-Ziff method) */
int main(){
  srand(time(NULL));
  int64_t lsize = 20; // size in each dimension
  uint8_t dimension = 3; // lattice dimension
  bool static_center = false; // central node (photon) can be lost (not static)
  bool periodic = true; // periodic boundaries
  bool get_size = true; // determine size of largest connected component
  Graph g = get_GHZ_stars_for_nd_simple_cubic(lsize, dimension, static_center, periodic, get_size);

  int sweep_len = 50; // number of different probabilities
  float* sweep_prob; // probability that photon survives (tuned)
  sweep_prob = malloc(sweep_len * sizeof(float));
  int64_t largest_component_size;
  for(int i=0;i<sweep_len;i++){
    sweep_prob[i] = 0.90 + 0.1*(0.5/sweep_len + 1.0*i/sweep_len);
    largest_component_size = fusion_and_loss_no_newman_ziff(&g, 0.5, sweep_prob[i]);
    printf("%f %f\n", sweep_prob[i], (2.0*dimension+1) * (float)largest_component_size / g.nnode);
  }

  free_graph(&g);
  free(sweep_prob);
  return 0;
}

