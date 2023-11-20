#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "inc/graph_type.h"
#include "inc/graph_construction.h"
#include "inc/percolation_main.h"

/* simulate size of largest component when applying photon loss to existing graph (no Newman-Ziff method) */
int main(){
  srand(time(NULL));
  int64_t lsize = 500; // size in each dimension
  bool get_size = true; // determine size of largest connected component
  bool periodic = false; // no periodic boundaries
  bool edge_list = true; // keep track of edges in list (needed for apply_loss)
  Graph g = get_2d_triangular(lsize, periodic, get_size, edge_list);

  int sweep_len = 200; // number of different probabilities
  float* sweep_prob; // site/bond probability (tuned)
  sweep_prob = malloc(sweep_len * sizeof(float)); // probability that photon survives (tuned)
  int64_t largest_component_size;
  for(int i=0;i<sweep_len;i++){
    sweep_prob[i] = 0.6 + 0.4*(0.5/sweep_len + 1.0*i/sweep_len);
    largest_component_size = apply_loss(&g, sweep_prob[i]);
    printf("%f %f\n", sweep_prob[i], (float)largest_component_size / g.nnode);
  }

  free(sweep_prob);
  free_graph(&g);
  return 0;
}
