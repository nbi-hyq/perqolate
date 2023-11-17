#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "inc/graph_type.h"
#include "inc/graph_construction.h"
#include "inc/percolation_main.h"

/* simulate size of largest component when applying photon loss to existing graph, sweep graph size (no Newman-Ziff method) */
int main(){
  srand(time(NULL));
  int64_t lsize; // size in each dimension
  bool get_size = true; // determine size of largest connected component
  bool periodic = false; // no periodic boundaries
  bool edge_list = true; // keep track of edges in list (needed for apply_loss)
  int sweep_len = 50; // number of different probabilities
  int n_rep = 30; // number of repetitions for averaging

  int64_t largest_component_size;
  float sweep_prob; // site/bond probability (tuned)
  for(lsize=100;lsize<201;lsize=lsize+100){
    Graph g = get_2d_triangular(lsize, periodic, get_size, edge_list);
    for(int rep=0;rep<n_rep;rep++){
      for(int i=0;i<sweep_len;i++){
        sweep_prob = 0.6 + 0.4*(0.5/sweep_len + 1.0*i/sweep_len);
        largest_component_size = apply_loss(&g, sweep_prob);
        printf("%f %f\n", sweep_prob, (float)largest_component_size / g.nnode);
      }
    }
    free_graph(&g);
  }
  return 0;
}
