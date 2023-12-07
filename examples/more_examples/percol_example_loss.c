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
  bool edge_list = true; // keep track of edges in list (needed for apply_loss)
  int sweep_len = 50; // number of different probabilities
  int rep = 10; // number of repetitions
  float sweep_prob; // 1 - loss probability (tuned)
  float percol_probability;
  Graph g = get_2d_triangular(lsize, periodic, get_size, edge_list);
  for(int i=0;i<sweep_len;i++){
    sweep_prob = 0.6 + 0.4*(0.5/sweep_len + 1.0*i/sweep_len);
    percol_probability = 0;
    for(int r=0;r<rep;r++){
      percol_probability += apply_loss(&g, sweep_prob);
    }
    printf("%f %f\n", sweep_prob, percol_probability / rep);
  }
  free_graph(&g);
  return 0;
}
