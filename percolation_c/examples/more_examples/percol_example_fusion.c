#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "inc/graph_type.h"
#include "inc/graph_construction.h"
#include "inc/percolation_main.h"

/* construct graph by fusion of micro-cluster states and check for percolation in one dimension (no Newman_Ziff method) */
int main(){
  srand(time(NULL));
  int64_t lsize = 30; // size in each dimension
  uint8_t dimension = 3; // lattice dimension
  bool static_center = false; // central node (photon) can be lost (not static)
  bool periodic = false; // no periodic boundaries
  bool get_size = false; // do not determine size of largest connected component (but check for percolation)
  
  int sweep_len = 50; // number of different probabilities
  int rep = 10; // number of repetitions
  float sweep_prob; // probability that photon survives (tuned)
  float percol_probability; // probability of percolation
  for(int i=0;i<sweep_len;i++){
    sweep_prob = 0.9 + 0.1*(0.5/sweep_len + 1.0*i/sweep_len);
    percol_probability = 0;
    for(int r=0;r<rep;r++){
      Graph g = get_GHZ_stars_for_nd_simple_cubic(lsize, dimension, static_center, periodic, get_size);
      percol_probability += fusion_and_loss_no_newman_ziff(&g, 0.5, sweep_prob);
      free_graph(&g);
    }
    printf("%f %f\n", sweep_prob, percol_probability / rep); // average gives percolation probability
  }
  return 0;
}
