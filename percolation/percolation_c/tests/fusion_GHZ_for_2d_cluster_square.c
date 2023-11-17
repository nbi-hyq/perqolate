#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"
#include "../inc/percolation_main.h"

int main(){
  srand(192832);
  int64_t lsize = 500; // size in each dimension
  int64_t largest_component_size; // size of largest connected component
  int r = 0;
  Graph g = get_3qubit_GHZ_for_2d_square(lsize, false, true, true);

  // below percolation threshold
  largest_component_size = fusion_and_loss_no_newman_ziff(&g, 0.64, 1.0); // 0.672 is threshold, assume no photon loss
  if(8.0*(float)largest_component_size / g.nnode > 0.1) r = -1;

  // above percolation threshold
  largest_component_size = fusion_and_loss_no_newman_ziff(&g, 0.70, 1.0); // 0.672 is threshold, assume no photon loss
  if(8.0*(float)largest_component_size / g.nnode < 0.5) r |= -1;

  free_graph(&g);

  return r;
}

