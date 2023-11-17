#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"
#include "../inc/percolation_main.h"

int main(){
  srand(122418);
  int64_t lsize = 500; // size in each dimension
  int64_t percolated; // percolated: 1, not percolated: 0
  int r = 0;
  Graph g = get_3qubit_GHZ_for_2d_square(lsize, false, false, false);

  // below percolation threshold
  percolated = fusion_and_loss_no_newman_ziff(&g, 0.64, 1.0); // 0.672 is threshold, assume no photon loss
  if(percolated) r = -1;

  // above percolation threshold
  percolated = fusion_and_loss_no_newman_ziff(&g, 0.70, 1.0); // 0.672 is threshold, assume no photon loss
  if(!percolated) r |= -1;

  free_graph(&g);

  return r;
}

