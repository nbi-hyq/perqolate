#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"
#include "../inc/percolation_main.h"

#define SWEEP_LEN 2

int main(){
  srand(654729);
  
  // 3D bcc bond percolation (no simple cubic edges)
  Graph g = get_nd_simple_cubic_modified(50, 3, false, false, true, false, 0, 0, false, false, true);
  float sweep_prob[SWEEP_LEN] = {0.16, 0.20};
  int64_t idxLambda;
  int64_t* percolated = percolate_bond(&g, 1.0, &idxLambda);
  float* expectation_value = get_expectation_value(sweep_prob, SWEEP_LEN, percolated, g.num_edges, 0, 1, true, true); // compute a bond-percolation curve (from bond sweep)
  // below percolation threshold
  int r = (expectation_value[0] > 0.1 ? -1 : 0);
  // above percolation threshold
  if(expectation_value[1] < 0.9) r |= -1;
  free_graph(&g);
  free(percolated);
  free(expectation_value);

  return r;
}

