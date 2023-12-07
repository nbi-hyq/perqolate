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
  srand(135679);

  // 4D diamond bond percolation
  Graph g = get_nd_simple_cubic_modified(35, 4, false, true, false, true, 0, 0, false, false, true);
  float sweep_prob[SWEEP_LEN] = {0.255, 0.285}; // 0.2715 expected from: https://doi.org/10.1142/S0129183198000431
  int64_t idxLambda;
  int64_t* percolated = percolate_bond(&g, 1.0, &idxLambda);
  float* expectation_value = get_expectation_value(sweep_prob, SWEEP_LEN, percolated, g.num_edges, 0, 1, true, true); // compute a bond-percolation curve (from bond sweep)
  // below percolation threshold
  int r = (expectation_value[0] > 0.1 ? -1 : 0);
  // above percolation threshold
  if(expectation_value[1] < 0.8) r |= -1;

  free_graph(&g);
  free(percolated);
  free(expectation_value);

  bool periodic = true;
  g = get_nd_simple_cubic_modified(10, 2, false, true, false, true, 0, 0, periodic, false, true);
  if(g.failed == true) r |= -1;
  free_graph(&g);

  g = get_nd_simple_cubic_modified(9, 2, false, true, false, true, 0, 0, periodic, false, true);
  if(g.failed == false) r |= -1; // periodic boundaries for diamond with uneven length does not work
  free_graph(&g);

  return r;
}

