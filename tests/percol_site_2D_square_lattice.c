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
  srand(126738);
  int64_t lsize = 1000; // size in each dimension
  int r;
  Graph g = get_2d_square(lsize, false, false, false);

  float sweep_prob[SWEEP_LEN] = {0.70, 0.80}; // note: expected site-percolation threshold at 0.75 in this example
  int64_t idxLambda;
  int64_t* percolated = percolate_site(&g, 0.726195, &idxLambda); // specify (fixed) bond probability
  float* expectation_value = get_expectation_value(sweep_prob, SWEEP_LEN, percolated, g.nnode, 0, 1, true, true);  // compute a site-percolation curve (from site sweep)
  r = expectation_value[0] < 0.001 ? 0 : -1;
  r |= expectation_value[1] > 0.999 ? 0 : -1;
  free(expectation_value);
  free(percolated);
 
  free_graph(&g);
  return r;
}

