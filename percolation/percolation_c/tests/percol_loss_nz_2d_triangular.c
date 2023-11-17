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
  srand(163838);
  int64_t lsize = 200; // size in each dimension
  bool get_size = false; // determine if percolated or not
  bool periodic = false; // no periodic boundaries
  int r;
  Graph g = get_2d_triangular(lsize, periodic, get_size, false);

  float sweep_prob[SWEEP_LEN] = {0.83, 0.94}; // note: threshold expected around 0.88
  int64_t numStatic; // number of static nodes (always 0 here)
  int64_t idxLambda;
  int64_t* percolated = apply_loss_nz(&g, 1.0, false, &numStatic, &idxLambda); // 1.0 bond probability, false: g->static_node would be ignored
  float* expectation_value = get_expectation_value(sweep_prob, SWEEP_LEN, percolated, g.nnode, numStatic, 1, true, true);  // compute a percolation curve
  r = expectation_value[0] < 0.1 ? 0 : -1;
  r |= expectation_value[1] > 0.5 ? 0 : -1;
  free(expectation_value);
  free(percolated);

  free_graph(&g);
  return r;
}
