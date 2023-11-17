#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "inc/graph_type.h"
#include "inc/graph_construction.h"
#include "inc/percolation_main.h"

/* site-percolation simulation given fixed bond probability (using Newman-Ziff method)
   test for percolation from one side to the other */
int main(){
  srand(time(NULL));
  int64_t lsize = 1000; // size in each dimension
  bool periodic = false; // no periodic boundaries
  bool get_size = false; // check percolation, not determine size of largest connected component
  Graph g = get_2d_square(lsize, periodic, get_size, false);

  int sweep_len = 100; // number of different probabilities
  float* sweep_prob; // site probability (tuned)
  sweep_prob = malloc(sweep_len * sizeof(float));
  for(int i=0;i<sweep_len;i++) sweep_prob[i] = 0.7 + 0.1*(0.5/sweep_len + 1.0*i/sweep_len);

  int64_t idxLambda; // not used here
  int64_t* percolated = percolate_site(&g, 0.726195, &idxLambda); // site-percolation given fixed bond probability (Newman-Ziff method)
  float* expectation_value = get_expectation_value(sweep_prob, sweep_len, percolated, g.nnode, 0, 1, true, false); // compute a site-percolation curve

  free(expectation_value);
  free(percolated);
  free(sweep_prob);
  free_graph(&g);
  return 0;
}
