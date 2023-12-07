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
  srand(494839);
  int64_t lsize = 1000; // size in each dimension
  int r;
  Graph g = get_2d_square(lsize, true, true, false);

  float sweep_prob[SWEEP_LEN] = {0.70, 0.80}; // note: expected site-percolation threshold at 0.75 in this example
  int64_t idxLambda;
  int64_t* largest_component = percolate_site(&g, 0.726195, &idxLambda); // specify (fixed) bond probability
  float* expectation_value = get_expectation_value(sweep_prob, SWEEP_LEN, largest_component, g.nnode, 0, g.nnode, true, true);  // compute a site-percolation curve (from site sweep)
  r = expectation_value[0] < 0.1 ? 0 : -1;
  r |= expectation_value[1] > 0.5 ? 0 : -1;
  free(expectation_value);
  free(largest_component);
  free_graph(&g);

  // 3D simple cubic with additional diagonals (NN and 2NN), classical site percolation
  g = get_nd_simple_cubic_modified(40, 3, false, true, false, false, 2, 2, false, false, false);
  sweep_prob[0] = 0.11; sweep_prob[1] = 0.17; // transition at about 0.137 (https://doi.org/10.1016/S0034-4877(12)60036-6)
  int64_t* arryPercolated = percolate_site(&g, 1.0, &idxLambda);
  expectation_value = get_expectation_value(sweep_prob, SWEEP_LEN, arryPercolated, g.nnode, 0, 1, true, true); // compute fusion/loss-percolation curve
  if(expectation_value[0] > 0.5) r = -1; // below percolation threshold
  if(expectation_value[1] < 0.5) r |= -1; // above percolation threshold
  free_graph(&g);
  free(arryPercolated);
  free(expectation_value);

  // 4D simple cubic with only bcc diagonals, classical site percolation
  g = get_nd_simple_cubic_modified(20, 4, false, false, true, false, 0, 0, false, false, false);
  sweep_prob[0] = 0.07; sweep_prob[1] = 0.13; // transition at about 0.106 (https://doi.org/10.1142/S0129183198000431)
  arryPercolated = percolate_site(&g, 1.0, &idxLambda);
  expectation_value = get_expectation_value(sweep_prob, SWEEP_LEN, arryPercolated, g.nnode, 0, 1, true, true); // compute fusion/loss-percolation curve
  if(expectation_value[0] > 0.5) r = -1; // below percolation threshold
  if(expectation_value[1] < 0.5) r |= -1; // above percolation threshold
  free_graph(&g);
  free(arryPercolated);
  free(expectation_value);

  // 4D diamond lattice, classical site percolation
  g = get_nd_simple_cubic_modified(20, 4, false, true, false, true, 0, 0, false, false, false);
  sweep_prob[0] = 0.24; sweep_prob[1] = 0.33; // transition at about 0.3 (https://doi.org/10.1142/S0129183198000431)
  arryPercolated = percolate_site(&g, 1.0, &idxLambda);
  expectation_value = get_expectation_value(sweep_prob, SWEEP_LEN, arryPercolated, g.nnode, 0, 1, true, true); // compute fusion/loss-percolation curve
  if(expectation_value[0] > 0.5) r = -1; // below percolation threshold
  if(expectation_value[1] < 0.5) r |= -1; // above percolation threshold
  free_graph(&g);
  free(arryPercolated);
  free(expectation_value);

  return r;
}

