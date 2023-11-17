#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"
#include "../inc/percolation_main.h"

int main(){
  srand(192083);
  int64_t lsize = 20; // size in each dimension
  int64_t largest_component_size; // size of largest connected component
  int r = 0;

  // without Newman-Ziff method (static center)
  Graph g = get_GHZ_stars_for_nd_simple_cubic(lsize, 3, true, true, true);
  // below percolation threshold
  largest_component_size = fusion_and_loss_no_newman_ziff(&g, 0.5, 0.92);
  if(7.0*(float)largest_component_size / g.nnode > 0.1) r = -1;
  // above percolation threshold
  largest_component_size = fusion_and_loss_no_newman_ziff(&g, 0.5, 0.97);
  if(7.0*(float)largest_component_size / g.nnode < 0.5) r |= -1;
  free_graph(&g);

  // with Newman-Ziff method (no static center)
  bool static_center = false; // nodes can be lost
  g = get_nd_simple_cubic(lsize, 3, static_center, true, true, false);
  int sweep_len = 2; // number of different probabilities
  float sweep_prob[2] = {0.93, 0.98}; // 1 - p_loss (tuned)
  float* expectation_value; // relative size of largest connected component
  int64_t numQbts; // assume 2 fusion qubits per edge -- overall number of qubits is more than number of nodes
  int64_t numStatic; // number of static nodes
  int64_t idxLambda;
  int64_t* largest_component = fusion_and_loss_nz(&g, 0.5, static_center, &numQbts, &numStatic, &idxLambda);
  expectation_value = get_expectation_value(sweep_prob, sweep_len, largest_component, numQbts, numStatic, g.nnode, true, true); // compute fusion/loss-percolation curve
  // below percolation threshold
  if(expectation_value[0] > 0.1) r = -1;
  // above percolation threshold
  if(expectation_value[1] < 0.5) r |= -1;
  free_graph(&g);
  free(largest_component);
  free(expectation_value);

  // with Newman-Ziff method (with static central node)
  static_center = true; // nodes cannot be lost
  g = get_nd_simple_cubic(lsize, 3, static_center, true, true, false);
  sweep_len = 2; // number of different probabilities
  sweep_prob[0] = 0.92; sweep_prob[1] = 0.97; // 1 - p_loss (tuned)
  largest_component = fusion_and_loss_nz(&g, 0.5, static_center, &numQbts, &numStatic, &idxLambda);
  expectation_value = get_expectation_value(sweep_prob, sweep_len, largest_component, numQbts, numStatic, g.nnode, true, true); // compute fusion/loss-percolation curve
  // below percolation threshold
  if(expectation_value[0] > 0.1) r = -1;
  // above percolation threshold
  if(expectation_value[1] < 0.5) r |= -1;
  free_graph(&g);
  free(largest_component);
  free(expectation_value);

  // simple cubic with additional diagonals, with Newman-Ziff method (with static central node)
  static_center = true; // nodes cannot be lost
  g = get_nd_simple_cubic_modified(50, 3, static_center, true, true, false, 0, 0, true, true, false);
  sweep_len = 2; // number of different probabilities
  sweep_prob[0] = 0.93; sweep_prob[1] = 0.96; // 1 - p_loss (tuned)
  largest_component = fusion_and_loss_nz(&g, 0.5, static_center, &numQbts, &numStatic, &idxLambda);
  expectation_value = get_expectation_value(sweep_prob, sweep_len, largest_component, numQbts, numStatic, g.nnode, true, true); // compute fusion/loss-percolation curve
  // below percolation threshold
  if(expectation_value[0] > 0.1) r = -1;
  // above percolation threshold
  if(expectation_value[1] < 0.2) r |= -1;
  free_graph(&g);
  free(largest_component);
  free(expectation_value);

  // 3D diamond with periodic boundaries, with Newman-Ziff method (with static central node)
  static_center = true; // nodes cannot be lost
  g = get_nd_simple_cubic_modified(40, 3, static_center, true, false, true, 0, 0, true, true, false);
  sweep_len = 2; // number of different probabilities
  sweep_prob[0] = 0.954; sweep_prob[1] = 0.974; // 1 - p_loss (tuned)
  largest_component = fusion_and_loss_nz(&g, 0.5, static_center, &numQbts, &numStatic, &idxLambda);
  expectation_value = get_expectation_value(sweep_prob, sweep_len, largest_component, numQbts, numStatic, g.nnode, true, true); // compute fusion/loss-percolation curve
  // below percolation threshold
  if(expectation_value[0] > 0.1) r = -1;
  // above percolation threshold
  if(expectation_value[1] < 0.2) r |= -1;
  free_graph(&g);
  free(largest_component);
  free(expectation_value);

  // try 3D diamond with periodic boundaries. Uneven sizes do not work
  static_center = true; // nodes cannot be lost
  g = get_nd_simple_cubic_modified(9, 3, static_center, true, false, true, 0, 0, true, true, false);
  if(g.failed != true) r |= -1;
  free_graph(&g);
  g = get_nd_simple_cubic_modified(10, 3, static_center, true, false, true, 0, 0, true, true, false);
  if(g.failed == true) r |= -1;
  free_graph(&g);
  g = get_nd_simple_cubic_modified(5, 3, static_center, true, false, true, 0, 0, true, false, false);
  if(g.failed != true) r |= -1;
  free_graph(&g);
  g = get_nd_simple_cubic_modified(6, 3, static_center, true, false, true, 0, 0, true, false, false);
  if(g.failed == true) r |= -1;
  free_graph(&g);

  return r;
}

