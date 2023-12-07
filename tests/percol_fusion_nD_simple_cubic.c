#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"
#include "../inc/percolation_main.h"

int main(){
  srand(163893);
  int64_t lsize = 20; // size in each dimension
  int64_t percolated; // percolated: 1, not percolated: 0
  bool staticCenter;
  int r = 0;

  // no Newman-Ziff method, static nodes
  staticCenter = true;
  Graph g = get_GHZ_stars_for_nd_simple_cubic(lsize, 3, staticCenter, false, false);
  // below percolation threshold
  percolated = fusion_and_loss_no_newman_ziff(&g, 0.5, 0.92);
  if(percolated) r = -1;
  // above percolation threshold
  percolated = fusion_and_loss_no_newman_ziff(&g, 0.5, 0.97);
  if(!percolated) r |= -1;
  free_graph(&g);

  // with Newman-Ziff method, no static nodes
  staticCenter = false; // nodes can be lost
  g = get_nd_simple_cubic(lsize, 3, staticCenter, false, false, false);
  int sweep_len = 2; // number of different probabilities
  float sweep_prob[2] = {0.93, 0.98}; // 1 - p_loss (tuned)
  float* expectation_value; // estimated percolation probability
  int64_t numQbts; // assume 2 fusion qubits per edge -- overall number of qubits is more than number of nodes
  int64_t numStatic; // number of static nodes
  int64_t idxLambda;
  int64_t* arryPercolated = fusion_and_loss_nz(&g, 0.5, staticCenter, &numQbts, &numStatic, &idxLambda);
  expectation_value = get_expectation_value(sweep_prob, sweep_len, arryPercolated, numQbts, numStatic, 1, true, true); // compute fusion/loss-percolation curve
  // below percolation threshold
  if(expectation_value[0] > 0.1) r = -1;
  // above percolation threshold
  if(expectation_value[1] < 0.5) r |= -1;
  free_graph(&g);
  free(arryPercolated);
  free(expectation_value);

  // with Newman-Ziff method, with static nodes
  staticCenter = true; // central nodes cannot be lost
  g = get_nd_simple_cubic(lsize, 3, staticCenter, false, false, false);
  sweep_prob[0] = 0.92; sweep_prob[1] = 0.97; // 1 - p_loss (tuned)
  arryPercolated = fusion_and_loss_nz(&g, 0.5, staticCenter, &numQbts, &numStatic, &idxLambda);
  expectation_value = get_expectation_value(sweep_prob, sweep_len, arryPercolated, numQbts, numStatic, 1, true, true); // compute fusion/loss-percolation curve
  // below percolation threshold
  if(expectation_value[0] > 0.1) r = -1;
  // above percolation threshold
  if(expectation_value[1] < 0.5) r |= -1;
  free_graph(&g);
  free(arryPercolated);
  free(expectation_value);

  // 3D simple cubic with additional diagonals (NN and 3NN), with Newman-Ziff method, with static nodes
  staticCenter = true; // central nodes cannot be lost
  g = get_nd_simple_cubic_modified(30, 3, staticCenter, true, true, false, 0, 0, false, false, false);
  sweep_prob[0] = 0.93; sweep_prob[1] = 0.96; // 1 - p_loss (tuned)
  arryPercolated = fusion_and_loss_nz(&g, 0.5, staticCenter, &numQbts, &numStatic, &idxLambda);
  expectation_value = get_expectation_value(sweep_prob, sweep_len, arryPercolated, numQbts, numStatic, 1, true, true); // compute fusion/loss-percolation curve
  // below percolation threshold
  if(expectation_value[0] > 0.1) r = -1;
  // above percolation threshold
  if(expectation_value[1] < 0.5) r |= -1;
  free_graph(&g);
  free(arryPercolated);
  free(expectation_value);

  // 3D diamond lattice based on simple cubic, with Newman-Ziff method, without static nodes
  staticCenter = false; // central nodes can be lost
  bool diamond = true;
  g = get_nd_simple_cubic_modified(30, 3, staticCenter, true, false, diamond, 0, 0, false, false, false);
  sweep_prob[0] = 0.963; sweep_prob[1] = 0.983; // 1 - p_loss (tuned)
  arryPercolated = fusion_and_loss_nz(&g, 0.5, staticCenter, &numQbts, &numStatic, &idxLambda);
  expectation_value = get_expectation_value(sweep_prob, sweep_len, arryPercolated, numQbts, numStatic, 1, true, true); // compute fusion/loss-percolation curve
  // below percolation threshold
  if(expectation_value[0] > 0.1) r = -1;
  // above percolation threshold
  if(expectation_value[1] < 0.5) r |= -1;
  free_graph(&g);
  free(arryPercolated);
  free(expectation_value);

  return r;
}

