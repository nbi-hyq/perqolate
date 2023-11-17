#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"
#include "../inc/percolation_main.h"

#define SWEEP_LEN 2

/* reference percolation simulations (using Newman-Ziff method) */
int main(){
  srand(247822);
  int64_t lsize = 60; // size in each dimension
  bool periodic = false; // no periodic boundaries
  bool get_size = false; // check percolation, not determine size of largest connected component
  bool static_center = false;
  uint8_t dimension = 3;
  bool list_edges = true; // must be true for bond-percolation
  int r = 0; // return

  /* NN + 4NN cases from
     https://journals.aps.org/pre/abstract/10.1103/PhysRevE.91.043301 (site percolation)
     https://journals.aps.org/pre/abstract/10.1103/PhysRevE.102.012102 (bond percolation) */
  int64_t num_vec = 6;
  int64_t* vecs = malloc(num_vec * dimension * sizeof(int64_t));
  memset(vecs, 0, num_vec * dimension * sizeof(int64_t));
  vecs[dimension*0+0] = 1; // NN
  vecs[dimension*1+1] = 1; // NN
  vecs[dimension*2+2] = 1; // NN
  vecs[dimension*3+0] = 2; // 4NN
  vecs[dimension*4+1] = 2; // 4NN
  vecs[dimension*5+2] = 2; // 4NN

  Graph g = get_lattice_from_nd_simple_cubic_and_vectors(lsize, dimension, vecs, num_vec, static_center, periodic, get_size, list_edges);

  float sweep_prob[SWEEP_LEN] = {0.144, 0.157};
  int64_t idxLambda;
  int64_t* percolated = percolate_site(&g, 1.0, &idxLambda); // site-percolation given fixed bond probability (Newman-Ziff method)
  float* expectation_value = get_expectation_value(sweep_prob, SWEEP_LEN, percolated, g.nnode, 0, 1, true, false); // compute a site-percolation curve
  if(expectation_value[0] > 0.5) r = -1; // below percolation threshold
  if(expectation_value[1] < 0.5) r |= -1; // above percolation threshold
  printf("%f\n", (double)idxLambda/g.nnode);
  free(expectation_value);
  free(percolated);

  sweep_prob[0] = 0.104; sweep_prob[1] = 0.110;
  percolated = percolate_bond(&g, 1.0, &idxLambda); // bond-percolation given fixed site probability (Newman-Ziff method)
  expectation_value = get_expectation_value(sweep_prob, SWEEP_LEN, percolated, g.num_edges, 0, 1, true, false); // compute a bond-percolation curve
  if(expectation_value[0] > 0.5) r = -1; // below percolation threshold
  if(expectation_value[1] < 0.5) r |= -1; // above percolation threshold
  printf("%f\n", (double)idxLambda/g.num_edges);
  free(expectation_value);
  free(percolated);

  free_graph(&g);
  free(vecs);

  /* 3NN + 4NN cases from
     https://journals.aps.org/pre/abstract/10.1103/PhysRevE.91.043301 (site percolation)
     https://journals.aps.org/pre/abstract/10.1103/PhysRevE.102.012102 (bond percolation) */
  num_vec = 7;
  vecs = malloc(num_vec * dimension * sizeof(int64_t));
  memset(vecs, 0, num_vec * dimension * sizeof(int64_t));
  vecs[dimension*0+0] = 1; vecs[dimension*0+1] = 1; vecs[dimension*0+2] = 1; // 3NN
  vecs[dimension*1+0] = 1; vecs[dimension*1+1] = -1; vecs[dimension*1+2] = 1; // 3NN
  vecs[dimension*2+0] = -1; vecs[dimension*2+1] = 1; vecs[dimension*2+2] = 1; // 3NN
  vecs[dimension*3+0] = -1; vecs[dimension*3+1] = -1; vecs[dimension*3+2] = 1; // 3NN
  vecs[dimension*4+0] = 2; // 4NN
  vecs[dimension*5+1] = 2; // 4NN
  vecs[dimension*6+2] = 2; // 4NN

  g = get_lattice_from_nd_simple_cubic_and_vectors(lsize, dimension, vecs, num_vec, static_center, periodic, get_size, list_edges);

  /* //TBD: check again
  sweep_prob[0] = 0.202; sweep_prob[1] = 0.208;
  percolated = percolate_site(&g, 1.0, &idxLambda); // site-percolation given fixed bond probability (Newman-Ziff method)
  expectation_value = get_expectation_value(sweep_prob, SWEEP_LEN, percolated, g.nnode, 0, 1, true, false); // compute a site-percolation curve
  if(expectation_value[0] > 0.5) r = -1; // below percolation threshold
  if(expectation_value[1] < 0.5) r |= -1; // above percolation threshold
  printf("%f\n", (double)idxLambda/g.nnode);
  free(expectation_value);
  free(percolated);
  */

  sweep_prob[0] = 0.098; sweep_prob[1] = 0.104;
  percolated = percolate_bond(&g, 1.0, &idxLambda); // bond-percolation given fixed site probability (Newman-Ziff method)
  expectation_value = get_expectation_value(sweep_prob, SWEEP_LEN, percolated, g.num_edges, 0, 1, true, false); // compute a bond-percolation curve
  if(expectation_value[0] > 0.5) r = -1; // below percolation threshold
  if(expectation_value[1] < 0.5) r |= -1; // above percolation threshold
  printf("%f\n", (double)idxLambda/g.num_edges);
  free(expectation_value);
  free(percolated);

  free_graph(&g);
  free(vecs);

  /* 2NN + 4NN cases from
     https://journals.aps.org/pre/abstract/10.1103/PhysRevE.91.043301 (site percolation)
     https://journals.aps.org/pre/abstract/10.1103/PhysRevE.102.012102 (bond percolation) */
  num_vec = 9;
  vecs = malloc(num_vec * dimension * sizeof(int64_t));
  memset(vecs, 0, num_vec * dimension * sizeof(int64_t));
  vecs[dimension*0+0] = 1; vecs[dimension*0+1] = 0; vecs[dimension*0+2] = 1; // 2NN
  vecs[dimension*1+0] = -1; vecs[dimension*1+1] = 0; vecs[dimension*1+2] = 1; // 2NN
  vecs[dimension*2+0] = 0; vecs[dimension*2+1] = 1; vecs[dimension*2+2] = 1; // 2NN
  vecs[dimension*3+0] = 0; vecs[dimension*3+1] = -1; vecs[dimension*3+2] = 1; // 2NN
  vecs[dimension*4+0] = 1; vecs[dimension*4+1] = 1; vecs[dimension*4+2] = 0; // 2NN
  vecs[dimension*5+0] = 1; vecs[dimension*5+1] = -1; vecs[dimension*5+2] = 0; // 2NN
  vecs[dimension*6+0] = 2; // 4NN
  vecs[dimension*7+1] = 2; // 4NN
  vecs[dimension*8+2] = 2; // 4NN

  g = get_lattice_from_nd_simple_cubic_and_vectors(lsize, dimension, vecs, num_vec, static_center, periodic, get_size, list_edges);

  /* //TBD: check again
  sweep_prob[0] = 0.157; sweep_prob[1] = 0.163;
  percolated = percolate_site(&g, 1.0, &idxLambda); // site-percolation given fixed bond probability (Newman-Ziff method)
  expectation_value = get_expectation_value(sweep_prob, SWEEP_LEN, percolated, g.nnode, 0, 1, true, false); // compute a site-percolation curve
  if(expectation_value[0] > 0.5) r = -1; // below percolation threshold
  if(expectation_value[1] < 0.5) r |= -1; // above percolation threshold
  printf("%f\n", (double)idxLambda/g.nnode);
  free(expectation_value);
  free(percolated);
  */

  sweep_prob[0] = 0.072; sweep_prob[1] = 0.078;
  percolated = percolate_bond(&g, 1.0, &idxLambda); // bond-percolation given fixed site probability (Newman-Ziff method)
  expectation_value = get_expectation_value(sweep_prob, SWEEP_LEN, percolated, g.num_edges, 0, 1, true, false); // compute a bond-percolation curve
  if(expectation_value[0] > 0.5) r = -1; // below percolation threshold
  if(expectation_value[1] < 0.5) r |= -1; // above percolation threshold
  printf("%f\n", (double)idxLambda/g.num_edges);
  free(expectation_value);
  free(percolated);

  free_graph(&g);
  free(vecs);

  return r;
}
