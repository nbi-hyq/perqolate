#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"
#include "../inc/percolation_main.h"

#define SWEEP_LEN 2

/* percolation simulation for NN + 2NN + 3NN case from (using Newman-Ziff method)
   bond-percolation: https://journals.aps.org/pre/abstract/10.1103/PhysRevE.102.012102
   site-percolation: https://doi.org/10.1016/S0034-4877(12)60036-6 */
int main(){
  srand(17284);
  int64_t lsize = 70; // size in each dimension
  bool periodic = false; // no periodic boundaries
  bool get_size = false; // check percolation, not determine size of largest connected component
  bool static_center = false;
  uint8_t dimension = 3;
  uint8_t max_step = 1; // maximum number of steps per dimension
  bool list_edges = true; // must be true for bond-percolation
  int r = 0; // return

  int64_t num_vec = 1;
  for(uint8_t j=0;j<dimension;j++) num_vec = num_vec * (2*(int64_t)max_step + 1);
  num_vec -= 1; // -1 because 0-vector is not considered
  num_vec /= 2; // if two vectors just differ by sign, just look at one of them
  int64_t* vecs = get_non_equivalent_vecs(dimension, max_step); // vectors specifying the lattice

  Graph g = get_lattice_from_nd_simple_cubic_and_vectors(lsize, dimension, vecs, num_vec, static_center, periodic, get_size, list_edges);

  /* site percolation */
  float sweep_prob[SWEEP_LEN] = {0.09, 0.11};
  int64_t idxLambda;
  int64_t* percolated = percolate_site(&g, 1.0, &idxLambda); // site-percolation given fixed bond probability (Newman-Ziff method)
  float* expectation_value = get_expectation_value(sweep_prob, SWEEP_LEN, percolated, g.nnode, 0, 1, true, false); // compute a site-percolation curve
  if(expectation_value[0] > 0.5) r = -1; // below percolation threshold
  if(expectation_value[1] < 0.5) r |= -1; // above percolation threshold
  printf("%f\n", (double)idxLambda/g.nnode);
  free(expectation_value);
  free(percolated);

  /* bond percolation */
  sweep_prob[0] = 0.045; sweep_prob[1] = 0.055;
  percolated = percolate_bond(&g, 1.0, &idxLambda); // bond-percolation given fixed site probability (Newman-Ziff method)
  expectation_value = get_expectation_value(sweep_prob, SWEEP_LEN, percolated, g.num_edges, 0, 1, true, false); // compute a bond-percolation curve
  if(expectation_value[0] > 0.5) r = -1; // below percolation threshold
  if(expectation_value[1] < 0.5) r |= -1; // above percolation threshold
  printf("%f\n", (double)idxLambda/g.num_edges);
  free(expectation_value);
  free(percolated);

  free(vecs);
  free_graph(&g);
  return r;
}
