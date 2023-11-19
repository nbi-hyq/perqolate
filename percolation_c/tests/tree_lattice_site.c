#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"
#include "../inc/percolation_main.h"

/* test: site percolation simulation of tree lattices */
int main(){
  srand(5689678);
  float pBond = 1.0; // all bonds are there: pure site-percolation
  bool static_center = false;
  bool edge_list = false;
  bool get_size = false;
  bool mute = true;
  int n_rep = 30; // number of repetitions for averaging
  int64_t idxLambda; // not used here
  int sweep_len = 100;
  float* sweep_prob = malloc(sweep_len * sizeof(float));
  for(int i=0; i < sweep_len;i++) sweep_prob[i] = 0.0 + 1.0*(0.5/sweep_len + 1.0*i/sweep_len);

  // Bethe lattice with vertex degree 4, expected site-percolation threshold: 1/(4-1)
  uint8_t depth_max = 9;
  uint8_t* num_child = malloc(depth_max+1);
  uint8_t vertex_degree = 4;
  num_child[0] = vertex_degree; // root node
  num_child[depth_max] = 0; // leave nodes have no children
  for(int i=1; i < depth_max;i++) num_child[i] = vertex_degree-1;

  Graph g = get_tree_lattice(depth_max, num_child, static_center, get_size, edge_list);

  int64_t* percolated_sum = NULL;
  int64_t* percolated;
  for(int rep=0; rep<n_rep; rep++){ //repeat everything for averaging
    percolated = percolate_site(&g, pBond, &idxLambda);
    if(rep==0){
     percolated_sum = percolated;
    } else {
      for(int64_t i=0; i<g.nnode+1;i++) percolated_sum[i] += percolated[i];
      free(percolated);
    }
  }
  float* expectation_value = get_expectation_value(sweep_prob, sweep_len, percolated_sum, g.nnode, 0, n_rep, true, mute);  // compute a site-percolation curve (from site sweep)

  int r = expectation_value[10] < 0.05 ? 0 : -1;
  r |= expectation_value[40] > 0.15 ? 0 : -1;
  free(expectation_value);
  free(num_child);
  free_graph(&g);

  // Bethe lattice with vertex degree 5, expected site-percolation threshold: 1/(5-1)
  depth_max = 8;
  num_child = malloc(depth_max+1);
  vertex_degree = 5;
  num_child[0] = vertex_degree; // root node
  num_child[depth_max] = 0; // leave nodes have no children
  for(int i=1; i < depth_max;i++) num_child[i] = vertex_degree-1;

  g = get_tree_lattice(depth_max, num_child, static_center, get_size, edge_list);

  for(int rep=0; rep<n_rep; rep++){ //repeat everything for averaging
    percolated = percolate_site(&g, pBond, &idxLambda);
    if(rep==0){
      free(percolated_sum);
      percolated_sum = percolated;
    } else {
      for(int64_t i=0; i<g.nnode+1;i++) percolated_sum[i] += percolated[i];
      free(percolated);
    }
  }
  expectation_value = get_expectation_value(sweep_prob, sweep_len, percolated_sum, g.nnode, 0, n_rep, true, mute);  // compute a site-percolation curve (from site sweep)

  r |= expectation_value[15] < 0.05 ? 0 : -1;
  r |= expectation_value[35] > 0.15 ? 0 : -1;
  free(expectation_value);
  free(percolated_sum);
  free(num_child);
  free_graph(&g);
  free(sweep_prob);

  return r;
}

