#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"
#include "../inc/percolation_main.h"

int main(){
  srand(73628);
  int64_t lsize = 200; // size in each dimension
  bool get_size = true; // determine size of largest connected component
  bool periodic = false; // no periodic boundaries
  bool edge_list = true; // keep track of edges in list (needed for apply_loss)
  int64_t largest_component_size;
  int r = 0;
  Graph g = get_2d_triangular(lsize, periodic, get_size, edge_list);

  // below percolation threshold
  largest_component_size = apply_loss(&g, 0.83);
  if((float)largest_component_size / g.nnode > 0.1) r = -1;

  // above percolation threshold
  largest_component_size = apply_loss(&g, 0.94);
  if((float)largest_component_size / g.nnode < 0.5) r |= -1;

  g.fusion_node = malloc(g.nnode * sizeof(bool)); // required for apply_loss_bfs

  // below percolation threshold
  largest_component_size = apply_loss_bfs(&g, 0.83);
  if((float)largest_component_size / g.nnode > 0.1) r = -1;

  // above percolation threshold
  largest_component_size = apply_loss_bfs(&g, 0.94);
  if((float)largest_component_size / g.nnode < 0.5) r |= -1;

  free_graph(&g);
  return r;
}
