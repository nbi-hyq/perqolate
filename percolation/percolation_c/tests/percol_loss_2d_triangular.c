#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"
#include "../inc/percolation_main.h"

int main(){
  srand(629209);
  int64_t lsize = 200; // size in each dimension
  bool get_size = false; // determine size of largest connected component
  bool edge_list = true; // keep track of edges in list (needed for apply_loss)
  bool periodic = false; // no periodic boundaries
  int64_t percolated;
  int r = 0;
  Graph g = get_2d_triangular(lsize, periodic, get_size, edge_list);

  // below percolation threshold
  percolated = apply_loss(&g, 0.83);
  if(percolated > 0) r = -1;

  // above percolation threshold
  percolated = apply_loss(&g, 0.93);
  if(percolated < 1) r |= -1;

  // below percolation threshold
  percolated = apply_loss_bfs(&g, 0.83);
  if(percolated > 0) r = -1;

  // above percolation threshold
  percolated = apply_loss_bfs(&g, 0.93);
  if(percolated < 1) r |= -1;

  free_graph(&g);
  return r;
}
