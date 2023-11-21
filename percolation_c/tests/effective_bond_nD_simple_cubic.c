#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"
#include "../inc/percolation_main.h"

/* bond percolation simulation of simple cubic lattice in several dimensions (by fusing GHZ-grphs with no photon loss) */
int main(){
  srand(247564);
  int64_t lsize; // size in each dimension
  bool static_center = true; // central node (QD) cannot be lost
  bool periodic = true; // periodic boundaries
  bool get_size = true; // determine size of largest connected component
  int64_t largest_component_size;
  uint8_t dimension; // lattice dimension
  int r = 0;
  
  dimension = 2;
  lsize = 100;
  Graph g = get_GHZ_stars_for_nd_simple_cubic(lsize, dimension, static_center, periodic, get_size);
  largest_component_size = fusion_and_loss_no_newman_ziff(&g, 0.4, 1.0); // below percolation threshold
  if((2.0*dimension+1) * (float)largest_component_size / g.nnode > 0.1) r = -1;
  largest_component_size = fusion_and_loss_no_newman_ziff(&g, 0.6, 1.0); // above percolation threshold
  if((2.0*dimension+1) * (float)largest_component_size / g.nnode < 0.5) r |= -1;
  free_graph(&g);

  dimension = 3;
  lsize = 20;
  g = get_GHZ_stars_for_nd_simple_cubic(lsize, dimension, static_center, periodic, get_size);
  largest_component_size = fusion_and_loss_no_newman_ziff(&g, 0.2, 1.0); // below percolation threshold
  if((2.0*dimension+1) * (float)largest_component_size / g.nnode > 0.1) r = -1;
  largest_component_size = fusion_and_loss_no_newman_ziff(&g, 0.3, 1.0); // above percolation threshold
  if((2.0*dimension+1) * (float)largest_component_size / g.nnode < 0.5) r |= -1;
  free_graph(&g);

  dimension = 4;
  lsize = 15;
  g = get_GHZ_stars_for_nd_simple_cubic(lsize, dimension, static_center, periodic, get_size);
  largest_component_size = fusion_and_loss_no_newman_ziff(&g, 0.13, 1.0); // below percolation threshold
  if((2.0*dimension+1) * (float)largest_component_size / g.nnode > 0.1) r = -1;
  largest_component_size = fusion_and_loss_no_newman_ziff(&g, 0.20, 1.0); // above percolation threshold
  if((2.0*dimension+1) * (float)largest_component_size / g.nnode < 0.5) r |= -1;
  free_graph(&g);

  dimension = 5;
  lsize = 10;
  g = get_GHZ_stars_for_nd_simple_cubic(lsize, dimension, static_center, periodic, get_size);
  largest_component_size = fusion_and_loss_no_newman_ziff(&g, 0.10, 1.0); // below percolation threshold
  if((2.0*dimension+1) * (float)largest_component_size / g.nnode > 0.1) r = -1;
  largest_component_size = fusion_and_loss_no_newman_ziff(&g, 0.15, 1.0); // above percolation threshold
  if((2.0*dimension+1) * (float)largest_component_size / g.nnode < 0.5) r |= -1;
  free_graph(&g);

  return r;
}
