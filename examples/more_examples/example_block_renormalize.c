#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "inc/graph_type.h"
#include "inc/graph_construction.h"
#include "inc/percolation_main.h"

/* construct graph by fusion of micro-cluster states, lattice renormalization (repreat and tune size of overall graph), no Newman-Ziff method */
int main(){
  srand(time(NULL));
  uint8_t dimension = 3; // lattice dimension
  int64_t lsize = 10; // number of blocks in each direction
  int lblock; // length of block in each dimension
  bool static_center = true; // central node (QD) cannot be lost

  float pSite; // efficiency (1 - p_loss)
  int num_pSite = 50;
  uint8_t nFusionPhotons = 2;
  float pFail = 0.5;

  bool periodic = true; // periodic boundaries
  bool get_size = true; // must be true for renormalization method used
  bool get_size_rn = true; // determine size of largest connected component (renomralized lattice)
  for(lblock=2; lblock<10; lblock++){
    for(int i=0; i<num_pSite; i++){
      pSite = 0.92 + 0.08*(0.5/num_pSite + 1.0*i/num_pSite);
      printf("%f ", pSite);
      Graph g = get_nd_simple_cubic_block(lsize, lblock, dimension, static_center, periodic, get_size, true);
      double bondFraction;
      int64_t* bestNode = NULL;
      Graph g_rn = renormalize_block_fusion_net(&g, pSite, nFusionPhotons, pFail, lblock, lsize, dimension, periodic, get_size_rn, &bondFraction, &bestNode);
      printf("  %f ", bondFraction); // fraction of existing bonds in renormalized lattice
      int64_t sizemax = bfs_full(&g_rn);
      printf("  %f\n", (double)sizemax/(double)g_rn.nnode); // relative size of largest component in renormalized lattice
      free_graph(&g_rn);
      free_graph(&g);
      free(bestNode);
    }
  }

  periodic = false; // no periodic boundaries for cluster spanning
  get_size = true; // must be true for renormalization method used
  get_size_rn = false; // check cluster spanning of renormalized lattice (renomralized lattice)
  for(lblock=9; lblock<10; lblock++){
    for(int i=0; i<num_pSite; i++){
      pSite = 0.92 + 0.08*(0.5/num_pSite + 1.0*i/num_pSite); 
      printf("%f ", pSite);
      Graph g = get_nd_simple_cubic_block(lsize, lblock, dimension, static_center, periodic, get_size, true);
      double bondFraction;
      int64_t* bestNode = NULL;
      Graph g_rn = renormalize_block_fusion_net(&g, pSite, nFusionPhotons, pFail, lblock, lsize, dimension, periodic, get_size_rn, &bondFraction, &bestNode);
      printf("  %f ", bondFraction); // fraction of existing bonds in renormalized lattice
      int64_t percol = bfs_full(&g_rn);
      printf("  %li\n", percol); // relative size of largest component in renormalized lattice
      free_graph(&g_rn);
      free_graph(&g);
      free(bestNode);
    }
  }
  return 0;
}

