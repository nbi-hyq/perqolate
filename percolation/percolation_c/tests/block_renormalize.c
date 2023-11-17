#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"
#include "../inc/percolation_main.h"

#define NUM_P 4

/* construct graph by fusion of micro-cluster states and do percolation simulation (repreat and tune size of overall graph), no Newman-Ziff method */
int main(){
  srand(5675879);
  int r = 0; // return
  uint8_t dimension = 3; // lattice dimension
  int64_t lsize = 10; // number of blocks in each direction
  int lblock = 5; // length of block in each dimension
  bool static_center = true; // central node (QD) cannot be lost

  float pSite[NUM_P] = {0.93, 0.97, 0.98, 0.99}; // efficiency (1 - p_loss)
  float boundLow[NUM_P] = {0.0, 0.4, 0.71, 0.92};
  float boundHigh[NUM_P] = {0.11, 0.5, 0.75, 0.95};
  uint8_t nFusionPhotons = 2;
  float pFail = 0.5;

  /* get size of renormalized lattice */
  bool periodic = true; // periodic boundaries
  bool get_size = true; // must be true for renormalization method used
  bool get_size_rn = true; // determine size of largest connected component (renormalized lattice)
  for(int i=0; i<NUM_P; i++){
    Graph g = get_nd_simple_cubic_block(lsize, lblock, dimension, static_center, periodic, get_size, true);
    double bondFraction;
    int64_t* bestNode = NULL;
    Graph g_rn = renormalize_block_fusion_net(&g, pSite[i], nFusionPhotons, pFail, lblock, lsize, dimension, periodic, get_size_rn, &bondFraction, &bestNode);
    r |= (bondFraction > boundLow[i] && bondFraction < boundHigh[i] ? 0 : 1);
    if(i==0){
      int64_t sizemax = bfs_full(&g_rn);
      r |= ((double)sizemax/(double)g_rn.nnode > 0.02);
    } else if(i==1){
      int64_t sizemax = bfs_full(&g_rn);
      r |= ((double)sizemax/(double)g_rn.nnode < 0.5);
    }
    free_graph(&g_rn);
    free_graph(&g);
    free(bestNode);
  }

  /* check cluster spanning of renormalized lattice */
  periodic = false; // no periodic boundaries
  get_size = true; // must be true for renormalization method used
  get_size_rn = false; // cluster spanning (renormalized lattice)
  for(int i=0; i<2; i++){
    Graph g = get_nd_simple_cubic_block(lsize, lblock, dimension, static_center, periodic, get_size, true);
    double bondFraction;
    int64_t* bestNode = NULL;
    Graph g_rn = renormalize_block_fusion_net(&g, pSite[i], nFusionPhotons, pFail, lblock, lsize, dimension, periodic, get_size_rn, &bondFraction, &bestNode);
    r |= (bondFraction > boundLow[i] && bondFraction < boundHigh[i] ? 0 : 1);
    int64_t percol = bfs_full(&g_rn);
    r |= (percol == i ? 0 : 1);
    free_graph(&g_rn);
    free_graph(&g);
    free(bestNode);
  }
  return r;
}

