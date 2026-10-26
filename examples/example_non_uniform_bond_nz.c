#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "inc/graph_type.h"
#include "inc/graph_construction.h"
#include "inc/percolation_main.h"

/* bond-percolation simulation with non-uniform edge probabilities (Newman-Ziff-like method) */
int main(int argc, char *argv[]){
  float p_bnd_fix;
  if (argc > 1) p_bnd_fix = atof(argv[1]); // bond probability in direction transversal to the one where cluster spanning is checked
  else {printf("no probability specified!\n"); return -1;}

  if (argc > 2) srand(atoi(argv[2]));
  else srand(time(NULL));

  uint64_t num_rep; // repetitions for averaging
  if (argc > 3) num_rep = atoi(argv[3]);
  else num_rep = 100;

  /* define unit graph of simple cubic lattice */
  uint8_t nedge = 3; // number of different edges in unit graph representation
  uint8_t n_dim = 3; // number of dimensions

  uint8_t* blk_edge_label = malloc((size_t)nedge);
  for(uint8_t i=0; i<nedge; i++) blk_edge_label[i] = i;
  uint8_t* blk_edges = malloc((size_t)nedge * n_dim);
  memset(blk_edges, 0, (size_t)nedge * n_dim);
  int8_t* blk_vec = malloc((size_t)nedge * n_dim);
  memset(blk_vec, 0, (size_t)nedge * n_dim);
  blk_vec[0] = 1; blk_vec[4] = 1; blk_vec[8] = 1;

  UnitGraph u = new_unit_graph(blk_edges, blk_vec, blk_edge_label, nedge, n_dim); // percolation/cluster spanning will be checked along direction of first vector when using this construction
  free(blk_edges);
  free(blk_vec);
  free(blk_edge_label);

  /* define the actual graph */
  bool periodic = false; // no periodic boundaries for start-stop percolation
  bool get_size = false; // do define start/stop set
  bool edge_list = true; // track edges

  /* run start-stop bond percolation */
  float* pBond = malloc((size_t)nedge * sizeof(float));
  for(int64_t size=20; size<=160; size*=2){
    Graph g = get_lattice_from_unit_graph(&u, size, n_dim, false, periodic, get_size, edge_list);
    pBond[1] = p_bnd_fix; pBond[2] = p_bnd_fix; // bonds transversal to spanning direction are fixed
    pBond[0] = -1; // bonds in spanning direction are added one-by-one in Newman-Ziff style
    double threshold = 0;
    for(uint64_t r=0;r<num_rep;r++){
      threshold += percolate_bond_nonuniform_nz(&g, pBond);
    }
    printf("%li %.15f\n", size, threshold / num_rep);
    free_graph(&g);
  }

  free_unit_graph(&u);
  free(pBond);
  return 0;
}
