#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"
#include "../inc/percolation_main.h"

/* bond percolation threshold of syndrome graphs form 6-ring construction in https://www.nature.com/articles/s41467-023-36493-1*/
int main(){
  srand(792927);
  int n_avg = 1000; // number of repetitions for averaging
  bool periodic = false; // no periodic boundaries
  bool get_size = false; // check percolation, not determine size of largest connected component
  bool static_center = false;
  uint8_t dimension = 3;
  bool list_edges = true; // must be true for bond-percolation

  int64_t num_vec = 6;
  int64_t* vecs = malloc(num_vec * dimension * sizeof(int64_t));
  memset(vecs, 0, num_vec * dimension * sizeof(int64_t));
  vecs[dimension*0+0] = 1;
  vecs[dimension*1+1] = 1;
  vecs[dimension*2+2] = 1;
  vecs[dimension*3+0] = -1; vecs[dimension*3+1] = 1;
  vecs[dimension*4+1] = -1; vecs[dimension*4+2] = 1;
  vecs[dimension*5+0] = -1; vecs[dimension*5+2] = 1;

  for(int lsize=50; lsize<251; lsize+=50){
    Graph g = get_lattice_from_nd_simple_cubic_and_vectors(lsize, dimension, vecs, num_vec, static_center, periodic, get_size, list_edges);

    int64_t idxLambda;
    double avg = 0;
    for(int i=0; i < n_avg; i++) {
      int64_t* percolated = percolate_site(&g, 1.0, &idxLambda); // site-percolation given fixed bond probability (Newman-Ziff method)
      free(percolated);
      avg += (double)idxLambda/g.nnode;
    }
    printf("%f ", avg/n_avg);
    free_graph(&g);
  }

  free(vecs);

  return 0;
}
