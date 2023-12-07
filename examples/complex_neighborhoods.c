#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"
#include "../inc/percolation_main.h"

/* reference percolation simulations (using Newman-Ziff method)
   https://journals.aps.org/pre/abstract/10.1103/PhysRevE.91.043301 (site percolation)
   https://journals.aps.org/pre/abstract/10.1103/PhysRevE.102.012102 (bond percolation) */
int main(){
  srand(234783);
  int64_t lsize = 249; // size in each dimension
  int n_avg = 10; // number of repetitions for averaging
  bool periodic = false; // no periodic boundaries
  bool get_size = false; // check percolation, not determine size of largest connected component
  bool static_center = false;
  uint8_t dimension = 3;
  bool list_edges = true; // must be true for bond-percolation

  /* NN + 4NN */
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

  int64_t idxLambda;
  double avg = 0;
  for(int i=0; i < n_avg; i++) {
    int64_t* percolated = percolate_site(&g, 1.0, &idxLambda); // site-percolation given fixed bond probability (Newman-Ziff method)
    free(percolated);
    avg += (double)idxLambda/g.nnode;
  }
  printf("%f ", avg/n_avg);
  printf("%f ", avg/n_avg/0.15040);

  avg = 0;
  for(int i=0; i < n_avg; i++) {
    int64_t* percolated = percolate_bond(&g, 1.0, &idxLambda); // bond-percolation given fixed site probability (Newman-Ziff method)
    free(percolated);
    avg += (double)idxLambda/g.num_edges;
  }
  printf("%f ", avg/n_avg);
  printf("%f ", avg/n_avg/0.1068263);
  printf("NN + 4NN (site/bond)\n");

  free_graph(&g);
  free(vecs);

  /* 3NN + 4NN */
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

  avg = 0;
  for(int i=0; i < n_avg; i++) {
    int64_t* percolated = percolate_site(&g, 1.0, &idxLambda); // site-percolation given fixed bond probability (Newman-Ziff method)
    free(percolated);
    avg += (double)idxLambda/g.nnode;
  }
  printf("%f ", avg/n_avg);
  printf("%f ", avg/n_avg/0.20490);

  avg = 0;
  for(int i=0; i < n_avg; i++) {
    int64_t* percolated = percolate_bond(&g, 1.0, &idxLambda); // bond-percolation given fixed site probability (Newman-Ziff method)
    free(percolated);
    avg += (double)idxLambda/g.num_edges;
  }
  printf("%f ", avg/n_avg);
  printf("%f ", avg/n_avg/0.1012133);
  printf("3NN + 4NN (site/bond)\n");

  free_graph(&g);
  free(vecs);

  /* NN + 3NN */
  num_vec = 7;
  vecs = malloc(num_vec * dimension * sizeof(int64_t));
  memset(vecs, 0, num_vec * dimension * sizeof(int64_t));
  vecs[dimension*0+0] = 1; // NN
  vecs[dimension*1+1] = 1; // NN
  vecs[dimension*2+2] = 1; // NN
  vecs[dimension*3+0] = 1; vecs[dimension*3+1] = 1; vecs[dimension*3+2] = 1; // 3NN
  vecs[dimension*4+0] = 1; vecs[dimension*4+1] = -1; vecs[dimension*4+2] = 1; // 3NN
  vecs[dimension*5+0] = -1; vecs[dimension*5+1] = 1; vecs[dimension*5+2] = 1; // 3NN
  vecs[dimension*6+0] = -1; vecs[dimension*6+1] = -1; vecs[dimension*6+2] = 1; // 3NN

  g = get_lattice_from_nd_simple_cubic_and_vectors(lsize, dimension, vecs, num_vec, static_center, periodic, get_size, list_edges);

  avg = 0;
  for(int i=0; i < n_avg; i++) {
    int64_t* percolated = percolate_site(&g, 1.0, &idxLambda); // site-percolation given fixed bond probability (Newman-Ziff method)
    free(percolated);
    avg += (double)idxLambda/g.nnode;
  }
  printf("%f ", avg/n_avg);
  printf("%f ", avg/n_avg/0.1420);

  avg = 0;
  for(int i=0; i < n_avg; i++) {
    int64_t* percolated = percolate_bond(&g, 1.0, &idxLambda); // bond-percolation given fixed site probability (Newman-Ziff method)
    free(percolated);
    avg += (double)idxLambda/g.num_edges;
  }
  printf("%f ", avg/n_avg);
  printf("%f ", avg/n_avg/0.0920213);
  printf("NN + 3NN (site/bond)\n");

  free_graph(&g);
  free(vecs);

  /* NN + 2NN */
  num_vec = 9;
  vecs = malloc(num_vec * dimension * sizeof(int64_t));
  memset(vecs, 0, num_vec * dimension * sizeof(int64_t));
  vecs[dimension*0+0] = 1; vecs[dimension*0+1] = 0; vecs[dimension*0+2] = 1; // 2NN
  vecs[dimension*1+0] = -1; vecs[dimension*1+1] = 0; vecs[dimension*1+2] = 1; // 2NN
  vecs[dimension*2+0] = 0; vecs[dimension*2+1] = 1; vecs[dimension*2+2] = 1; // 2NN
  vecs[dimension*3+0] = 0; vecs[dimension*3+1] = -1; vecs[dimension*3+2] = 1; // 2NN
  vecs[dimension*4+0] = 1; vecs[dimension*4+1] = 1; vecs[dimension*4+2] = 0; // 2NN
  vecs[dimension*5+0] = 1; vecs[dimension*5+1] = -1; vecs[dimension*5+2] = 0; // 2NN
  vecs[dimension*6+0] = 1; // NN
  vecs[dimension*7+1] = 1; // NN
  vecs[dimension*8+2] = 1; // NN

  g = get_lattice_from_nd_simple_cubic_and_vectors(lsize, dimension, vecs, num_vec, static_center, periodic, get_size, list_edges);

  avg = 0;
  for(int i=0; i < n_avg; i++) {
    int64_t* percolated = percolate_site(&g, 1.0, &idxLambda); // site-percolation given fixed bond probability (Newman-Ziff method)
    free(percolated);
    avg += (double)idxLambda/g.nnode;
  }
  printf("%f ", avg/n_avg);
  printf("%f ", avg/n_avg/0.1372);

  avg = 0;
  for(int i=0; i < n_avg; i++) {
    int64_t* percolated = percolate_bond(&g, 1.0, &idxLambda); // bond-percolation given fixed site probability (Newman-Ziff method)
    free(percolated);
    avg += (double)idxLambda/g.num_edges;
  }
  printf("%f ", avg/n_avg);
  printf("%f ", avg/n_avg/0.0752326);
  printf("NN + 2NN (site/bond)\n");

  free_graph(&g);
  free(vecs);

  /* 2NN + 4NN */
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

  avg = 0;
  for(int i=0; i < n_avg; i++) {
    int64_t* percolated = percolate_site(&g, 1.0, &idxLambda); // site-percolation given fixed bond probability (Newman-Ziff method)
    free(percolated);
    avg += (double)idxLambda/g.nnode;
  }
  printf("%f ", avg/n_avg);
  printf("%f ", avg/n_avg/0.15950);

  avg = 0;
  for(int i=0; i < n_avg; i++) {
    int64_t* percolated = percolate_bond(&g, 1.0, &idxLambda); // bond-percolation given fixed site probability (Newman-Ziff method)
    free(percolated);
    avg += (double)idxLambda/g.num_edges;
  }
  printf("%f ", avg/n_avg);
  printf("%f ", avg/n_avg/0.0751589);
  printf("2NN + 4NN (site/bond)\n");

  free_graph(&g);
  free(vecs);

  /* 2NN + 3NN */
  num_vec = 10;
  vecs = malloc(num_vec * dimension * sizeof(int64_t));
  memset(vecs, 0, num_vec * dimension * sizeof(int64_t));
  vecs[dimension*0+0] = 1; vecs[dimension*0+1] = 0; vecs[dimension*0+2] = 1; // 2NN
  vecs[dimension*1+0] = -1; vecs[dimension*1+1] = 0; vecs[dimension*1+2] = 1; // 2NN
  vecs[dimension*2+0] = 0; vecs[dimension*2+1] = 1; vecs[dimension*2+2] = 1; // 2NN
  vecs[dimension*3+0] = 0; vecs[dimension*3+1] = -1; vecs[dimension*3+2] = 1; // 2NN
  vecs[dimension*4+0] = 1; vecs[dimension*4+1] = 1; vecs[dimension*4+2] = 0; // 2NN
  vecs[dimension*5+0] = 1; vecs[dimension*5+1] = -1; vecs[dimension*5+2] = 0; // 2NN
  vecs[dimension*6+0] = 1; vecs[dimension*6+1] = 1; vecs[dimension*6+2] = 1; // 3NN
  vecs[dimension*7+0] = 1; vecs[dimension*7+1] = -1; vecs[dimension*7+2] = 1; // 3NN
  vecs[dimension*8+0] = -1; vecs[dimension*8+1] = 1; vecs[dimension*8+2] = 1; // 3NN
  vecs[dimension*9+0] = -1; vecs[dimension*9+1] = -1; vecs[dimension*9+2] = 1; // 3NN

  g = get_lattice_from_nd_simple_cubic_and_vectors(lsize, dimension, vecs, num_vec, static_center, periodic, get_size, list_edges);

  avg = 0;
  for(int i=0; i < n_avg; i++) {
    int64_t* percolated = percolate_site(&g, 1.0, &idxLambda); // site-percolation given fixed bond probability (Newman-Ziff method)
    free(percolated);
    avg += (double)idxLambda/g.nnode;
  }
  printf("%f ", avg/n_avg);
  printf("%f ", avg/n_avg/0.1036);

  avg = 0;
  for(int i=0; i < n_avg; i++) {
    int64_t* percolated = percolate_bond(&g, 1.0, &idxLambda); // bond-percolation given fixed site probability (Newman-Ziff method)
    free(percolated);
    avg += (double)idxLambda/g.num_edges;
  }
  printf("%f ", avg/n_avg);
  printf("%f ", avg/n_avg/0.0629283);
  printf("2NN + 3NN (site/bond)\n");

  free_graph(&g);
  free(vecs);

  /* NN + 3NN + 4NN */
  num_vec = 10;
  vecs = malloc(num_vec * dimension * sizeof(int64_t));
  memset(vecs, 0, num_vec * dimension * sizeof(int64_t));
  vecs[dimension*0+0] = 1; // NN
  vecs[dimension*1+1] = 1; // NN
  vecs[dimension*2+2] = 1; // NN
  vecs[dimension*3+0] = 1; vecs[dimension*3+1] = 1; vecs[dimension*3+2] = 1; // 3NN
  vecs[dimension*4+0] = 1; vecs[dimension*4+1] = -1; vecs[dimension*4+2] = 1; // 3NN
  vecs[dimension*5+0] = -1; vecs[dimension*5+1] = 1; vecs[dimension*5+2] = 1; // 3NN
  vecs[dimension*6+0] = -1; vecs[dimension*6+1] = -1; vecs[dimension*6+2] = 1; // 3NN
  vecs[dimension*7+0] = 2; // 4NN
  vecs[dimension*8+1] = 2; // 4NN
  vecs[dimension*9+2] = 2; // 4NN

  g = get_lattice_from_nd_simple_cubic_and_vectors(lsize, dimension, vecs, num_vec, static_center, periodic, get_size, list_edges);

  avg = 0;
  for(int i=0; i < n_avg; i++) {
    int64_t* percolated = percolate_site(&g, 1.0, &idxLambda); // site-percolation given fixed bond probability (Newman-Ziff method)
    free(percolated);
    avg += (double)idxLambda/g.nnode;
  }
  printf("%f ", avg/n_avg);
  printf("%f ", avg/n_avg/0.11920);

  avg = 0;
  for(int i=0; i < n_avg; i++) {
    int64_t* percolated = percolate_bond(&g, 1.0, &idxLambda); // bond-percolation given fixed site probability (Newman-Ziff method)
    free(percolated);
    avg += (double)idxLambda/g.num_edges;
  }
  printf("%f ", avg/n_avg);
  printf("%f ", avg/n_avg/0.0624379);
  printf("NN + 3NN + 4NN (site/bond)\n");

  free_graph(&g);
  free(vecs);

  /* NN + 2NN + 4NN */
  num_vec = 12;
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
  vecs[dimension*9+0] = 1; // NN
  vecs[dimension*10+1] = 1; // NN
  vecs[dimension*11+2] = 1; // NN

  g = get_lattice_from_nd_simple_cubic_and_vectors(lsize, dimension, vecs, num_vec, static_center, periodic, get_size, list_edges);

  avg = 0;
  for(int i=0; i < n_avg; i++) {
    int64_t* percolated = percolate_site(&g, 1.0, &idxLambda); // site-percolation given fixed bond probability (Newman-Ziff method)
    free(percolated);
    avg += (double)idxLambda/g.nnode;
  }
  printf("%f ", avg/n_avg);
  printf("%f ", avg/n_avg/0.11440);

  avg = 0;
  for(int i=0; i < n_avg; i++) {
    int64_t* percolated = percolate_bond(&g, 1.0, &idxLambda); // bond-percolation given fixed site probability (Newman-Ziff method)
    free(percolated);
    avg += (double)idxLambda/g.num_edges;
  }
  printf("%f ", avg/n_avg);
  printf("%f ", avg/n_avg/0.0533056);
  printf("NN + 2NN + 4NN (site/bond)\n");

  free_graph(&g);
  free(vecs);

  /* NN + 2NN + 3NN */
  num_vec = 13;
  vecs = malloc(num_vec * dimension * sizeof(int64_t));
  memset(vecs, 0, num_vec * dimension * sizeof(int64_t));
  vecs[dimension*0+0] = 1; vecs[dimension*0+1] = 0; vecs[dimension*0+2] = 1; // 2NN
  vecs[dimension*1+0] = -1; vecs[dimension*1+1] = 0; vecs[dimension*1+2] = 1; // 2NN
  vecs[dimension*2+0] = 0; vecs[dimension*2+1] = 1; vecs[dimension*2+2] = 1; // 2NN
  vecs[dimension*3+0] = 0; vecs[dimension*3+1] = -1; vecs[dimension*3+2] = 1; // 2NN
  vecs[dimension*4+0] = 1; vecs[dimension*4+1] = 1; vecs[dimension*4+2] = 0; // 2NN
  vecs[dimension*5+0] = 1; vecs[dimension*5+1] = -1; vecs[dimension*5+2] = 0; // 2NN
  vecs[dimension*6+0] = 1; vecs[dimension*6+1] = 1; vecs[dimension*6+2] = 1; // 3NN
  vecs[dimension*7+0] = 1; vecs[dimension*7+1] = -1; vecs[dimension*7+2] = 1; // 3NN
  vecs[dimension*8+0] = -1; vecs[dimension*8+1] = 1; vecs[dimension*8+2] = 1; // 3NN
  vecs[dimension*9+0] = -1; vecs[dimension*9+1] = -1; vecs[dimension*9+2] = 1; // 3NN
  vecs[dimension*10+0] = 1; // NN
  vecs[dimension*11+1] = 1; // NN
  vecs[dimension*12+2] = 1; // NN

  g = get_lattice_from_nd_simple_cubic_and_vectors(lsize, dimension, vecs, num_vec, static_center, periodic, get_size, list_edges);

  avg = 0;
  for(int i=0; i < n_avg; i++) {
    int64_t* percolated = percolate_site(&g, 1.0, &idxLambda); // site-percolation given fixed bond probability (Newman-Ziff method)
    free(percolated);
    avg += (double)idxLambda/g.nnode;
  }
  printf("%f ", avg/n_avg);
  printf("%f ", avg/n_avg/0.0976);

  avg = 0;
  for(int i=0; i < n_avg; i++) {
    int64_t* percolated = percolate_bond(&g, 1.0, &idxLambda); // bond-percolation given fixed site probability (Newman-Ziff method)
    free(percolated);
    avg += (double)idxLambda/g.num_edges;
  }
  printf("%f ", avg/n_avg);
  printf("%f ", avg/n_avg/0.0497080);
  printf("NN + 2NN + 3NN (site/bond)\n");

  free_graph(&g);
  free(vecs);

  /* 2NN + 3NN + 4NN */
  num_vec = 13;
  vecs = malloc(num_vec * dimension * sizeof(int64_t));
  memset(vecs, 0, num_vec * dimension * sizeof(int64_t));
  vecs[dimension*0+0] = 1; vecs[dimension*0+1] = 0; vecs[dimension*0+2] = 1; // 2NN
  vecs[dimension*1+0] = -1; vecs[dimension*1+1] = 0; vecs[dimension*1+2] = 1; // 2NN
  vecs[dimension*2+0] = 0; vecs[dimension*2+1] = 1; vecs[dimension*2+2] = 1; // 2NN
  vecs[dimension*3+0] = 0; vecs[dimension*3+1] = -1; vecs[dimension*3+2] = 1; // 2NN
  vecs[dimension*4+0] = 1; vecs[dimension*4+1] = 1; vecs[dimension*4+2] = 0; // 2NN
  vecs[dimension*5+0] = 1; vecs[dimension*5+1] = -1; vecs[dimension*5+2] = 0; // 2NN
  vecs[dimension*6+0] = 1; vecs[dimension*6+1] = 1; vecs[dimension*6+2] = 1; // 3NN
  vecs[dimension*7+0] = 1; vecs[dimension*7+1] = -1; vecs[dimension*7+2] = 1; // 3NN
  vecs[dimension*8+0] = -1; vecs[dimension*8+1] = 1; vecs[dimension*8+2] = 1; // 3NN
  vecs[dimension*9+0] = -1; vecs[dimension*9+1] = -1; vecs[dimension*9+2] = 1; // 3NN
  vecs[dimension*10+0] = 2; // 4NN
  vecs[dimension*11+1] = 2; // 4NN
  vecs[dimension*12+2] = 2; // 4NN

  g = get_lattice_from_nd_simple_cubic_and_vectors(lsize, dimension, vecs, num_vec, static_center, periodic, get_size, list_edges);

  avg = 0;
  for(int i=0; i < n_avg; i++) {
    int64_t* percolated = percolate_site(&g, 1.0, &idxLambda); // site-percolation given fixed bond probability (Newman-Ziff method)
    free(percolated);
    avg += (double)idxLambda/g.nnode;
  }
  printf("%f ", avg/n_avg);
  printf("%f ", avg/n_avg/0.11330);

  avg = 0;
  for(int i=0; i < n_avg; i++) {
    int64_t* percolated = percolate_bond(&g, 1.0, &idxLambda); // bond-percolation given fixed site probability (Newman-Ziff method)
    free(percolated);
    avg += (double)idxLambda/g.num_edges;
  }
  printf("%f ", avg/n_avg);
  printf("%f ", avg/n_avg/0.0474609);
  printf("2NN + 3NN + 4NN (site/bond)\n");

  free_graph(&g);
  free(vecs);

  /* NN + 2NN + 3NN + 4NN */
  num_vec = 16;
  vecs = malloc(num_vec * dimension * sizeof(int64_t));
  memset(vecs, 0, num_vec * dimension * sizeof(int64_t));
  vecs[dimension*0+0] = 1; vecs[dimension*0+1] = 0; vecs[dimension*0+2] = 1; // 2NN
  vecs[dimension*1+0] = -1; vecs[dimension*1+1] = 0; vecs[dimension*1+2] = 1; // 2NN
  vecs[dimension*2+0] = 0; vecs[dimension*2+1] = 1; vecs[dimension*2+2] = 1; // 2NN
  vecs[dimension*3+0] = 0; vecs[dimension*3+1] = -1; vecs[dimension*3+2] = 1; // 2NN
  vecs[dimension*4+0] = 1; vecs[dimension*4+1] = 1; vecs[dimension*4+2] = 0; // 2NN
  vecs[dimension*5+0] = 1; vecs[dimension*5+1] = -1; vecs[dimension*5+2] = 0; // 2NN
  vecs[dimension*6+0] = 1; vecs[dimension*6+1] = 1; vecs[dimension*6+2] = 1; // 3NN
  vecs[dimension*7+0] = 1; vecs[dimension*7+1] = -1; vecs[dimension*7+2] = 1; // 3NN
  vecs[dimension*8+0] = -1; vecs[dimension*8+1] = 1; vecs[dimension*8+2] = 1; // 3NN
  vecs[dimension*9+0] = -1; vecs[dimension*9+1] = -1; vecs[dimension*9+2] = 1; // 3NN
  vecs[dimension*10+0] = 2; // 4NN
  vecs[dimension*11+1] = 2; // 4NN
  vecs[dimension*12+2] = 2; // 4NN
  vecs[dimension*13+0] = 1; // NN
  vecs[dimension*14+1] = 1; // NN
  vecs[dimension*15+2] = 1; // NN

  g = get_lattice_from_nd_simple_cubic_and_vectors(lsize, dimension, vecs, num_vec, static_center, periodic, get_size, list_edges);

  avg = 0;
  for(int i=0; i < n_avg; i++) {
    int64_t* percolated = percolate_site(&g, 1.0, &idxLambda); // site-percolation given fixed bond probability (Newman-Ziff method)
    free(percolated);
    avg += (double)idxLambda/g.nnode;
  }
  printf("%f ", avg/n_avg);
  printf("%f ", avg/n_avg/0.10000);

  avg = 0;
  for(int i=0; i < n_avg; i++) {
    int64_t* percolated = percolate_bond(&g, 1.0, &idxLambda); // bond-percolation given fixed site probability (Newman-Ziff method)
    free(percolated);
    avg += (double)idxLambda/g.num_edges;
  }
  printf("%f ", avg/n_avg);
  printf("%f ", avg/n_avg/0.0392312);
  printf("NN + 2NN + 3NN + 4NN (site/bond)\n");

  free_graph(&g);
  free(vecs);

  /* 4NN (same as NN / simple cubic lattice), consistency check, lsize must be odd as the graph consist of several mutually unconnected partitions */
  num_vec = 3;
  vecs = malloc(num_vec * dimension * sizeof(int64_t));
  memset(vecs, 0, num_vec * dimension * sizeof(int64_t));
  vecs[dimension*0+0] = 2; // 4NN
  vecs[dimension*1+1] = 2; // 4NN
  vecs[dimension*2+2] = 2; // 4NN

  g = get_lattice_from_nd_simple_cubic_and_vectors(lsize, dimension, vecs, num_vec, static_center, periodic, get_size, list_edges);

  avg = 0;
  for(int i=0; i < n_avg; i++) {
    int64_t* percolated = percolate_site(&g, 1.0, &idxLambda); // site-percolation given fixed bond probability (Newman-Ziff method)
    free(percolated);
    avg += (double)idxLambda/g.nnode;
  }
  printf("%f ", avg/n_avg);
  printf("%f ", avg/n_avg/0.3116);

  avg = 0;
  for(int i=0; i < n_avg; i++) {
    int64_t* percolated = percolate_bond(&g, 1.0, &idxLambda); // bond-percolation given fixed site probability (Newman-Ziff method)
    free(percolated);
    avg += (double)idxLambda/g.num_edges;
  }
  printf("%f ", avg/n_avg);
  printf("%f ", avg/n_avg/0.2488);
  printf("4NN (site/bond)\n");

  free_graph(&g);
  free(vecs);

  return 0;
}
