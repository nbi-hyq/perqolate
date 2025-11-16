#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"
#include "../inc/percolation_main.h"

#define MAX_LINE 1024

static int scan_edges(char* line, uint8_t* arry){
  char* p = line;
  int i = 0;
  while(p < line+MAX_LINE) {
    char* end;
    long int value = strtol(p, &end, 10);
    if(value == 0L && end == p){
      break;
    } else {
      arry[i] = (uint8_t)value;
      i++;
    }
    p = end;
  }
  return i;
}

static int scan_edge_vecs(char* line, int8_t* arry){
  char* p = line;
  int i = 0;
  while(p < line+MAX_LINE) {
    char* end;
    long int value = strtol(p, &end, 10);
    if(value == 0L && end == p){
      break;
    } else {
      arry[i] = (int8_t)value;
      i++;
    }
    p = end;
  }
  return i;
}

/* construct periodic graph from unit-graph */
int main(int argc, char **argv){
  /* set simulation parameters */
  srand(2572457);
  int64_t lsize = 200; // size in each dimension
  int n_avg = 200; // number of repetitions for averaging
  bool periodic = false; // no periodic boundaries
  bool get_size = false; // check percolation, not determine size of largest connected component
  bool static_center = false;
  uint8_t dimension = 2; // dimension of lattice
  bool list_edges = true; // must be true for bond-percolation

  /* check number of args, count lines */
  if(argc != 2) return -1;
  FILE* fp = fopen(argv[1], "r");
  if (fp == NULL){
      printf("Could not open file %s", argv[1]);
      return -2;
  }
  int line_count = 0;
  char c;
  for (c = getc(fp); c != EOF; c = getc(fp))
      if (c == '\n') // count new lines
          line_count = line_count + 1;
  rewind(fp);

  /* simulate all lattices */
  char* line = malloc(sizeof(char) * MAX_LINE);
  uint8_t* blk_edges = malloc(sizeof(uint8_t) * MAX_LINE);
  int8_t* blk_vec = malloc(sizeof(uint8_t) * MAX_LINE);
  for (int i_lattice=0; i_lattice < line_count/2; i_lattice++){ 
    /* load and create the unit graph */ 
    fgets(line, MAX_LINE, fp); // 1st line is the egdes of the unit graph
    uint8_t nedge = scan_edges(line, blk_edges) / 2;
    fgets(line, MAX_LINE, fp); // 2nd line is the edge vectors of the unit graph
    scan_edge_vecs(line, blk_vec);

    UnitGraph unt = new_unit_graph(blk_edges, blk_vec, nedge, dimension);

    /* create periodic graph */
    Graph g = get_lattice_from_unit_graph(&unt, lsize, dimension, static_center, periodic, get_size, list_edges);

    /* run simulation */
    int64_t idxLambda;
    double avg = 0;
    for(int i=0; i < n_avg; i++) {
      int64_t* percolated = percolate_site(&g, 1.0, &idxLambda); // site-percolation given fixed bond probability (Newman-Ziff method)
      free(percolated);
      avg += (double)idxLambda/g.nnode;
    }
    printf("%i %f ", i_lattice, avg/n_avg);

    avg = 0;
    for(int i=0; i < n_avg; i++) {
      int64_t* percolated = percolate_bond(&g, 1.0, &idxLambda); // bond-percolation given fixed site probability (Newman-Ziff method)
      free(percolated);
      avg += (double)idxLambda/g.num_edges;
    }
    printf("%f\n", avg/n_avg);

    free_graph(&g);
    free_unit_graph(&unt);
  }

  free(blk_edges);
  free(blk_vec);
  fclose(fp);
  free(line);
  return 0;
}
