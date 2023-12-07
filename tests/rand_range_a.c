#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"
#include "../inc/percolation_main.h"

int main(){
  srand(73628);
  int rt = 0; // return

  /* check random number uniformity */
  uint64_t a_rg[5] = {4, 16, 21, 63, 64};
  for (int j=0; j<5; j++){
    int* frequ = malloc(a_rg[j] * sizeof(int64_t));
    memset(frequ, 0, a_rg[j] * sizeof(int64_t));
    int num_rep = 100000;
    for(int r=0; r<num_rep; r++){
      uint64_t k = rand_range(a_rg[j]);
      frequ[k] += 1;
    }
    for (uint64_t i=0; i < a_rg[j]; i++){
      int h = num_rep / a_rg[j];
      if (frequ[i] < h - 5*sqrt(h) || frequ[i] > h + 5*sqrt(h)){ // frequ[i] is Poisson distributed
        rt = 1;
        printf("%f %i %f\n ", h - 5*sqrt(h), frequ[i], h + 5*sqrt(h));
      }
    }
    free(frequ);
  }

  printf("%i", rt);
  return rt;
}
