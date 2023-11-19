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

  /* check permutation */
  int num_rep = 5000;
  int64_t len = 100;
  int64_t* sum = malloc(len * sizeof(int64_t));
  memset(sum, 0, len * sizeof(int64_t));
  for (int r=0; r<num_rep; r++){
    int64_t* arry = permutation(len);
    for (int64_t i=0; i<len; i++) sum[i] += arry[i];
    free(arry);
  }

  for (int64_t i=0; i<len; i++){
    double h = (double)sum[i] / ((double)(len - 1)/2.0) / num_rep;
    printf("%f, ", h);
    if (h < 1 - 4/sqrt(num_rep) || h > 1 + 4/sqrt(num_rep)) rt = 1;
  }
  free(sum);

  /* check random number uniformity */
  uint64_t a_rg[5] = {4, 16, 21, 63, 64};
  for (int j=0; j<5; j++){
    int* frequ = malloc(a_rg[j] * sizeof(int64_t));
    memset(frequ, 0, a_rg[j] * sizeof(int64_t));
    num_rep = 100000;
    for(int r=0; r<num_rep; r++){
      uint64_t k = rand_range(a_rg[j]);
      frequ[k] += 1;
    }
    for (uint64_t i=0; i < a_rg[j]; i++){
      int h = num_rep / a_rg[j];
      if (frequ[i] < h - 4*sqrt(h) || frequ[i] > h + 4*sqrt(h)) rt = 1;
      printf("%f %i %f\n ", h - 5*sqrt(h), frequ[i], h + 5*sqrt(h));
    }
    free(frequ);
  }

  return rt;
}
