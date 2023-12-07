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
    double h = (double)sum[i] / num_rep;
    double m = (double)(len - 1)/2.0;
    if (h < m - 2.0*m/sqrt(12)/sqrt(num_rep)*5 || h > m + 2.0*m/sqrt(12)/sqrt(num_rep)*5) rt = 1; // h is uniformly distributed for num_rep==1
  }
  free(sum);

  printf("%i", rt);
  return rt;
}
