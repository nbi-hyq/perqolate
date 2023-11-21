#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"
#include "../inc/percolation_main.h"
#include "../inc/global.h"

int main(){
  srand(73628);
  int rt = 0; // return
  int n_max = CSTD_RAND_MAX;

  /* simple check for random number uniformity over maximum range */
  for(int r=0;r<200;r++){
    double mean = 0;
    int cnt_low = 0;
    int num_rep = 100000;
    int* num = malloc(num_rep * sizeof(int));
    for (int j=0; j<num_rep; j++){
      num[j] = rand() & n_max;
      mean += (double)num[j] / num_rep;
      if (num[j] < n_max/2) cnt_low++;
    }
    if (cnt_low > num_rep/2 + 5*sqrt(num_rep)/2 || cnt_low < num_rep/2 - 5*sqrt(num_rep)/2){
      rt = 1;
      printf("%f %i %f\n", num_rep/2 - 5*sqrt(num_rep)/2, cnt_low, num_rep/2 + 5*sqrt(num_rep)/2);
    }
    if (mean > n_max/2.0 + n_max/sqrt(12)/sqrt(num_rep)*5 || mean < n_max/2.0 - n_max/sqrt(12)/sqrt(num_rep)*5){ // mean is uniformly distributed for num_rep==1
      rt = 1;
      printf("%f\n %f\n %f\n", n_max/2.0 - n_max/sqrt(12)/sqrt(num_rep)*5, mean, n_max/2.0 + n_max/sqrt(12)/sqrt(num_rep)*5);
    }
    free(num);
  }

  printf("%i", rt);
  return rt;
}
