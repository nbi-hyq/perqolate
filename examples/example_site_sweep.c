#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "inc/graph_type.h"
#include "inc/graph_construction.h"
#include "inc/percolation_main.h"

/* site-percolation simulation given fixed bond probability (using Newman-Ziff method) */
int main(){
  srand(time(NULL));
  bool periodic = false; // no periodic boundaries
  bool get_size = false; // percolation
  int sweep_len = 100; // number of different probabilities
  int n_rep = 3; // number of repetitions for averaging
  float* sweep_prob; // site probability (tuned)
  sweep_prob = malloc(sweep_len * sizeof(float));
  for(int i=0;i<sweep_len;i++) sweep_prob[i] = 0.5/sweep_len + 1.0*i/sweep_len;

  for(int64_t lsize=100;lsize<401;lsize=lsize+100){
    Graph g = get_2d_square(lsize, periodic, get_size, false);
    int64_t* percolated_sum  = NULL;
    for(int rep=0;rep<n_rep;rep++){
      int64_t idxLambda; // number of bonds where percolation appears
      int64_t* percolated = percolate_site(&g, 0.726195, &idxLambda); // site-percolation given fixed bond probability (Newman-Ziff method)
      printf("%f\n", ((double)idxLambda - 0.5)/ (double)g.nnode);
      if(rep==0){
        percolated_sum = percolated;
      } else{
        for(int64_t i=0; i<g.nnode+1;i++) percolated_sum[i] += percolated[i];
        free(percolated);
      }
    }
    float* expectation_value = get_expectation_value(sweep_prob, sweep_len, percolated_sum, g.nnode, 0, n_rep, true, false); // compute a site-percolation curve
    free(expectation_value);
    free(percolated_sum);
    free_graph(&g);
  }
  free(sweep_prob);
  return 0;
}
