#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "inc/graph_type.h"
#include "inc/graph_construction.h"
#include "inc/percolation_main.h"

/* construct graph by fusing corresponding GHZ-states with photon loss (using Newman-Ziff method)
   add edges step by step and check if the loss threshold gets lower */
int main(){
  srand(time(NULL));
  int64_t lsize = 70; // size in each dimension
  uint8_t dimension = 3; // lattice dimension
  bool static_center = true; // central node cannot be lost
  bool periodic = false; // no periodic boundaries
  bool get_size = false; // look for percolation (cluster spanning), not largest component size

  uint8_t max_step = 1; // maximum number of steps per dimension
  Graph g;
  
  int64_t num_vec_max = 1;
  for(uint8_t j=0;j<dimension;j++) num_vec_max = num_vec_max * (2*(int64_t)max_step + 1);
  num_vec_max -= 1; // -1 because 0-vector is not considered
  num_vec_max /= 2; // if two vectors just differ by sign, just look at one of them

  int repetitions_per_run = 2*num_vec_max; // make larger, if routine stops before reaching local optium
  int number_of_runs = 3; // overall number of repetitions
  double lambda_min = 1.0; // result of best percolation simulation in a run
  double lambda_min_best_run = 1.0; // minimum for all runs
  int64_t cnt_fail = 0; // count how often reducing the number of vectors made it worse

  int64_t* vecs; // vectors specifying the lattice
  int64_t num_vec; // number of vectors for the lattice
  int64_t* vecs_best_run = malloc(sizeof(int64_t) * (int64_t)dimension * num_vec_max); // vectors of best run
  int64_t num_vec_best_run = 0; // number of vectors for best run (initial value has no meaning here)

  /* get best lattice by greedy optimization */
  for(int run=0; run<number_of_runs; run++){
    vecs = get_non_equivalent_vecs(dimension, max_step);
    num_vec = 1; // start with just a single vector
    for(int rep=0; rep<repetitions_per_run; rep++){
      for(int i=0; i<num_vec; i++){
        for(uint8_t j=0; j<dimension; j++) printf("%li ", vecs[i*(int64_t)dimension + j]);
        printf("%s\n", " ");
      }

      g = get_lattice_from_nd_simple_cubic_and_vectors(lsize, dimension, vecs, num_vec, static_center, periodic, get_size, false);

      int64_t numQbts; // 2 fusion qubits per edge -- overall number of qubits is more than number of nodes
      int64_t numStatic; // counts static nodes (is set 0 when static_center==false)
      int64_t idxLambda;
      int64_t* percolated = fusion_and_loss_nz(&g, 0.5, static_center, &numQbts, &numStatic, &idxLambda);
      free_graph(&g);
      double lambda = ((double)idxLambda -(double)numStatic - 0.5) / ((double)numQbts -(double)numStatic); // estimation of percolation threshold
      printf("%s%f\n", "----> lambda = ", lambda);

      if(rep == 0 || lambda < lambda_min){ // then: go on with one vector more
        lambda_min = lambda;
        num_vec += 1;
        cnt_fail = 0;
      }else { // then: add a different vector and remove the most recently added one
        cnt_fail++;
        if(num_vec - 1 + cnt_fail > num_vec_max - 1){
          free(percolated);
          break; //adding a vector always makes the lattice worse
        }
        int64_t temp[dimension];
        memcpy(temp, vecs + (num_vec - 1 + cnt_fail) * (int64_t)dimension, sizeof(int64_t) * (int64_t)dimension);
        memcpy(vecs + (num_vec - 1 + cnt_fail) * (int64_t)dimension, vecs + (num_vec - 1) * (int64_t)dimension, sizeof(int64_t) * (int64_t)dimension);
        memcpy(vecs + (num_vec - 1) * (int64_t)dimension, temp, sizeof(int64_t) * (int64_t)dimension);
      }
      free(percolated);
    }
    num_vec--; // keep one vector less for final result
    if(lambda_min < lambda_min_best_run){ // save the results
      num_vec_best_run = num_vec;
      memcpy(vecs_best_run, vecs, sizeof(int64_t) * (int64_t)dimension * num_vec_max);
      lambda_min_best_run = lambda_min;
    }
    printf("----> lambda_min: %f\n", lambda_min);
    printf("-------- next run --------\n");
  }
  free(vecs);

  /* print the final result */
  printf("final result: -------\n");
  for(int i=0; i<num_vec_best_run; i++){
    for(uint8_t j=0; j<dimension; j++) printf("%li ", vecs_best_run[i*(int64_t)dimension + j]);
    printf("%s\n", " ");
  }
  printf("----> lambda_min_best_run: %f\n", lambda_min_best_run);

  /* range for simulations of best lattice found */
  int n_rep = 3; // number of repetitions
  int sweep_len = 400; // number of different probabilities
  float* sweep_prob; // 1 - p_loss (tuned)
  sweep_prob = malloc(sweep_len * sizeof(float));
  for(int i=0;i<sweep_len;i++) sweep_prob[i] = 0.90 + 0.1*(0.5/sweep_len + 1.0*i/sweep_len);

  /* simulate percolation with best lattice found */
  printf("%s\n", "---percolation simulation----");
  for(lsize=20;lsize<81;lsize=lsize+20){
    g = get_lattice_from_nd_simple_cubic_and_vectors(lsize, dimension, vecs_best_run, num_vec_best_run, static_center, periodic, get_size, false);
    int64_t* percolated_sum = NULL;
    int64_t numQbts; // 2 fusion qubits per edge -- overall number of qubits is more than number of nodes
    int64_t numStatic; // counts static nodes (is set 0 when static_center==false)
    for(int rep=0;rep<n_rep;rep++){
      int64_t idxLambda;
      int64_t* percolated = fusion_and_loss_nz(&g, 0.5, static_center, &numQbts, &numStatic, &idxLambda);
      printf("%f\n", ((double)idxLambda -(double)numStatic - 0.5) / ((double)numQbts -(double)numStatic));
      if(rep==0){
        percolated_sum = percolated;
      } else{
        for(int64_t i=0; i<numQbts+1;i++) percolated_sum[i] += percolated[i];
        free(percolated);
      }
    }
    float* expectation_value = get_expectation_value(sweep_prob, sweep_len, percolated_sum, numQbts, numStatic, n_rep, true, false); // compute fusion/loss-percolation curve
    free(percolated_sum);
    free(expectation_value);
    free_graph(&g);
  }

  /* simulate largest component size with best lattice found */
  printf("%s\n", "---largest component simulation----");
  periodic = true; // periodic boundaries for this simulation
  get_size = true; // look for largest component size
  for(lsize=20;lsize<81;lsize=lsize+20){
    g = get_lattice_from_nd_simple_cubic_and_vectors(lsize, dimension, vecs_best_run, num_vec_best_run, static_center, periodic, get_size, false);
    int64_t* largest_component_sum = NULL;
    int64_t numQbts; // 2 fusion qubits per edge -- overall number of qubits is more than number of nodes
    int64_t numStatic; // counts static nodes (is set 0 when static_center==false)
    for(int rep=0;rep<n_rep;rep++){
      int64_t idxLambda;
      int64_t* largest_component = fusion_and_loss_nz(&g, 0.5, static_center, &numQbts, &numStatic, &idxLambda);
      if(rep==0){
        largest_component_sum = largest_component;
      } else{
        for(int64_t i=0; i<numQbts+1;i++) largest_component_sum[i] += largest_component[i];
        free(largest_component);
      }
    }
    float* expectation_value = get_expectation_value(sweep_prob, sweep_len, largest_component_sum, numQbts, numStatic, g.nnode * n_rep, true, false); // compute fusion/loss-percolation curve
    free(largest_component_sum);
    free(expectation_value);
    free_graph(&g);
  }

  free(sweep_prob);
  free(vecs_best_run);
  return 0;
}

