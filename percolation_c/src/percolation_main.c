#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include "../inc/global.h"
#include "../inc/graph_type.h"
#include "../inc/graph_construction.h"
#define CSTD_RAND_MAX (0x7FFF)

// find root and path compression
int64_t findroot(Graph* g, int64_t i){
  if (g->ptr[i]<0) return i;  // return index of root node
  return g->ptr[i] = findroot(g, g->ptr[i]);  // recursively go to root node, plus: do path-compression on the fly
}

/* merge two graph fragments in g->ptr representation, return new root node index */
int64_t merge_root(Graph* g, int64_t r1, int64_t r2, int64_t* big){
  if (g->ptr[r1] > g->ptr[r2]){
    if (g->ptr[r2] > START) g->ptr[r2] += g->ptr[r1]; // add size of smaller component to bigger one (unless one is START/STOP for indicating cluster spanning)
    g->ptr[r1] = r2; // attach component with root r1 to larger component with root r2
    r1 = r2; // update root
  } else {
    if (g->ptr[r1] > START) g->ptr[r1] += g->ptr[r2]; // add size of smaller component to bigger one (unless one is START/STOP for indicating cluster spanning)
    g->ptr[r2] = r1; // attach component with root r2 to larger component with root r1
  }
  if (g->get_size && -g->ptr[r1] > *big) *big = -g->ptr[r1]; // update largest component size
  return r1;
}

/* helper function: give every (undirected) edge a unique number */
int64_t get_edge_idx(Graph* g, int64_t s1, int64_t s2){
  uint8_t k;
  if (s1 < s2){
    for(k=0;k<g->len_nb[s1];k++) if (g->nn[s1*g->num_nb_max + k] == s2) break;
    return s1*g->num_nb_max + k;
  } else {
    for(k=0;k<g->len_nb[s2];k++) if (g->nn[s2*g->num_nb_max + k] == s1) break;
    return s2*g->num_nb_max + k;
  }
}

/* update g->ptr when a new node has been added */
static bool update_ptr_when_new_node(Graph* g, int64_t s1, int64_t* big, bool* bond){
  int64_t s2, r1, r2;
  bool percolated = false;
  r1 = s1; // start as own root
  for (uint8_t j=0;j<g->len_nb[s1];j++){ // go through all (possible) neighbors of s1
    s2 = g->nn[s1*g->num_nb_max + j];
    if (g->ptr[s2]>MEASUREZ && bond[get_edge_idx(g, s1, s2)]){ // do nothing if neighbor node is empty (LOST or MEASUREZ), or bond is missing
      r2 = findroot(g, s2); // find root of neighbor and do path compression
      if (r2 != r1) { // if connected to same root via another neighbor before, do nothing
        if(!g->get_size && ((g->ptr[r1] == START && g->ptr[r2] == STOP) || (g->ptr[r2] == START && g->ptr[r1] == STOP))){
          percolated = true;
          break; // percolation (no change when adding more nodes)
        }
        r1 = merge_root(g, r1, r2, big); // update largest component size, update root
      }
    }
  }
  return percolated;
}

/* update g->ptr when a new node has been added, use different data-structure to check that bond is present (compared to update_ptr_when_new_node) */
static bool update_ptr_when_new_node_b(Graph* g, int64_t s1, int64_t* big, uint8_t* arr_n_fusion){
  int64_t s2, r1, r2;
  bool percolated = false;
  r1 = s1; // start as own root
  for (uint8_t j=0;j<g->len_nb[s1];j++){ // go through all (possible) neighbors of s1
    s2 = g->nn[s1*g->num_nb_max + j];
    if (g->ptr[s2]>MEASUREZ && arr_n_fusion[s1*g->num_nb_max + j] < 255){ // do nothing if neighbor node is empty (LOST or MEASUREZ), or bond is missing (arr_n_fusion[] == 255)
      r2 = findroot(g, s2); // find root of neighbor and do path compression
      if (r2 != r1) { // if connected to same root via another neighbor before, do nothing
        if(!g->get_size && ((g->ptr[r1] == START && g->ptr[r2] == STOP) || (g->ptr[r2] == START && g->ptr[r1] == STOP))){
          percolated = true;
          break; // percolation (no change when adding more nodes)
        }
        r1 = merge_root(g, r1, r2, big); // update largest component size, update root
      }
    }
  }
  return percolated;
}

/* get random number with 64-bit precision (CSTD_RAND_MAX=0x7FFF: RAND_MAX guaranteed by C standard) */
static uint64_t rand64(void) {
  const uint64_t b0 = rand() & CSTD_RAND_MAX;
  const uint64_t b1 = rand() & CSTD_RAND_MAX;
  const uint64_t b2 = rand() & CSTD_RAND_MAX;
  const uint64_t b3 = rand() & CSTD_RAND_MAX;
  const uint64_t b4 = rand();
  const uint64_t k = (b0) | (b1 << 15) | (b2 << 30) | (b3 << 45) | (b4 << 60);
  return k;
}

/* get random number with more than 32-bit precision, uniform in range [0, n) */
uint64_t rand_range(uint64_t n) {
  uint64_t max = (1ul << 63) - (1ul << 63) % n - 1;
  uint64_t r = rand64();
  while (r > max) r = rand64();
  return r % n;
}

/* helper function: get array with shuffled indices */
int64_t* permutation(int64_t len){
  int64_t* array;
  array = malloc(len * sizeof(int64_t));
  int64_t i,j,temp;
  for (i=0;i<len;i++) array[i] = i;
  for (i=0;i<len;i++){
    j = i + rand_range(len-i);
    temp = array[i];
    array[i] = array[j];
    array[j] = temp;
  }
  return array;
}

/* helper function: shuffle an array, if(staticFirst) put all indices with arryStatic[i]==True in the beginning */
int64_t* permutation_static_first(int64_t len, bool staticFirst, bool* arryStatic, int64_t* numStatic){
  *numStatic = 0; // counter of static nodes
  if(staticFirst == false){
    return permutation(len);
  } else {
    int64_t* idxStatic = malloc(len * sizeof(int64_t));
    int64_t* idxNonStatic = malloc(len * sizeof(int64_t));
    int64_t numNonStatic = 0;
    for(int64_t i=0; i<len; i++){
      if(arryStatic[i]){
        idxStatic[(*numStatic)++] = i;
      } else {
        idxNonStatic[numNonStatic++] = i;
      }
    }
    if(numNonStatic > 0){ // 0: all nodes are set static -- nothing needs to be done
      int64_t* order = permutation(numNonStatic);
      for(int64_t k=0; k<numNonStatic; k++){
        idxStatic[(*numStatic)+k] = idxNonStatic[order[k]]; // attach shuffled nonStatic indices at end
      }
      free(order);
    }
    free(idxNonStatic);
    return idxStatic;
  }
}

/* classical: percolation simulation, sweep site occupation (pBond fixed) */
int64_t* percolate_site(Graph* g, float pBond, int64_t* idxLambda){ // TBD: support static nodes
  int64_t i,s1; // s1: start node index
  bool stop_loop;
  int64_t big = 0; // size of largest component
  bool* bond = malloc(g->nnode * (size_t)g->num_nb_max * sizeof(bool));
  for(i=0;i<g->nnode * g->num_nb_max;i++) bond[i] = (pBond > (float)rand()/(float)RAND_MAX); // randomly remove edges
  int64_t* order = permutation(g->nnode);  // random ordering of node indices
  int64_t* largest_component = malloc((g->nnode+1) * sizeof(int64_t)); // size of largest component (number of nodes is position of array)
  largest_component[0] = 0; // without nodes size is 0 and there is no percolation (also 0)
  for (i=0;i<g->nnode;i++) g->ptr[i] = LOST; // start without any nodes
  for (i=0;i<g->nnode;i++){
    s1 = order[i];
    if(g->get_size) g->ptr[s1] = -1; // isolated node, own root with size 1
    else if(g->start[s1]) g->ptr[s1] = START; // label new node as start node for percolation
    else if(g->stop[s1]) g->ptr[s1] = STOP; // label new node as stop node for percolation
    else g->ptr[s1] = -1; // isolated node, own root with size 1
    stop_loop = update_ptr_when_new_node(g, s1, &big, bond);
    if(stop_loop) break; // stop loop when percolation is achieved
    if(g->get_size) largest_component[i+1] = big;
    else largest_component[i+1] = 0; // no percolation yet
  }
  if(!g->get_size){
    for(int64_t k=i;k<g->nnode;k++) largest_component[k+1] = 1; // set to 1 where percolation is achieved
    *idxLambda = i+1;
  }
  free(bond);
  free(order);
  return largest_component;
}

// classical: percolation simulation, sweep bond occupation (pSite fixed), using Newman-Ziff algorithm
int64_t* percolate_bond(Graph* g, float pSite, int64_t* idxLambda){ // TBD: support static bonds, option to remove specific nodes
  int64_t i,s1,s2,r1,r2; // s1,s2: start node index; r1,r2: root node index
  int64_t big = 0; // size of largest component
  int64_t* order = permutation(g->num_edges);  // random ordering of edge indices
  int64_t* largest_component = malloc((g->num_edges+1) * sizeof(int64_t)); // size of largest component resp. 0/1 for percolation or not (number of edges is position of array)
  largest_component[0] = (g->get_size ? 1 : 0); // without edges size is 1 (isolated node), and no percolation (0)
  for(i=0;i<g->nnode;i++){
    if (pSite < (float)rand()/(float)RAND_MAX && g->static_node[i]==false) g->ptr[i] = LOST; // remove this node
    else if (!g->get_size && g->start[i]) g->ptr[i] = START; // start nodes for pecolation
    else if (!g->get_size && g->stop[i]) g->ptr[i] = STOP; // stop nodes for pecolation
    else g->ptr[i] = -1; // make node isolated (own root)
  }
  for (i=0;i<g->num_edges;i++){
    s1 = g->edge_node[order[i]]; // get node index 1 of the edge
    s2 = g->nn[s1*g->num_nb_max + g->edge_idx[order[i]]]; // get node index 2 of the edge
    if(g->ptr[s2]!=LOST) r2 = findroot(g, s2); //find root and do path compression
    if(g->ptr[s1]!=LOST) r1 = findroot(g, s1); //find root and do path compression
    if (g->ptr[s1]!=LOST && g->ptr[s2]!=LOST && r2!=r1){ // if both nodes of the edge exist and root is not the same
      if(!g->get_size && ((g->ptr[r1] == START && g->ptr[r2] == STOP) || (g->ptr[r2] == START && g->ptr[r1] == STOP))) break; // percolation (no change when adding more bonds)
      merge_root(g, r1, r2, &big); // update largest component size, update root (return not needed)
    }
    if(g->get_size) largest_component[i+1] = big; // update size of largest component
    else largest_component[i+1] = 0; // no percolation yet
  }
  if(!g->get_size){
    for(int64_t k=i;k<g->num_edges;k++) largest_component[k+1] = 1; // set to 1 where percolation is achieved
    *idxLambda = i+1;
  }
  free(order);
  return largest_component;
}

/* quantum: given lost qubit and fusions, determine which other qubits must be measured in Z */
static void crawl_nb_label_measure_z(Graph* g, int64_t node_idx){
  int64_t nb_idx;
  for(uint8_t i=0;i<g->len_nb[node_idx];i++){
    nb_idx = g->nn[node_idx*g->num_nb_max + i];
    if(g->ptr[nb_idx] != LOST && g->ptr[nb_idx] != MEASUREZ){
      g->ptr[nb_idx] = MEASUREZ;
      if(g->fusion_node[nb_idx] && g->fusion_success[nb_idx]){
        crawl_nb_label_measure_z(g, g->fusion_partner[nb_idx]); // crawl further starting from fusion parter of the neighbor
      }
    }
  }
}

/* crawl connected graph component starting from a given node, don't take paths with LOST, ALLPLYZ
   return: size of graph component including the start node, return 0/1 for percolation no/yes if g->get_size==false
   application1: fusion networks using g->fusion_node (resource states not connected in g->ptr)
   application2: loss on given cluster states (note: all g->fusion_node[i] must be 0) */
static int64_t crawl_graph_bfs(Graph* g, int64_t* node_queue, bool* visited, int64_t node_idx){
  int64_t cnt, node_queue_pos, node_queue_len, nb_idx;
  node_queue_pos = 0;
  cnt = 0;
  node_queue_len = 1;
  node_queue[0] = node_idx;
  while(node_queue_pos < node_queue_len){
    node_idx = node_queue[node_queue_pos];
    for(uint8_t i=0;i<g->len_nb[node_idx];i++){ // go through all neighbors
      nb_idx = g->nn[node_idx*g->num_nb_max + i];
      if(visited[nb_idx]==false && g->ptr[nb_idx]!=LOST && g->ptr[nb_idx]!=MEASUREZ){
        if(!g->get_size && g->stop[nb_idx] && !g->fusion_node[nb_idx]) return -1;  // percolation achieved
        node_queue[node_queue_len] = nb_idx;
        node_queue_len++;
        visited[nb_idx] = true;
      }
    }
    if(g->fusion_node[node_idx]){
      if(g->fusion_success[node_idx] && visited[g->fusion_partner[node_idx]]==false){  // no need to check if fusion_partner is LOST or APPLY (since fusion and lost nodes not called initially, if one fusion node is Z -- the other is surrounded by Z)
        node_queue[node_queue_len] = g->fusion_partner[node_idx]; // add the fusion partner to the queue
        node_queue_len++;
        visited[g->fusion_partner[node_idx]] = true;
      }
    } else {
      cnt++; // count only the nodes that are not used for fusion
    }
    node_queue_pos++;
  }
  return cnt;
}

/* traverse full graph state or fusion network breadth first
   return: largest graph component or information about cluster spanning
   application1: fusion networks using g->fusion_node (resource states not connected in g->ptr)
   application2: loss on given cluster states (note: all g->fusion_node[i] must be 0) */
int64_t bfs_full(Graph* g){
  int64_t big = 0;
  int64_t cnt = 0;
  int64_t* node_queue = malloc(g->nnode * sizeof(int64_t));
  bool* visited = malloc(g->nnode * sizeof(bool));
  memset(visited, 0, g->nnode * sizeof(bool));
  if(!g->get_size) {  // test for percolation from one side to the other
    for(int64_t i=0;i<g->nnode;i++){
      if(g->start[i] && g->ptr[i] != LOST && g->ptr[i] != MEASUREZ && visited[i]==false && g->fusion_node[i]==false){
        visited[i] = true;
        if(crawl_graph_bfs(g, node_queue, visited, i) == -1){
          big=1; // percolation achieved
          break; // no need to look at more nodes
        }
      }
    }
  } else {  // get size of largest connected component
    for(int64_t i=0;i<g->nnode;i++){
      if(g->ptr[i] != LOST && g->ptr[i] != MEASUREZ && visited[i]==false && g->fusion_node[i]==false){
        visited[i] = true;
        cnt = crawl_graph_bfs(g, node_queue, visited, i);
        if(cnt > big) big = cnt;
      }
    }
  }
  free(node_queue);
  free(visited);
  return big;
}

/* quantum: create large graph by fusions (inluding loss) and determine largest component (or percolation) */
int64_t fusion_and_loss_no_newman_ziff(Graph* g, float pFusion, float pNoLoss){
  int64_t i;
  // 1a) initalize fusion success and loss probabilistically
  for(i=0;i<g->nnode;i++) g->ptr[i] = -1; // initialize, number means nothing (but not lost) here
  for(i=0;i<g->nnode;i++){
    if(pNoLoss < (float)rand()/(float)RAND_MAX && g->static_node[i]==false){
      g->ptr[i] = LOST;
      if(g->fusion_node[i]){
        g->ptr[g->fusion_partner[i]] = LOST; // set the fusion partner lost (can also affect static_node), too (since both are looped, this can happen 2x)
      }
    }
  }
  // 1b) initialize probabilistic fusion success
  for(i=0;i<g->nnode;i++){
    if(g->fusion_node[i] && g->fusion_partner[i]>i){
      if(pFusion > (float)rand()/(float)RAND_MAX){
        g->fusion_success[i] = true;
        g->fusion_success[g->fusion_partner[i]] = true;
      } else {
        g->fusion_success[i] = false;
        g->fusion_success[g->fusion_partner[i]] = false;
      }
    }
  }

  // 2) crawl neighborhood of lost qubits and distribute MEASUREZ
  for(i=0;i<g->nnode;i++) if(g->ptr[i] == LOST) crawl_nb_label_measure_z(g, i);

  // 3) breadth-first traversal of the graph to get its largest connected component (or check for cluster spanning)
  int64_t big = bfs_full(g);
  return big;
}

/* quantum case only: apply loss to graph by picking random qubits for loss, no Newman-Ziff method (alternative: apply_loss_bfs)
   pSite: probability that photon survives
   TBD: missing bond option (if loss is applied after fusion failures) */
int64_t apply_loss(Graph* g, float pSite){
  int64_t i;
  /* setup g->ptr */
  for(i=0;i<g->nnode;i++){
    if (!g->get_size && g->start[i]) g->ptr[i] = START; // start nodes for pecolation
    else if (!g->get_size && g->stop[i]) g->ptr[i] = STOP; // stop nodes for pecolation
    else g->ptr[i] = -1; // make node isolated (own root)
  }
  /* apply losses: do in 2nd step to avoid that MEASUREZ can be overwritten */
  for(i=0;i<g->nnode;i++){
    if(pSite < (float)rand()/(float)RAND_MAX && g->static_node[i]==false){
      g->ptr[i] = LOST;
      for (uint8_t j=0;j<g->len_nb[i];j++){
        int64_t idx_nb = g->nn[i*g->num_nb_max + j];
        if(g->ptr[idx_nb] != LOST) g->ptr[idx_nb] = MEASUREZ; // mark all non-lost neighbors with Z
      }
    }
  }
  /* get largest component size or information on percolation/cluster-spanning
     in constrast to BFS, ptr could be used when further edges are added by fusion */
  int64_t s1,s2,r1,r2; // s1,s2: start node index; r1,r2: root node index
  int64_t big = 0; // size of largest component
  for (i=0;i<g->num_edges;i++){
    s1 = g->edge_node[i]; // get node index 1 of the edge
    s2 = g->nn[s1*g->num_nb_max + g->edge_idx[i]]; // get node index 2 of the edge
    if(g->ptr[s2]>MEASUREZ) r2 = findroot(g, s2); //find root and do path compression
    if(g->ptr[s1]>MEASUREZ) r1 = findroot(g, s1); //find root and do path compression
    if (g->ptr[s1]>MEASUREZ && g->ptr[s2]>MEASUREZ && r2!=r1){ // if both nodes of the edge exist and root is not the same
      if(!g->get_size && ((g->ptr[r1] == START && g->ptr[r2] == STOP) || (g->ptr[r2] == START && g->ptr[r1] == STOP))) return 1; // percolation achieved
      merge_root(g, r1, r2, &big); // update largest component size, update root (return not needed)
    }
  }
  if(g->get_size) return big; // return size of larget connected component
  else return 0; // 0 means no percolation achieved
}

/* quantum case only: apply loss to graph by picking random qubits for loss, no Newman-Ziff method (alternative to apply_loss but not really faster)
   pSite: probability that photon survives
   note: all g->fusion_node[i] must be set to 0 */
int64_t apply_loss_bfs(Graph* g, float pSite){
  int64_t i;
  /* setup g->ptr */
  for(i=0;i<g->nnode;i++){
    if (!g->get_size && g->start[i]) g->ptr[i] = START; // start nodes for pecolation
    else if (!g->get_size && g->stop[i]) g->ptr[i] = STOP; // stop nodes for pecolation
    else g->ptr[i] = -1; // make node isolated (own root)
  }
  /* apply losses: do in 2nd step to avoid that MEASUREZ can be overwritten */
  for(i=0;i<g->nnode;i++){
    if(pSite < (float)rand()/(float)RAND_MAX && g->static_node[i]==false){
      g->ptr[i] = LOST;
      for (uint8_t j=0;j<g->len_nb[i];j++){
        int64_t idx_nb = g->nn[i*g->num_nb_max + j];
        if(g->ptr[idx_nb] != LOST) g->ptr[idx_nb] = MEASUREZ; // mark all non-lost neighbors with Z
      }
    }
  }
  /* get largest component size or information on percolation/cluster-spanning */
  memset(g->fusion_node, 0, g->nnode * sizeof(bool));  // must be zero so that crawl_graph_bfs can be used for this application
  return bfs_full(g);
}

/* apply loss to existing graph state, using modified Newman-Ziff algorithm
   g: graph representing a graph state
   pBond: probability of bond (only meaningful if bond is missing due to imperfect fusion etc. and loss happens afterwards)
   considerStaticNode: if true, nodes specified as static by g->static_node cannot be lost
   numStaticNode: pointer to number of static nodes
   idxLambda: if g->get_size==false, this points to the index where the percolation transition in the return (largest_component) happens */
int64_t* apply_loss_nz(Graph* g, float pBond, bool considerStaticNode, int64_t* numStaticNode, int64_t* idxLambda){
  int64_t i,s1,s2; // s1,s2: start node index
  bool stop_loop = false;
  int64_t big = 0; // size of largest component
  int64_t* order = permutation_static_first(g->nnode, considerStaticNode, g->static_node, numStaticNode);  // random ordering of node indices
  bool* bond = malloc(g->nnode * (size_t)g->num_nb_max * sizeof(bool));
  for(i=0;i<g->nnode * g->num_nb_max;i++) bond[i] = (pBond > (float)rand()/(float)RAND_MAX); // randomly remove edges
  int64_t* largest_component = malloc((g->nnode+1) * sizeof(int64_t)); // size of largest component resp. percolation (number of nodes is position of array)
  largest_component[0] = 0; // with everything LOST size is 0 and no percolation (also 0)
  for (i=0;i<g->nnode;i++) g->ptr[i] = LOST; // start without any nodes
  for (i=0;i<g->nnode;i++){
    s1 = order[i];
    g->ptr[s1] = MEASUREZ; // add the node by changing it from LOST to MEASUREZ as default
    for (uint8_t j=0;j<g->len_nb[s1]+1;j++){ // go though all neighbors that might be affected, +1 -- j==len_nb[s1] is node itself
      s2 = (j<g->len_nb[s1] ? g->nn[s1*g->num_nb_max + j] : s1); // neighbor or node itself
      if (g->ptr[s2]==MEASUREZ && (s2==s1 || bond[get_edge_idx(g, s1, s2)])){ // neighbor must not be LOST nor exist already && bond to neighbor exists or node itself (s1==s2)
        bool nbLost = false;
        for (uint8_t k=0; k<g->len_nb[s2]; k++){ // check if any neighbor of the neighbor is LOST
          if(g->ptr[g->nn[s2*g->num_nb_max + k]]==LOST) {nbLost = true; break;}
        }
        if(!nbLost){ // no neighbor LOST
          if(g->get_size) g->ptr[s2] = -1; // isolated node, own root with size 1
          else if(g->start[s2]) g->ptr[s2] = START; // label new node as start node for percolation
          else if(g->stop[s2]) g->ptr[s2] = STOP; // label new node as stop node for percolation
          else g->ptr[s2] = -1; // isolated node, own root with size 1
          stop_loop = update_ptr_when_new_node(g, s2, &big, bond);
        }
      }
    }
    if(stop_loop) break; // stop loop when percolation is achieved
    else if(g->get_size) largest_component[i+1] = big;
    else largest_component[i+1] = 0; // no percolation yet
  }
  if(!g->get_size){
    for(int64_t k=i;k<g->nnode;k++) largest_component[k+1] = 1; // set to 1 where percolation is achieved
    *idxLambda = i+1;
  }
  free(order);
  free(bond);
  return largest_component;
}

/* apply loss to fusion network, return information about failed fusions/edges
   g: connected graph that represents fusion network of star-shaped resource states
   pSite: probability that photon survives
   nFusionPhotons: number of photons in fusion, 2 is standard fusion
   pFail: failure probability of fusion (should be adapted to nFusionPhotons) */
bool* apply_loss_to_fusion_net(Graph* g, float pSite, uint8_t nFusionPhotons, float pFail){
  int64_t i;
  /* setup g->ptr */
  for(i=0;i<g->nnode;i++){
    if (!g->get_size && g->start[i]) g->ptr[i] = START; // start nodes for pecolation
    else if (!g->get_size && g->stop[i]) g->ptr[i] = STOP; // stop nodes for pecolation
    else g->ptr[i] = -1; // make node isolated (own root)
  }
  /* treat losses on fusions (do in 2nd step to avoid that MEASUREZ can be overwritten by g->ptr[i] = START/STOP) */
  float pNoEdgeLoss = 1.0; // probability of no loss on edge
  for(uint8_t k=0; k<nFusionPhotons; k++) pNoEdgeLoss *= pSite;
  bool* fusionS = malloc(g->nnode * (size_t)g->num_nb_max * sizeof(bool)); // success of fusion
  for (i=0; i<g->nnode * (int64_t)g->num_nb_max; i++) fusionS[i] = (pFail < (float)rand()/(float)RAND_MAX);
  for (i=0;i<g->num_edges;i++){ // take care of losses in fusions
    if(pNoEdgeLoss < (float)rand()/(float)RAND_MAX){
      int64_t s1 = g->edge_node[i]; // get node 1 of the edge
      int64_t s2 = g->nn[s1*g->num_nb_max + g->edge_idx[i]]; // get node 2 of the edge
      g->ptr[s1] = MEASUREZ;
      g->ptr[s2] = MEASUREZ;
    }
  }
  /* take care of losses on sites (central qubits) and effect of losses on neighbors */
  for(i=0;i<g->nnode;i++){
    if(pSite < (float)rand()/(float)RAND_MAX && g->static_node[i]==false){
      g->ptr[i] = LOST;
      for (uint8_t j=0;j<g->len_nb[i];j++){
        int64_t idx_nb = g->nn[i*g->num_nb_max + j];
        if(fusionS[get_edge_idx(g, i, idx_nb)]){ // only if fusion succeeds, neighbor can be affected by site loss
          if(g->ptr[idx_nb] != LOST) g->ptr[idx_nb] = MEASUREZ; // mark all non-lost neighbors with Z
        }
      }
    }
  }
  return fusionS;
}

/* apply loss to fusion network and then analyze percolation, no Newman-Ziff speedup (very similar to apply_loss)
   g: connected graph that represents fusion network of star-shaped resource states
   pSite: probability that photon survives
   nFusionPhotons: number of photons in fusion, 2 is standard fusion
   pFail: failure probability of fusion (should be adapted to nFusionPhotons to be physical) */
int64_t analyse_fusion_net(Graph* g, float pSite, uint8_t nFusionPhotons, float pFail){ // TBD: missing bond option (if loss is applied after fusion failures)
  /* apply loss, get largest connected component after this step */
  bool* fusionS = apply_loss_to_fusion_net(g, pSite, nFusionPhotons, pFail);

  /* get largest component size or information about percolation/cluster-spanning
     (using and updating ptr, likely slower than BFS but ptr could be used when further edges are added by later fusions) */
  int64_t s1,s2,r1,r2; // s1,s2: start node index; r1,r2: root node index
  int64_t big = 0; // size of largest component
  for (int64_t i=0;i<g->num_edges;i++){
    s1 = g->edge_node[i]; // get node index 1 of the edge
    s2 = g->nn[s1*g->num_nb_max + g->edge_idx[i]]; // get node index 2 of the edge
    if (!fusionS[get_edge_idx(g, s1, s2)]) continue; // no connection due to fusion failure
    if(g->ptr[s2]>MEASUREZ) r2 = findroot(g, s2); //find root and do path compression
    if(g->ptr[s1]>MEASUREZ) r1 = findroot(g, s1); //find root and do path compression
    if (g->ptr[s1]>MEASUREZ && g->ptr[s2]>MEASUREZ && r2!=r1){ // if both nodes of the edge exist and root is not the same
      if(!g->get_size && ((g->ptr[r1] == START && g->ptr[r2] == STOP) || (g->ptr[r2] == START && g->ptr[r1] == STOP))){
        free(fusionS);
        return 1; // percolation achieved
      }
      merge_root(g, r1, r2, &big); // update largest component size, update root
    }
  }
  free(fusionS);
  if(g->get_size) return big; // return size of larget connected component
  else return 0; // 0 means no percolation (cluster spanning) achieved
}

/* apply loss to fusion network made of star-shaped resource states (using a modified Newman-Ziff algorithm), max 255 neighbors per node allowed
   g: connected graph that represents fusion network
   pFusion: probability that fusion secceeds
   considerStaticNode: if true, nodes specified as static by g->static_node cannot be lost
   numQbts: count all qubits in fusion network
   numStaticNode: count all static qubits that cannot be lost
   idxLambda: if g->get_size==false, this points to the index where the percolation transition in the return (largest_component) happens */
int64_t* fusion_and_loss_nz(Graph* g, float pFusion, bool considerStaticNode, int64_t* numQbts, int64_t* numStaticNode, int64_t* idxLambda){
  int64_t i, s1, s2;
  int64_t cntQbts = 0; // number of all qubits (also what is used for fusions)
  for(i=0;i<g->nnode;i++) cntQbts += (g->len_nb[i] + 1); // +1 for central node itself
  bool* static_node;
  if(considerStaticNode){
    static_node = malloc(cntQbts * sizeof(bool)); // g->static_node cannot be used as cntQbts > g->nnode
    memset(static_node, 0, cntQbts * sizeof(bool));
  }
  int64_t* centerIdx = malloc(cntQbts * sizeof(int64_t)); // index of central qubit
  uint8_t* nbPos = malloc(cntQbts); // position where the neighbor node is in g->nn
  int64_t* largest_component = malloc((cntQbts+1) * sizeof(int64_t));
  largest_component[0] = 0; // with all lost size is 0 and there is no percolation (also 0)
  int64_t cnt = 0;
  for(i=0;i<g->nnode;i++){
    for(uint8_t k=0; k<g->len_nb[i]; k++){
      centerIdx[cnt] = i;
      nbPos[cnt] = k;
      cnt++;
    }
    centerIdx[cnt] = i;
    nbPos[cnt] = 255; // 255 refers to central node itself
    if(considerStaticNode){static_node[cnt] = g->static_node[i];}
    cnt++;
  }
  int64_t* order;
  if(considerStaticNode){
    order = permutation_static_first(cntQbts, considerStaticNode, static_node, numStaticNode); // random ordering of indices (static nodes come first)
    free(static_node);
  } else{
    order = permutation(cntQbts); // random ordering of photon indices
    *numStaticNode = 0;
  }

  for (i=0;i<g->nnode;i++) g->ptr[i] = LOST; // start without any central nodes (add also static nodes one-by-one)
  bool* fusionQ = malloc(g->nnode * (size_t)g->num_nb_max * sizeof(bool)); // is this fusion qubits lost?
  for(i=0;i<g->nnode * g->num_nb_max;i++) fusionQ[i] = 0; // start with all fusion qubits lost (0)?
  bool* fusionS = malloc(g->nnode * (size_t)g->num_nb_max * sizeof(bool)); // success of fusion
  for(i=0;i<g->nnode * g->num_nb_max;i++) fusionS[i] = (pFusion > (float)rand()/(float)RAND_MAX); // randomly set fusion failures

  // add the qubits one-by-one
  int64_t* listNodes = malloc(255 * sizeof(int64_t)); // list all neighbors etc. that might have changed
  uint8_t numNodes; // number of nodes that might have changed
  bool stop_loop = false; // stop when percolation achieved
  int64_t big = 0; // size of largest component
  for (i=0;i<cntQbts;i++){
    // apply the change: make the corresponding qubit not lost, define list of all qubits that might have changed
    s1 = centerIdx[order[i]];
    uint8_t s1nbPos = nbPos[order[i]];
    if(s1nbPos == 255){ // center node
      g->ptr[s1] = MEASUREZ; // center qubit not LOST anymore
      numNodes = g->len_nb[s1] + 1; // all neighbors and node itself (+1) might change
      listNodes[0] = s1;
      for(uint8_t k=1; k<numNodes; k++) listNodes[k] = g->nn[s1 * g->num_nb_max + k-1];
    } else { // photon on an edge
      fusionQ[s1 * (size_t)g->num_nb_max + s1nbPos] = 1; // fusion qubit not lost anymore
      listNodes[0] = s1;
      listNodes[1] = g->nn[s1 * g->num_nb_max + s1nbPos];
      numNodes = 2; // only the two nodes that share the edge on which the lost photon is might be affected
    }

    // go through all qubits that might have been affected and update ptr
    for (uint8_t j=0;j<numNodes;j++){ // go though all nodes that might be affected
      s2 = listNodes[j];
      if (g->ptr[s2]==MEASUREZ){ // if neighbor is not LOST nor exist already, then something might change
        bool nbLost = false;
        for (uint8_t k=0;k<g->len_nb[s2];k++){ // check all edges (with their fusion qubits) one by one
          int64_t nb = g->nn[s2 * (size_t)g->num_nb_max + k];
          if(fusionQ[s2 * (size_t)g->num_nb_max + k]==0 || (g->ptr[nb]==LOST && fusionS[get_edge_idx(g, s2, nb)])){ // check fusion qubit lost or successful connection to LOST node? (if 2nd fusion qubit lost -> does not matter if this is true)
            nbLost = true; break;
          } else { // 2nd fusion photon on edge lost? (the one that belongs to the neighboring resource state)
            uint8_t s;
            for(s=0; g->nn[nb * (size_t)g->num_nb_max + s] != s2; s++) {} // find fusion partner, TBD (maybe): faster with modified data structure
            if (fusionQ[nb * (size_t)g->num_nb_max + s]==0) {
              nbLost = true; break;
            }
          }
        }
        if(!nbLost){ // no neighbor of fusion qubit LOST
          if(g->get_size) g->ptr[s2] = -1; // isolated node, own root with size 1
          else if(g->start[s2]) g->ptr[s2] = START; // label new node as start node for percolation
          else if(g->stop[s2]) g->ptr[s2] = STOP; // label new node as stop node for percolation
          else g->ptr[s2] = -1; // isolated node, own root with size 1
          stop_loop = update_ptr_when_new_node(g, s2, &big, fusionS);
        }
      }
    }

    if(stop_loop) break; // stop loop when percolation is achieved
    else if(g->get_size) largest_component[i+1] = big;
    else largest_component[i+1] = 0; // no percolation yet
  }
  if(!g->get_size){
    for(int64_t k=i;k<cntQbts;k++) largest_component[k+1] = 1; // set to 1 where percolation is achieved
    *idxLambda = i+1;
  }

  free(listNodes);
  free(order);
  free(centerIdx);
  free(nbPos);
  free(fusionQ);
  free(fusionS);
  *numQbts = cntQbts;
  return largest_component;
}

/* give a random number specifying after how many attempts a fusion succeeds (if there is not loss)
   the frequency distribution follows and exponential decay with (1-pFusion)^n, 255 means that none of the fusions succeeds
   pFusion: fusion success probability, n_fusion_max: maximum number of fusions */
static uint8_t get_fusion_success_level(float pFusion, uint8_t n_fusion_max){
  for(uint8_t n=0; n<n_fusion_max; n++){
    if(pFusion > (float)rand()/(float)RAND_MAX) return n+1;
  }
  return 255; // 255 means that no fusion succeeds
}

/* apply loss to fusion network made of star-shaped resource states (using a modified Newman-Ziff algorithm),
   allow locally adaptive repetition of fusion (makes only sense when the central node is a quantum emitter), max 255 neighbors allowed, max 254 fusion repetitions allowed
   g: connected graph that represents fusion network
   pFusion: probability that a single fusion secceeds
   n_fusion_max: maximum number of fusion attempts
   numQbts: count all qubits in fusion network
   numStaticNode: count all static qubits that cannot be lost
   idxLambda: if g->get_size==false, this points to the index where the percolation transition in the return (largest_component) happens */
int64_t* fusion_repeat_and_loss_nz(Graph* g, float pFusion, uint8_t n_fusion_max, int64_t* numQbts, int64_t* numStaticNode, int64_t* idxLambda){
  int64_t i, s1, s2;
  int64_t cntQbts = g->nnode; // number of all qubits (including repetition of fusions), add fusion qubits later
  uint8_t* arr_n_fusion = malloc(g->nnode * (int64_t)g->num_nb_max);  // store how often the fusion involving a particular fusion qubit is repeated
  memset(arr_n_fusion, 0, g->nnode * (int64_t)g->num_nb_max);
  uint8_t* n_no_loss = malloc(g->nnode * (int64_t)g->num_nb_max); // number of repeated fusion qubits that are not lost
  memset(n_no_loss, 0, g->nnode * (int64_t)g->num_nb_max); // assume all fusion qubits lost initially

  /* set the number of fusion attempts for the case of no photon loss */
  for(i=0;i<g->nnode;i++){
    for(uint8_t k=0; k<g->len_nb[i]; k++){
      int64_t idx = i * (size_t)g->num_nb_max + k;
      if(arr_n_fusion[idx] == 0){ // if number of fusion attempts is not set for this edge yet
        arr_n_fusion[idx] = get_fusion_success_level(pFusion, n_fusion_max);
        uint8_t s;
        int64_t nb = g->nn[idx]; // the qubit connected to the fusion partner
        for(s=0; g->nn[nb * (size_t)g->num_nb_max + s] != i; s++) {} // find fusion partner, TBD (maybe): faster with modified data structure
        arr_n_fusion[nb * (size_t)g->num_nb_max + s] = arr_n_fusion[idx]; // set the number of repetitions to same value for the fusion partner
        if(arr_n_fusion[idx] < 255){ // 255 means that all n_fusion_max fusions don't succeed
          cntQbts += (int64_t)arr_n_fusion[idx] * 2; // count both qubits on edge times number of fusion repetitions
        } else {
          cntQbts += (int64_t)n_fusion_max * 2; // qubits from all failed fusion attempts
        }
      }
    }
  }

  /* data structure for looping through all lossy qubits in random order */
  bool* static_node = malloc(cntQbts * sizeof(bool)); // g->static_node cannot be used as cntQbts > g->nnode
  memset(static_node, 0, cntQbts * sizeof(bool));
  int64_t* centerIdx = malloc(cntQbts * sizeof(int64_t)); // index of central qubit
  uint8_t* nbPos = malloc(cntQbts); // position where the neighbor node is in g->nn
  int64_t* largest_component = malloc((cntQbts+1) * sizeof(int64_t));
  largest_component[0] = 0; // with all lost size is 0 and there is no percolation (also 0)
  int64_t cnt = 0;
  for(i=0;i<g->nnode;i++){
    for(uint8_t k=0; k<g->len_nb[i]; k++){ // loop all fusion qubits of the central node
      for(uint8_t rep=0; rep < arr_n_fusion[i * g->num_nb_max + k] && rep < n_fusion_max; rep++){ // loop all repetitions of the fusion qubit
        centerIdx[cnt] = i;
        nbPos[cnt] = k;
        cnt++;
      }
    }
    centerIdx[cnt] = i;
    nbPos[cnt] = 255; // 255 refers to central node itself
    static_node[cnt] = g->static_node[i];
    cnt++;
  }
  int64_t* order = permutation_static_first(cntQbts, true, static_node, numStaticNode); // random ordering of indices (static nodes come first)
  free(static_node);

  for (i=0;i<g->nnode;i++) g->ptr[i] = LOST; // start without any central nodes (add also static nodes one-by-one)

  /* add the qubits one-by-one */
  int64_t* listNodes = malloc(255 * sizeof(int64_t)); // list all neighbors etc. that might have changed
  uint8_t numNodes; // number of nodes that might have changed
  bool stop_loop = false; // stop when percolation achieved
  int64_t big = 0; // size of largest component
  for (i=0;i<cntQbts;i++){
    /* apply the change: make the corresponding qubit not lost, define list of all qubits that might have changed */
    s1 = centerIdx[order[i]];
    uint8_t s1nbPos = nbPos[order[i]];
    if(s1nbPos == 255){ // center node
      g->ptr[s1] = MEASUREZ; // center qubit not LOST anymore
      numNodes = g->len_nb[s1] + 1; // all neighbors and node itself (+1) might change
      listNodes[0] = s1;
      for(uint8_t k=1; k<numNodes; k++) listNodes[k] = g->nn[s1 * g->num_nb_max + k-1];
    } else { // photon on an edge
      int64_t idx = s1 * (size_t)g->num_nb_max + s1nbPos;
      n_no_loss[idx] += 1; // one repeated fusion qubit less is lost
      if(n_no_loss[idx] < arr_n_fusion[idx] && n_no_loss[idx] < n_fusion_max){
        largest_component[i+1] = (g->get_size ? big : 0); // size/percolation cannot change
        continue; // nothing can change if there is still one repeated fusion qubit lost
      } else {
        uint8_t s;
        int64_t nb = g->nn[idx]; // the qubit connected to the fusion partner
        for(s=0; g->nn[nb * (size_t)g->num_nb_max + s] != s1; s++) {} // find fusion partner, TBD (maybe): faster with modified data structure
        int64_t idx_nb = nb * (size_t)g->num_nb_max + s;
        if(n_no_loss[idx_nb] < arr_n_fusion[idx_nb] && n_no_loss[idx_nb] < n_fusion_max){
          largest_component[i+1] = (g->get_size ? big : 0); // size/percolation cannot change
          continue; // nothing can change if there is still one repeated fusion qubit lost
        } else {
          listNodes[0] = s1;
          listNodes[1] = nb;
          numNodes = 2; // only the two nodes that share the edge on which the lost photon is might be affected
        }
      }
    }

    /* go through all qubits that might have been affected and update ptr */
    for (uint8_t j=0;j<numNodes;j++){ // go though all nodes that might be affected
      s2 = listNodes[j];
      if (g->ptr[s2]==MEASUREZ){ // neighbor must not be LOST nor exist already
        bool nbLost = false;
        for (uint8_t k=0;k<g->len_nb[s2];k++){ // check all edges (with their fusion qubits) one by one
          int64_t idx = s2 * (size_t)g->num_nb_max + k;
          int64_t nb = g->nn[idx];
          if((n_no_loss[idx] < arr_n_fusion[idx] && n_no_loss[idx] < n_fusion_max) || (g->ptr[nb]==LOST && arr_n_fusion[idx] < 255)){ // check fusion qubit lost or successful connection to LOST node? (<255 means successful fusion)
            nbLost = true; break;
          } else { // 2nd fusion photon on edge lost? (the one that belongs to the neighboring resource state)
            uint8_t s;
            for(s=0; g->nn[nb * (size_t)g->num_nb_max + s] != s2; s++) {} // find fusion partner, TBD (maybe): faster with modified data structure
            if (n_no_loss[nb * (size_t)g->num_nb_max + s] < arr_n_fusion[nb * (size_t)g->num_nb_max + s] && n_no_loss[nb * (size_t)g->num_nb_max + s] < n_fusion_max) {
              nbLost = true; break;
            }
          }
        }
        if(!nbLost){ // no neighbor or fusion qubit LOST
          if(g->get_size) g->ptr[s2] = -1; // isolated node, own root with size 1
          else if(g->start[s2]) g->ptr[s2] = START; // label new node as start node for percolation
          else if(g->stop[s2]) g->ptr[s2] = STOP; // label new node as stop node for percolation
          else g->ptr[s2] = -1; // isolated node, own root with size 1
          stop_loop = update_ptr_when_new_node_b(g, s2, &big, arr_n_fusion);
        }
      }
    }

    if(stop_loop) break; // stop loop when percolation is achieved
    else if(g->get_size) largest_component[i+1] = big;
    else largest_component[i+1] = 0; // no percolation yet
  }
  if(!g->get_size){
    for(int64_t k=i;k<cntQbts;k++) largest_component[k+1] = 1; // set to 1 where percolation is achieved
    *idxLambda = i+1;
   }

  free(listNodes);
  free(arr_n_fusion);
  free(n_no_loss);
  free(order);
  free(centerIdx);
  free(nbPos);
  *numQbts = cntQbts;
  return largest_component;
}

/* get expectation value of largest connected component (see Eq. 10) */
float* get_expectation_value(float* sweep_prob, int sweep_len, int64_t* largest_component, int64_t N, int64_t nStatic, int64_t size_max, bool approx, bool mute){
  // sweep_prob: probabilities that a node/bond/photon/fusion exists, all sweep_prob[i] must be in range (0-1)
  // sweep_len: length of array with probabilities
  // largest_component: largest component size for different numbers of edges/nodes etc.
  // N: length of largest_component
  // nStatic: number of static nodes/bonds that are always present (offset)
  // size_max: maximum value of largest_component
  // approx: cut loop when bnnp becomes too small to contribute
  // mute: if true it does not print anything
  if(nStatic>0){
    N = N - nStatic;
    largest_component += nStatic; // add offeset to pointer
  }
  int64_t mean, n, n_min, n_max;
  float sum;
  float* bnnp;
  bnnp = malloc((N+1) * sizeof(float));
  float* expectation_value;
  expectation_value = malloc(sweep_len * sizeof(float));
  memset(expectation_value, 0, sweep_len * sizeof(float));
  for(int i=0;i<sweep_len;i++){
    mean = (int64_t)(N * sweep_prob[i]);
    bnnp[mean] = 1;
    sum = 0;
    if(approx){
      for(n=mean+1;n<N+1;n++){
        bnnp[n] = bnnp[n-1] * (N-n+1)/n*sweep_prob[i]/(1-sweep_prob[i]);
        if(bnnp[n] * N < bnnp[mean] / 1000) break; // heuristic cutoff
      }
      n_max = n - 1;
      for(n=mean-1;n>=0;n--){
        bnnp[n] = bnnp[n+1] * (n+1)/(N-n)*(1-sweep_prob[i])/sweep_prob[i];
        if(bnnp[n] * N < bnnp[mean] / 1000) break; // heuristic cutoff
      }
      n_min = n + 1;
      for(n=n_min;n<=n_max;n++) sum += bnnp[n];
      for(n=n_min;n<=n_max;n++) expectation_value[i] += bnnp[n]/sum * largest_component[n]/size_max;
    } else{
      for(n=mean+1;n<N+1;n++) bnnp[n] = bnnp[n-1] * (N-n+1)/n*sweep_prob[i]/(1-sweep_prob[i]);
      for(n=mean-1;n>=0;n--) bnnp[n] = bnnp[n+1] * (n+1)/(N-n)*(1-sweep_prob[i])/sweep_prob[i];
      for(n=0;n<N+1;n++) sum += bnnp[n];
      for(n=0;n<N+1;n++) expectation_value[i] += bnnp[n]/sum * largest_component[n]/size_max;
    }
    if(mute==false) printf("%f %f\n", sweep_prob[i], expectation_value[i]);
  }
  free(bnnp);
  return expectation_value;
}

