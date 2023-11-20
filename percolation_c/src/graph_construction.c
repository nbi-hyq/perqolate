#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "../inc/global.h"
#include "../inc/graph_type.h"
#include "../inc/percolation_main.h"

/* append an edge to array (for keeping track of edges when looping) */
static void append_edge(Graph* g, int64_t i, uint8_t idx){
  if (i<g->nn[i*g->num_nb_max + idx]){ // avoid double counting
    g->edge_node[g->num_edges] = i;
    g->edge_idx[g->num_edges] = idx;
    g->num_edges++;
  }
}

/* create array representation with all the edges (call when g->nn is created) */
static void make_edge_array(Graph* g){
  for (int64_t i=0; i<g->nnode; i++){
    for (uint8_t idx=0;idx<g->len_nb[i];idx++) append_edge(g, i, idx);
  }
}

/* helper function to convert index to arbitrary digital system representation */
static void index_to_digital_vec(int64_t idx, int64_t* vec, uint8_t n_dim, int64_t lsize){
  // idx: decimal index, vec: stores the digital system representation, n_dim: number of dimensions, lsize: size in each dimension
  int64_t factor = 1;
  memset(vec, 0, (size_t)n_dim*sizeof(int64_t));
  for(int pos=n_dim-1;pos>=0;pos--){
    vec[pos] = (idx % (factor * lsize) / factor);
    idx -= vec[pos] * factor;
    factor *= lsize;
  }
}

/* helper function to convert index to arbitrary digital system representation (one digital position is skipped) */
static void index_to_digital_vec_skip(int64_t idx, int64_t* vec, uint8_t n_dim, uint8_t pos_fix, int64_t val_fix, int64_t lsize){
  // idx: decimal index, vec: stores the digital system representation, n_dim: number of dimensions, lsize: size in each dimension, pos_fix: skipped position, pos_fix: value at skipped position
  int64_t factor = 1;
  memset(vec, 0, (size_t)n_dim*sizeof(int64_t));
  for(int pos=n_dim-1;pos>=0;pos--){
    if(pos != pos_fix){
      vec[pos] = (idx % (factor * lsize) / factor);
      idx -= vec[pos] * factor;
      factor *= lsize;
    } else {
      vec[pos] = val_fix;
    }
  }
}

/* helper function to convert digital system representation back to decimal int */
static int64_t digital_vec_to_index(int64_t* vec, uint8_t n_dim, int64_t lsize){
  int64_t idx = 0;
  int64_t factor = 1;
  for(int pos=n_dim-1;pos>=0;pos--){
    idx += vec[pos] * factor;
    factor *= lsize;
  }
  return idx;
}

/* create 2d square lattice */
Graph get_2d_square(int64_t lsize, bool periodic, bool get_size, bool edge_list){
  int64_t i;
  uint8_t num_nb_max = 4;
  Graph g = new_graph(lsize*lsize, num_nb_max, get_size, edge_list, false);
  if (periodic) { // construct graph with periodic boundaries
    for (i=0; i<g.nnode; i++){
      g.nn[i*g.num_nb_max + 0] = (i+1)%g.nnode;
      g.nn[i*g.num_nb_max + 1] = (i+g.nnode-1)%g.nnode;
      g.nn[i*g.num_nb_max + 2] = (i+lsize)%g.nnode;
      g.nn[i*g.num_nb_max + 3] = (i+g.nnode-lsize)%g.nnode;
      if (i%lsize == 0) g.nn[i*g.num_nb_max + 1] = i+lsize-1;
      if ((i+1)%lsize == 0) g.nn[i*g.num_nb_max + 0] = i-lsize+1;
      g.len_nb[i] = num_nb_max;
    }
  } else { // construct graph without periodic boundaries
    int64_t x, y;
    uint8_t idx;
    for (x=0; x<lsize; x++){
      for (y=0; y<lsize; y++){
      idx = 0;
      i = x * lsize + y;
      if(y<lsize-1) {g.nn[i*num_nb_max + idx] = i + 1; idx++;}
      if(y>0) {g.nn[i*num_nb_max + idx] = i - 1; idx++;}
      if(x<lsize-1) {g.nn[i*num_nb_max + idx] = i+lsize; idx++;}
      if(x>0) {g.nn[i*num_nb_max + idx] = i-lsize; idx++;}
      g.len_nb[i] = idx;
      }
    }
  }
  if(edge_list) make_edge_array(&g);
  memset(g.static_node, 0, sizeof(bool)*g.nnode); // make all nodes non-static
  if(!get_size){
    memset(g.start,0,g.nnode*sizeof(bool));
    memset(g.stop,0,g.nnode*sizeof(bool));
    for (i=0; i<lsize; i++) g.start[i] = 1;
    for (i=(lsize-1)*lsize; i<g.nnode; i++) g.stop[i] = 1;
  }
  return g;
}

/* create 2d triangular lattice (no periodic boudaries) */
Graph get_2d_triangular(int64_t lsize, bool periodic, bool get_size, bool edge_list){
  int64_t x, y, i;
  uint8_t idx;
  uint8_t num_nb_max = 6;
  Graph g = new_graph(lsize * lsize, num_nb_max, get_size, edge_list, false);
  for (x=0; x<lsize; x++){
    for (y=0; y<lsize; y++){
      idx = 0;
      i = x * lsize + y;
      if(y<lsize-1){ g.nn[i*num_nb_max + idx] = i + 1; idx++;}
      else if(periodic){ g.nn[i*num_nb_max + idx] = i - (lsize-1); idx++;}
      if(y>0){ g.nn[i*num_nb_max + idx] = i - 1; idx++;}
      else if(periodic){ g.nn[i*num_nb_max + idx] = i + (lsize-1); idx++;}
      if(x<lsize-1){ g.nn[i*num_nb_max + idx] = i+lsize; idx++;}
      else if(periodic){ g.nn[i*num_nb_max + idx] = y; idx++;}
      if(x>0){ g.nn[i*num_nb_max + idx] = i-lsize; idx++;}
      else if(periodic){ g.nn[i*num_nb_max + idx] = (lsize-1) * lsize + y; idx++;}
      if(x<lsize-1 && y<lsize-1){ g.nn[i*num_nb_max + idx] = i+lsize+1; idx++;}
      else if(periodic){ g.nn[i*num_nb_max + idx] = ((x+1)%lsize) * lsize + ((y+1)%lsize); idx++;}
      if(x>0 && y>0){ g.nn[i*num_nb_max + idx] = i-lsize-1; idx++;}
      else if(periodic){ g.nn[i*num_nb_max + idx] = ((x-1)%lsize) * lsize + ((y-1)%lsize); idx++;}
      g.len_nb[i] = idx;
    }
  }
  for (x=0; x<lsize; x++){
    for (y=0; y<lsize; y++){
      i = x * lsize + y;
      if(edge_list){
        for(idx=0;idx<g.len_nb[i];idx++) append_edge(&g, i, idx);
      }
      g.static_node[i] = false; // make all nodes non-static
    }
  }
  if(!get_size){
    memset(g.start,0,g.nnode*sizeof(bool));
    memset(g.stop,0,g.nnode*sizeof(bool));
    for (i=0; i<lsize; i++) g.start[i] = 1;
    for (i=(lsize-1)*lsize; i<g.nnode; i++) g.stop[i] = 1;
  }
  return g;
}

/* recursively add node and connection to child nodes */
static void add_node(Graph* g, uint8_t* num_child, int64_t* node_idx_next, uint8_t depth, uint8_t num_nb_max){
  uint8_t idx_offset = (depth ? 1 : 0); // all nodes except root already have one connection from above
  int64_t node_idx_here = *node_idx_next; // save current node
  if(g->get_size == false && num_child[depth] == 0){
    g->stop[node_idx_here] = 1; // stop at leave nodes
  }
  uint8_t idx; // child number index in g->nn
  for(idx=idx_offset; idx<num_child[depth]+idx_offset; idx++){
    *node_idx_next += 1; // child is next node
    g->nn[node_idx_here * num_nb_max + idx] = *node_idx_next;
    g->nn[(*node_idx_next) * num_nb_max] = node_idx_here;
    add_node(g, num_child, node_idx_next, depth+1, num_nb_max);
  }
  g->len_nb[node_idx_here] = idx;
}

/* create tree lattice (can be for instance the Bethe lattice) */
Graph get_tree_lattice(uint8_t depth_max, uint8_t* num_child, bool static_center, bool get_size, bool edge_list){
  int64_t nnode = 0;
  int64_t nnode_at_depth = 1;
  uint8_t num_nb_max = 0;
  for(uint8_t i=0; i<=depth_max; i++){
    nnode += nnode_at_depth;
    nnode_at_depth *= num_child[i];
    if(num_nb_max < num_child[i] + (i>0)){
      num_nb_max = num_child[i] + (i>0); // +1 for link to node above (except root node at depth==0)
    }
  }
  Graph g = new_graph(nnode, num_nb_max, get_size, edge_list, false);
  if(!get_size){
    memset(g.start,0,g.nnode*sizeof(bool));
    memset(g.stop,0,g.nnode*sizeof(bool));
    g.start[0] = 1; // start at root node
  }
  int64_t node_idx_next = 0;
  add_node(&g, num_child, &node_idx_next, 0, num_nb_max); // create lattice recursively (depth first)
  for (int64_t i=0; i<nnode; i++){
    if(edge_list){
      for(uint8_t idx=0;idx<g.len_nb[i];idx++) append_edge(&g, i, idx);
    }
    g.static_node[i] = (static_center ? true : false);
  }
  return g;
}

/* create an n-dimensional simple cubic lattice (option: periodic boundary conditions) */
Graph get_nd_simple_cubic(int64_t lsize, uint8_t n_dim, bool static_center, bool periodic, bool get_size, bool edge_list){
  int64_t numNodes = 1;
  uint8_t num_nb_max = 2*n_dim; // maximum number of neighbors
  for(uint8_t j=0;j<n_dim;j++) numNodes *= lsize; // number of nodes
  Graph g = new_graph(numNodes, num_nb_max, get_size, edge_list, false);
  int64_t i_vec[n_dim]; // digital system representation
  uint8_t digit_pos, nbIdx;
  int64_t iNb; // index of neighbor
  for(int64_t i_node=0;i_node<numNodes;i_node++){
    index_to_digital_vec(i_node, i_vec, n_dim, lsize); // get digital system representation (i_vec) of i_node
    g.static_node[i_node] = (static_center ? true : false);
    nbIdx = 0;
    for(digit_pos=0;digit_pos<n_dim;digit_pos++){ // make connections to all the neighbors
      if(i_vec[digit_pos] < lsize - 1){
        i_vec[digit_pos] += 1;
        iNb = digital_vec_to_index(i_vec, n_dim, lsize);
        i_vec[digit_pos] -= 1;
        g.nn[i_node * num_nb_max + (nbIdx++)] = iNb; // connect to neighbor
        if(edge_list) append_edge(&g, i_node, nbIdx-1);
      } else if (periodic) {
        int64_t save = i_vec[digit_pos];
        i_vec[digit_pos] = 0;
        iNb = digital_vec_to_index(i_vec, n_dim, lsize);
        i_vec[digit_pos] = save;
        g.nn[i_node * num_nb_max + (nbIdx++)] = iNb; // connect to neighbor
        if(edge_list) append_edge(&g, i_node, nbIdx-1);
      }
      if(i_vec[digit_pos] > 0){
        i_vec[digit_pos] -= 1;
        iNb = digital_vec_to_index(i_vec, n_dim, lsize);
        i_vec[digit_pos] += 1;
        g.nn[i_node * num_nb_max + (nbIdx++)] = iNb; // connect to neighbor
        if(edge_list) append_edge(&g, i_node, nbIdx-1);
      } else if (periodic) {
        int64_t save = i_vec[digit_pos];
        i_vec[digit_pos] = lsize-1;
        iNb = digital_vec_to_index(i_vec, n_dim, lsize);
        i_vec[digit_pos] = save;
        g.nn[i_node * num_nb_max + (nbIdx++)] = iNb; // connect to neighbor
        if(edge_list) append_edge(&g, i_node, nbIdx-1);
      }
    }
    g.len_nb[i_node] = nbIdx; // number of neighbors
  }
  if(!get_size){
    memset(g.start,0,g.nnode*sizeof(bool));
    memset(g.stop,0,g.nnode*sizeof(bool));
    for (int64_t i_node=0; i_node<numNodes/lsize; i_node++) g.start[i_node] = 1;
    for (int64_t i_node=numNodes/lsize*(lsize-1); i_node<numNodes; i_node++) g.stop[i_node] = 1;
  }
  return g;
}

/* create an n-dimensional simple cubic lattice with block structure (between the blocks there is one layer of qubits that is only connected along the normal of the interface)
   between neighboring blocks there is a layer of qubits that is only connected in the direction perpendicular to the common face
   lsize: number of blocks in one dimension
   lblock: length of one block (without 1 layer between the blocks) */
Graph get_nd_simple_cubic_block(int64_t lsize, int lblock, uint8_t n_dim, bool static_center, bool periodic, bool get_size, bool edge_list){
  int64_t numNodes = 1;
  int64_t lTotal = (lblock + 1)*lsize; // totat length of lattice in one dimension
  uint8_t num_nb_max = 2*n_dim; // maximum number of neighbors
  for(uint8_t j=0; j<n_dim; j++) numNodes *= lTotal; // number of nodes
  Graph g = new_graph(numNodes, num_nb_max, get_size, edge_list, false);
  int64_t i_vec[n_dim]; // digital system representation
  for(int64_t i_node=0;i_node<numNodes;i_node++){
    index_to_digital_vec(i_node, i_vec, n_dim, lTotal); // get digital system representation (i_vec) of i_node
    uint8_t nbIdx = 0;
    uint8_t d_start = 0;
    uint8_t d_stop = n_dim-1;
    bool face = false;
    for(uint8_t d=0; d<n_dim; d++){
      if((i_vec[d]+1) % (lblock+1) == 0){
        if(face == false){
          d_start = d; d_stop = d; // make connection only normal to face
          face = true;
        } else {
          d_start = 1; d_stop = 0; // make no connections where two or more faces cross
          break;
        }
      }
    }
    for(uint8_t d = d_start; d <= d_stop; d++){ // make connections to all the neighbors
      if(i_vec[d] < lTotal - 1){
        i_vec[d] += 1;
        g.nn[i_node * num_nb_max + (nbIdx++)] = digital_vec_to_index(i_vec, n_dim, lTotal); // connect to neighbor
        i_vec[d] -= 1;
      } else if (periodic) {
        int64_t save = i_vec[d];
        i_vec[d] = 0;
        g.nn[i_node * num_nb_max + (nbIdx++)] = digital_vec_to_index(i_vec, n_dim, lTotal); // connect to neighbor
        i_vec[d] = save;
      }
      if(i_vec[d] > 0){
        i_vec[d] -= 1;
        g.nn[i_node * num_nb_max + (nbIdx++)] = digital_vec_to_index(i_vec, n_dim, lTotal); // connect to neighbor
        i_vec[d] += 1;
      } else if (periodic) {
        int64_t save = i_vec[d];
        i_vec[d] = lTotal - 1;
        g.nn[i_node * num_nb_max + (nbIdx++)] = digital_vec_to_index(i_vec, n_dim, lTotal); // connect to neighbor
        i_vec[d] = save;
      }
    }
    g.len_nb[i_node] = nbIdx; // number of neighbors
  }
  if(edge_list) make_edge_array(&g);
  if(static_center){
    memset(g.static_node, 1, sizeof(bool)*g.nnode);
  } else{
    memset(g.static_node, 0, sizeof(bool)*g.nnode);
  }
  if(!get_size){
    memset(g.start,0,g.nnode*sizeof(bool));
    memset(g.stop,0,g.nnode*sizeof(bool));
    for (int64_t i_node=0; i_node<numNodes/lTotal; i_node++) g.start[i_node] = 1;
    for (int64_t i_node=numNodes/lTotal*(lTotal-1); i_node<numNodes; i_node++) g.stop[i_node] = 1;
  }
  return g;
}

/* helper function to convert vector specifying the block and vector specifying node position locally in block to global node index */
static int64_t get_global_idx(int64_t* i_vec_rn, int64_t* i_vec_blk, uint8_t n_dim, int64_t lblock, int64_t lTotal){
  int64_t i_vec_global[n_dim]; // global node position
  for(uint8_t i=0; i < n_dim; i++) i_vec_global[i] = i_vec_blk[i] + (lblock+1)*i_vec_rn[i];
  int64_t iNode = digital_vec_to_index(i_vec_global, n_dim, lTotal);
  return iNode;
}

/* renormalize fusion network created by "get_nd_simple_cubic_block" into a faulty simple cubic lattice (experimental)
   g: graph that is renormalized (g->get_size must be true)
   pSite: probability that photon survives
   nFusionPhotons: number of photons in fusion, 2 is standard fusion
   pFail: failure probability of fusion (should be adapted to nFusionPhotons)
   lblock: size of block for renormalization (block is larger by +1 for interface connecting to other block)
   lsize: size in blocks
   periodic: periodic boundary conditions
   get_size: largest component size (1) or cluster spanning (0) in renormalized lattice
   bondFraction: fraction of existing bonds in the renormalized lattice (can be used to estimate if renormalized lattice percolates)
   bestNode: node in block that becomes part of renormalized lattice (indexed by index of the block) */
Graph renormalize_block_fusion_net(Graph* g, float pSite, uint8_t nFusionPhotons, float pFail, int lblock, int64_t lsize, uint8_t n_dim, bool periodic, bool get_size, double* bondFraction, int64_t** bestNode){
  bool* fusionS = apply_loss_to_fusion_net(g, pSite, nFusionPhotons, pFail);
  int64_t nnodeBlk = 1; // number of nodes in one block
  for(uint8_t i=0; i<n_dim; i++) nnodeBlk *= (lblock+1);
  int64_t lTotal = lsize*(lblock+1); // length of the full lattice

  /* create graph for renormalized lattice */
  int64_t nnode_rn = 1; // number of nodes in renormalized lattice
  for(uint8_t j=0; j<n_dim; j++) nnode_rn *= lsize;
  bool edge_list = true; // if one wants to loop edges later
  bool fusion_list = true; // if one wants to loop edges later
  Graph rn = new_graph(nnode_rn, g->num_nb_max, get_size, edge_list, fusion_list); // create renormalized lattice
  memset(rn.len_nb, 0, rn.nnode); // start with empty graph
  memset(rn.fusion_node, 0, sizeof(bool)*rn.nnode); // only relevant for some functions applied to renormalized graph afterwards
  if(!rn.get_size){
    memset(rn.start,0,rn.nnode*sizeof(bool));
    memset(rn.stop,0,rn.nnode*sizeof(bool));
    for (int64_t i_node=0; i_node<rn.nnode/lsize; i_node++) rn.start[i_node] = 1;
    for (int64_t i_node=rn.nnode/lsize*(lsize-1); i_node<rn.nnode; i_node++) rn.stop[i_node] = 1;
  }
  for(int64_t i=0;i < rn.nnode;i++){
    if (!rn.get_size && rn.start[i]) rn.ptr[i] = START; // start nodes for pecolation
    else if (!rn.get_size && rn.stop[i]) rn.ptr[i] = STOP; // stop nodes for pecolation
    else rn.ptr[i] = -1; // make node isolated (own root)
  }
  memset(rn.fusion_node, 0, rn.nnode * sizeof(bool)); // must be 0 when using bfs_full()

  /* renormalization strategy: find largest connected component within blocks and keep the node with largest vertex degree */
  *bestNode = malloc(sizeof(int64_t) * rn.nnode);
  int64_t* sizeMax = malloc(sizeof(int64_t) * rn.nnode); // size of the largest connected component in the corresponding block
  memset(sizeMax, 0, sizeof(int64_t) * rn.nnode);
  int64_t i_vec_rn[n_dim]; // block position
  int64_t i_vec_blk[n_dim]; // node position within block
  for(int64_t iNode_rn=0; iNode_rn < rn.nnode; iNode_rn++){ // loop over all blocks
    /* find largest connected component within block */
    index_to_digital_vec(iNode_rn, i_vec_rn, n_dim, lsize);
    for(int64_t iNode_blk=0; iNode_blk < nnodeBlk; iNode_blk++){
      index_to_digital_vec(iNode_blk, i_vec_blk, n_dim, lblock+1);
      int64_t iNode = get_global_idx(i_vec_rn, i_vec_blk, n_dim, lblock, lTotal);
      if(g->ptr[iNode] <= MEASUREZ || g->len_nb[iNode] <= 2) continue; // g->len_nb[idx_nb] <= 2: qubits at the block interfaces are not counted
      int64_t r1 = findroot(g, iNode); // might be connected to neighbor from previous run
      for (uint8_t j=0; j<g->len_nb[iNode]; j++){
        int64_t idx_nb = g->nn[iNode*g->num_nb_max + j];
        if(g->len_nb[idx_nb] > 2 && g->ptr[idx_nb] > MEASUREZ && fusionS[get_edge_idx(g, iNode, idx_nb)]){
          int64_t r2 = findroot(g, idx_nb); // find root of neighbor and do path compression
          if (r2 != r1) r1 = merge_root(g, r1, r2, sizeMax + iNode_rn); // update sizeMax and root (if connected to same root already, do nothing)
        }
      }
    }
    /* find node within largest component of block with highest vertex degree */
    uint8_t maxDegreeBlk = 0; // largest vertex degree in block
    for(int64_t iNode_blk=0; iNode_blk < nnodeBlk; iNode_blk++){
      index_to_digital_vec(iNode_blk, i_vec_blk, n_dim, lblock+1);
      int64_t iNode = get_global_idx(i_vec_rn, i_vec_blk, n_dim, lblock, lTotal);
      if(-g->ptr[findroot(g, iNode)] != sizeMax[iNode_rn]) continue; // only consider nodes inside the largest connected component
      uint8_t degree = 0;
      for (uint8_t j=0; j<g->len_nb[iNode]; j++){
        int64_t idx_nb = g->nn[iNode*g->num_nb_max + j];
        if(g->len_nb[idx_nb] > 2 && g->ptr[idx_nb] > MEASUREZ && fusionS[get_edge_idx(g, iNode, idx_nb)]) degree++; // do not count edges at block surface (>2)
      }
      if(degree > maxDegreeBlk){
        maxDegreeBlk = degree;
        (*bestNode)[iNode_rn] = iNode;
        rn.static_node[iNode_rn] = g->static_node[iNode];
      }
    }
  }

  /* loop blocks and check if connection through common face to largest connected component of neighboring block is possible, if yes: add bond in renormalized graph */
  int numAllPossibleConnection = 0;
  int64_t nFace = 1; // number of nodes on a face
  for(uint8_t i=0; i<n_dim-1; i++) nFace *= lblock;
  for(int64_t iNode_rn=0; iNode_rn < rn.nnode; iNode_rn++){
    index_to_digital_vec(iNode_rn, i_vec_rn, n_dim, lsize);
    for(uint8_t d=0;d<n_dim;d++){
      if(periodic || i_vec_rn[d] < lsize-1){ // check connection to other blocks in forward direction
        numAllPossibleConnection++;
        for(int64_t i_face=0; i_face<nFace; i_face++){
          /* node inside block */
          index_to_digital_vec_skip(i_face, i_vec_blk, n_dim, d, lblock-1, lblock);
          int64_t iNodeIn = get_global_idx(i_vec_rn, i_vec_blk, n_dim, lblock, lTotal);
          if(-g->ptr[findroot(g, iNodeIn)] != sizeMax[iNode_rn]) continue; // no connection if node is not part of the largest component of block (smaller component, LOST, or MEASUREZ)
          /* node on other side of face (neighboring block) */
          i_vec_blk[d] = 0;
          int64_t save = i_vec_rn[d];
          i_vec_rn[d] = (i_vec_rn[d] + 1) % lsize;
          int64_t iNode_rn_nb = digital_vec_to_index(i_vec_rn, n_dim, lsize); // index of neighboring block
          int64_t iNodeOut = get_global_idx(i_vec_rn, i_vec_blk, n_dim, lblock, lTotal);
          i_vec_rn[d] = save;
          if(-g->ptr[findroot(g, iNodeOut)] != sizeMax[iNode_rn_nb]) continue; // no connection if node is not part of the largest component of neighboring block (smaller component, LOST, or MEASUREZ)
          /* node at face */
          i_vec_blk[d] = lblock;
          int64_t iNodeFace = get_global_idx(i_vec_rn, i_vec_blk, n_dim, lblock, lTotal);
          /* create connection */
          if(g->ptr[iNodeFace] == -1 && fusionS[get_edge_idx(g, iNodeFace, iNodeIn)] && fusionS[get_edge_idx(g, iNodeFace, iNodeOut)]){
            rn.num_edges++;
            rn.nn[iNode_rn*rn.num_nb_max + rn.len_nb[iNode_rn]] = iNode_rn_nb;
            rn.len_nb[iNode_rn]++;
            rn.nn[iNode_rn_nb*rn.num_nb_max + rn.len_nb[iNode_rn_nb]] = iNode_rn;
            rn.len_nb[iNode_rn_nb]++;
            break;
          }
        }
      }
    }
  }
  *bondFraction =  (double)rn.num_edges / (double)numAllPossibleConnection;
  free(fusionS);
  free(sizeMax);
  return rn;
}

/* create Raussendorf lattice (lsize: number of 6-qubit unit cell in one dimension) */
Graph get_raussendorf_lattice(int64_t lsize, bool static_center, bool periodic, bool get_size, bool edge_list){
  uint8_t n_dim = 3;
  int64_t nUnit = 2*n_dim; // size of the ring graph that makes the unit cell
  int64_t nCell = 1; // number of all unit cells
  uint8_t num_nb_max = 4; // number of neighbors
  for(uint8_t j=0;j<n_dim;j++) nCell *= lsize; // number of nodes
  Graph g = new_graph(nCell*nUnit, num_nb_max, get_size, edge_list, false);
  int64_t i_vec[n_dim]; // digital system representation (of unit cell index)
  memset(g.len_nb, 0, nCell*nUnit);
  for(int64_t i_cell=0;i_cell<nCell;i_cell++){
    index_to_digital_vec(i_cell, i_vec, n_dim, lsize); // get digital system representation (i_vec) of the unit cell
    for(int i_node=0;i_node<2*n_dim;i_node++){ // i_node is the local index in the cell
      int64_t idxNode1 = i_cell * nUnit + i_node; // global index of first node
      int64_t idxNode2 = i_cell * nUnit + (i_node + 1) % (2*n_dim); // global index of second node
      g.nn[idxNode1 * num_nb_max + (g.len_nb[idxNode1]++)] = idxNode2;
      g.nn[idxNode2 * num_nb_max + (g.len_nb[idxNode2]++)] = idxNode1;
      bool sign = i_node % 2; // 0: going from face qubit to edge, 1: going from edge qubit to face (in the unit cell ring)
      uint8_t d = i_node % n_dim; // dimension into which the connection in the ring unit cell goes
      if(sign == 0 && (periodic || i_vec[d] < lsize - 1)){
        int64_t save = i_vec[d];
        i_vec[d] = (i_vec[d] + 1) % lsize;
        int64_t i_cell_nb = digital_vec_to_index(i_vec, n_dim, lsize);
        i_vec[d] = save;
        idxNode2 = i_cell_nb * nUnit + (i_node + 1) % (2*n_dim); // global index of second node
        g.nn[idxNode1 * num_nb_max + (g.len_nb[idxNode1]++)] = idxNode2;
        g.nn[idxNode2 * num_nb_max + (g.len_nb[idxNode2]++)] = idxNode1;
      } else if(sign == 1 && (periodic || i_vec[d] > 0)){
        int64_t save = i_vec[d];
        i_vec[d] = (i_vec[d] > 0 ? i_vec[d] - 1 : lsize - 1);
        int64_t i_cell_nb = digital_vec_to_index(i_vec, n_dim, lsize);
        i_vec[d] = save;
        idxNode2 = i_cell_nb * nUnit + (i_node + 1) % (2*n_dim); // global index of second node
        g.nn[idxNode1 * num_nb_max + (g.len_nb[idxNode1]++)] = idxNode2;
        g.nn[idxNode2 * num_nb_max + (g.len_nb[idxNode2]++)] = idxNode1;
      }
    }
  }
  if(edge_list) make_edge_array(&g);
  if(static_center){
    memset(g.static_node, 1, sizeof(bool)*g.nnode);
  } else{
    memset(g.static_node, 0, sizeof(bool)*g.nnode);
  }
  if(!get_size){
    memset(g.start,0,g.nnode*sizeof(bool));
    memset(g.stop,0,g.nnode*sizeof(bool));
    for (int64_t i_cell=0; i_cell<nCell/lsize; i_cell++) memset(g.start + i_cell*nUnit, 1, nUnit); // make all nodes cells at edge start nodes
    for (int64_t i_cell=nCell/lsize*(lsize-1); i_cell<nCell; i_cell++) memset(g.stop + i_cell*nUnit, 1, nUnit); // make all nodes cells at edge stop nodes
  }
  return g;
}

/* helper function for get_nd_simple_cubic_modified (note: can overflow for very high dimensions (n_dim)):
   get the number of all diagonal connections with step number per diagonal between [num_diag_min, num_diag_max] */
static uint64_t get_additional_edge_num(uint8_t n_dim, uint8_t num_diag_min, uint8_t num_diag_max){
  uint64_t sum = 0;
  uint64_t bincoeff;
  uint64_t pow2;
  for(uint8_t n=num_diag_min; n<=num_diag_max; n++){
   bincoeff = 1;
   pow2 = 1;
   for(uint8_t k=0; k<n; k++){
     bincoeff = bincoeff * (n_dim - k);
     pow2 *= 2;
   }
   for(uint8_t k=2; k<=n; k++) bincoeff = bincoeff / k;
   sum = sum + pow2 * bincoeff;
  }
  return sum;
}

/* create an n-dimensional simple cubic lattice with diagonal connections (above 7 dimensions, num_nb_max overflows)
   lsize, n_dim: lattice size in one dimension, number of dimensions 
   static_center: central qubits cannot be lost if true
   vertical: turn on/off simple cubic connections
   bcc: turn on/off bcc connections
   diamond: leave out connections for diamond brickwork lattice
   num_diag_min, num_diag_max: specify range of more general lattice connection vectors that are considered. Both values refer to the number of non-zero vector elements
   periodic: periodic boundaries if true
   get_size: get size of largest connected component if true, check for cluster spanning otherwise
   edge_list: get list of edges for e.g. bond percolation if true */
Graph get_nd_simple_cubic_modified(int64_t lsize, uint8_t n_dim, bool static_center, bool vertical, bool bcc, bool diamond, uint8_t num_diag_min, uint8_t num_diag_max, bool periodic, bool get_size, bool edge_list){
  int64_t numNodes = 1; // count the nodes
  uint8_t num_bcc = 1; // maximum number of bcc diagonals
  uint32_t num_all = 1; // count all vectors in {-1,0,+1}^n_dim
  for(uint8_t j=0;j<n_dim;j++){
    numNodes *= lsize; // number of nodes
    if(bcc) num_bcc *= 2; // count the bcc diagonal connections
    num_all *= 3; // -1,0,+1 per dimension
  }
  if(num_diag_min < 2) num_diag_min = 2; // avoid that simple cubic connections are added twice
  uint8_t num_nb_max = 0; // maximum number of neighbors
  if(vertical) num_nb_max += 2*n_dim; // count simple cubic connections
  if(bcc && num_diag_max < n_dim) num_nb_max += num_bcc; // count bcc diagonal connections (if not part of general diagonals)
  if(num_diag_max > 1) num_nb_max += get_additional_edge_num(n_dim, num_diag_min, num_diag_max); // count all other diagonal edges

  Graph g = new_graph(numNodes, num_nb_max, get_size, edge_list, false);

  if(diamond && periodic && lsize % 2 == 1){
    printf("%s\n", "leave out every second edge, periodic, odd lsize: would lead to unidirectional edges at the boundaries");
    g.failed = true; // indicates that the created graph is not valid
    return g; // no graph can be created
  }

  int64_t i_vec[n_dim]; // digital system representation (of node index)
  uint8_t digit_pos, nbIdx;
  int64_t iNb; // index of neighbor
  for(int64_t i_node=0;i_node<numNodes;i_node++){
    index_to_digital_vec(i_node, i_vec, n_dim, lsize); // get digital system representation (i_vec) of i_node
    g.static_node[i_node] = (static_center ? true : false);
    nbIdx = 0;
    int64_t sum_i_vec = 0;
    for(int i=0; i<n_dim; i++) sum_i_vec += i_vec[i];

    // add simple cubic connections
    if(vertical){
      for(digit_pos=0;digit_pos<n_dim;digit_pos++){ // loop all simple cubic connections
        // step in positive direction
        if(diamond==false || digit_pos==0 || (sum_i_vec % 2 == 0)){ // if diamond=true, there is just one dimension with steps in both directions
          if(i_vec[digit_pos] < lsize - 1){
            i_vec[digit_pos] += 1;
            iNb = digital_vec_to_index(i_vec, n_dim, lsize);
            i_vec[digit_pos] -= 1;
            g.nn[i_node * num_nb_max + (nbIdx++)] = iNb; // connect to neighbor
            if(edge_list) append_edge(&g, i_node, nbIdx-1);
          } else if (periodic) {
            int64_t save = i_vec[digit_pos];
            i_vec[digit_pos] = 0;
            iNb = digital_vec_to_index(i_vec, n_dim, lsize);
            i_vec[digit_pos] = save;
            g.nn[i_node * num_nb_max + (nbIdx++)] = iNb; // connect to neighbor
            if(edge_list) append_edge(&g, i_node, nbIdx-1);
          }
        }
        // step in negative direction
        if(diamond==false || digit_pos==0 || (sum_i_vec % 2 == 1)){ // if diamond=true, there is just one dimension with steps in both directions
          if(i_vec[digit_pos] > 0){
            i_vec[digit_pos] -= 1;
            iNb = digital_vec_to_index(i_vec, n_dim, lsize);
            i_vec[digit_pos] += 1;
            g.nn[i_node * num_nb_max + (nbIdx++)] = iNb; // connect to neighbor
            if(edge_list) append_edge(&g, i_node, nbIdx-1);
          } else if (periodic) {
            int64_t save = i_vec[digit_pos];
            i_vec[digit_pos] = lsize-1;
            iNb = digital_vec_to_index(i_vec, n_dim, lsize);
            i_vec[digit_pos] = save;
            g.nn[i_node * num_nb_max + (nbIdx++)] = iNb; // connect to neighbor
            if(edge_list) append_edge(&g, i_node, nbIdx-1);
          }
        }
      }
    }

    // add general diagonal connections
    if(num_diag_max > 0){
      int64_t vec_diagonal[n_dim];
      for(uint64_t i=0;i<num_all; i++){ // loop all diagonal connections
        index_to_digital_vec(i_node, i_vec, n_dim, lsize); // get digital system representation (i_vec) of i_node
        index_to_digital_vec(i, vec_diagonal, n_dim, 3); // get binary representation (vec_diagonal) of the diagonal (representation: 0:-1, 1:0, 2:+1)
        int sum = 0;
        for(uint8_t k=0;k<n_dim;k++){
          if(vec_diagonal[k]!=1) sum = sum + 1; // !=1 means that there is a +- step in that dimension
        }
        if(sum < num_diag_min || sum > num_diag_max) continue; // the number of steps corresponds to a diagonal that is not considered
        bool stop = false;
        for(digit_pos=0; digit_pos<n_dim; digit_pos++){
          i_vec[digit_pos] += (vec_diagonal[digit_pos] - 1); // -1, 0, +1
          if(i_vec[digit_pos] == -1){
            if(periodic){i_vec[digit_pos] = lsize-1;}
            else{stop=true; break;}
          } else if (i_vec[digit_pos] == lsize){
            if(periodic){i_vec[digit_pos] = 0;}
            else{stop=true; break;}
          }
        }
        if(stop == false){
          iNb = digital_vec_to_index(i_vec, n_dim, lsize); // compute neighbor index
          g.nn[i_node * num_nb_max + (nbIdx++)] = iNb; // connect to neighbor
          if(edge_list) append_edge(&g, i_node, nbIdx-1);
        }
      }
    }

    // add bcc diagonal connections (can be done with general diagonal connections, but this here is faster)
    if(bcc && num_diag_max < n_dim){ // 2nd check to avoid making edges twice
      int64_t vec_diagonal[n_dim];
      for(int i=0;i<num_bcc; i++){ // loop all bcc diagonal connections
        index_to_digital_vec(i_node, i_vec, n_dim, lsize); // get digital system representation (i_vec) of i_node
        index_to_digital_vec(i, vec_diagonal, n_dim, 2); // get binary representation (vec_diagonal) of the bcc diagonal (representation: 0:-1, 1:+1)
        bool stop = false;
        for(digit_pos=0; digit_pos<n_dim; digit_pos++){
          i_vec[digit_pos] += (vec_diagonal[digit_pos] ? 1 : -1);
          if(i_vec[digit_pos] == -1){
            if(periodic){i_vec[digit_pos] = lsize-1;}
            else{stop=true; break;}
          } else if (i_vec[digit_pos] == lsize){
            if(periodic){i_vec[digit_pos] = 0;}
            else{stop=true; break;}
          }
        }
        if(stop == false){
          iNb = digital_vec_to_index(i_vec, n_dim, lsize); // compute neighbor index
          g.nn[i_node * num_nb_max + (nbIdx++)] = iNb; // connect to neighbor
          if(edge_list) append_edge(&g, i_node, nbIdx-1);
        }   
      }
    }

    g.len_nb[i_node] = nbIdx; // number of neighbors
  }
  if(!get_size){
    memset(g.start,0,g.nnode*sizeof(bool));
    memset(g.stop,0,g.nnode*sizeof(bool));
    for (int64_t i_node=0; i_node<numNodes/lsize; i_node++) g.start[i_node] = 1;
    for (int64_t i_node=numNodes/lsize*(lsize-1); i_node<numNodes; i_node++) g.stop[i_node] = 1;
  }
  return g;
}

/* helper function: get linearly independent vectors pointing from one hypercubic lattice point to the next (vectors get shuffled)
   n_dim: number of dimensions
   max_step: maximum number of integer steps per coordinate */
int64_t* get_non_equivalent_vecs(uint8_t n_dim, uint8_t max_step){
  int64_t num_vec = 1;
  for(uint8_t j=0;j<n_dim;j++) num_vec = num_vec * (2*(int64_t)max_step + 1);
  num_vec -= 1; // -1 because 0-vector is not considered
  num_vec /= 2; // if two vectors just differ by sign, just look at one of them
  int64_t* vecs = malloc(sizeof(int64_t) * (int64_t)n_dim * num_vec); // every block of size n_dim is one vector

  int64_t pos = 0; // counts number of vectors added
  int64_t vec[n_dim];
  for(int64_t i=0;i<num_vec*2 + 1; i++){
    index_to_digital_vec(i, vec, n_dim, 2*(int64_t)max_step + 1); // elements of vec are in range [0, 2*max_step]
    for(uint8_t j=0;j<n_dim;j++){
      if(vec[j] - (int64_t)max_step > 0){
        for(uint8_t k=0;k<n_dim;k++) vecs[pos * (int64_t)n_dim + (int64_t)k] = vec[k] - (int64_t)max_step;
        pos++;
        break;
      } else if(vec[j] - (int64_t)max_step < 0){
        break;
      }
    }
  }
  // shuffle vectors
  for(int64_t i=0;i<num_vec; i++){
    pos = i + (int64_t)round((num_vec-1-i) * (double)rand() / (double)RAND_MAX);
    memcpy(vec, vecs + i * (int64_t)n_dim, sizeof(int64_t) * (int64_t)n_dim);
    memcpy(vecs + i * (int64_t)n_dim, vecs + pos * (int64_t)n_dim, sizeof(int64_t) * (int64_t)n_dim);
    memcpy(vecs + pos * (int64_t)n_dim, vec, sizeof(int64_t) * (int64_t)n_dim);
  }
  return vecs;
}

/* construct a lattice starting from vertices on a hypercubic lattice and integer vectors connecting each vertex to its neighbors */
Graph get_lattice_from_nd_simple_cubic_and_vectors(int64_t lsize, uint8_t n_dim, int64_t* vecs, int64_t num_vecs, bool static_center, bool periodic, bool get_size, bool edge_list){
  int64_t numNodes = 1;
  for(uint8_t j=0;j<n_dim;j++) numNodes *= lsize; // number of nodes
  if(2 * num_vecs > 256 - 1) printf("%s\n", "error: too many neighbors");
  uint8_t num_nb_max = 2 * num_vecs; // maximum number of neighbors
  Graph g = new_graph(numNodes, num_nb_max, get_size, edge_list, false);

  int64_t i_vec[n_dim]; // digital system representation (of node index)
  for(int64_t i_node=0;i_node < numNodes;i_node++){
    g.static_node[i_node] = (static_center ? true : false);
    uint8_t nbIdx = 0;
    for(int64_t n_vec=0; n_vec < num_vecs; n_vec++){ // loop over all vectors connecting to neighbors
      for(uint8_t direction=0; direction<2; direction++){ // look at vector with two different signs
        index_to_digital_vec(i_node, i_vec, n_dim, lsize); // get/reset digital system representation i_vec for i_node
        bool stop = false; // stops when it becomes clear that i_vec would be outside the lattice
        for(int digit_pos=0; digit_pos<n_dim; digit_pos++){
          i_vec[digit_pos] += (1 - 2*direction) * vecs[n_vec * (int64_t)n_dim + digit_pos];
          if(i_vec[digit_pos] < 0){
            if(periodic){i_vec[digit_pos] = lsize + i_vec[digit_pos];}
            else{stop=true; break;}
          } else if (i_vec[digit_pos] > lsize - 1){
            if(periodic){i_vec[digit_pos] = i_vec[digit_pos] % lsize;}
            else{stop=true; break;}
          }
        }
        if(stop == false){
          g.nn[i_node * num_nb_max + (nbIdx++)] = digital_vec_to_index(i_vec, n_dim, lsize);  // compute neighbor index with digital_vec_to_index
          if(edge_list) append_edge(&g, i_node, nbIdx-1);
        }
      }
    }
    g.len_nb[i_node] = nbIdx; // number of neighbors
  }

  if(!get_size){
    memset(g.start,0,g.nnode*sizeof(bool));
    memset(g.stop,0,g.nnode*sizeof(bool));
    for (int64_t i_node=0; i_node<numNodes/lsize; i_node++) g.start[i_node] = 1;
    for (int64_t i_node=numNodes/lsize*(lsize-1); i_node<numNodes; i_node++) g.stop[i_node] = 1;
  }
  return g;
}

/* detached micro-clusters to create an n-dimensional simple cubic lattice (periodic boundary conditions) */
Graph get_GHZ_stars_for_nd_simple_cubic(int64_t lsize, uint8_t n_dim, bool static_center, bool periodic, bool get_size){
  int64_t num_cluster = 1;
  int64_t num_mu = 2*n_dim + 1; // size of micro-cluster: 2 neighbors in star per dimension plus central node
  uint8_t num_nb_max = 2*n_dim; // maximum number of neighbors
  for(uint8_t j=0;j<n_dim;j++) num_cluster *= lsize; // number of star micro-clusters
  Graph g = new_graph(num_cluster * num_mu, num_nb_max, get_size, false, true);
  int64_t i_vec_cluster[n_dim]; // digital system representation (of cluster index)
  uint8_t digit_pos;
  int64_t i_nb_cluster, i_node; // index of neighboring cluster, node index
  for(int64_t i_cluster=0;i_cluster<num_cluster;i_cluster++){
    index_to_digital_vec(i_cluster, i_vec_cluster, n_dim, lsize); // get digital system representation (i_vec_cluster) of i_cluster

    // central node
    i_node = i_cluster * num_mu;
    g.static_node[i_node] = (static_center ? true : false);
    g.len_nb[i_node] = num_nb_max; // maximum  number of neighbors for central node
    g.fusion_node[i_node] = false;
    for(uint8_t k=0;k<num_nb_max;k++) g.nn[i_node * num_nb_max + k] = i_node + k + 1; // connect to all leaf-nodes of star

    // fusion nodes: connect to neighboring micro-clusters, looping over all dimensions
    for(digit_pos=0;digit_pos<n_dim;digit_pos++){
      i_node = i_cluster * num_mu + 1 + digit_pos*2 + 1;
      if(i_vec_cluster[digit_pos] < lsize - 1){
        i_vec_cluster[digit_pos] += 1; // go forward in this dimension
        i_nb_cluster = digital_vec_to_index(i_vec_cluster, n_dim, lsize);
        i_vec_cluster[digit_pos] -= 1;
        g.fusion_node[i_node] = true;
        g.fusion_partner[i_node] = i_nb_cluster * num_mu + 1 + digit_pos*2 + 0;
      } else if (periodic){ // start from other side --> periodic boundary
        int64_t save = i_vec_cluster[digit_pos];
        i_vec_cluster[digit_pos] = 0; // start from other side
        i_nb_cluster = digital_vec_to_index(i_vec_cluster, n_dim, lsize);
        i_vec_cluster[digit_pos] = save; // reset
        g.fusion_node[i_node] = true;
        g.fusion_partner[i_node] = i_nb_cluster * num_mu + 1 + digit_pos*2 + 0;
      } else {
        g.fusion_node[i_node] = false; // don't treat as fusion node when at edge
      }
      g.static_node[i_node] = false; // make all fusion nodes non-static
      g.nn[i_node * num_nb_max] = i_cluster * num_mu; // connect to central node
      g.len_nb[i_node] = 1; // only one neighbor (fusion node)

      i_node = i_cluster * num_mu + 1 + digit_pos*2 + 0;
      if(i_vec_cluster[digit_pos] > 0){
        i_vec_cluster[digit_pos] -= 1; // go backward in this dimension
        i_nb_cluster = digital_vec_to_index(i_vec_cluster, n_dim, lsize);
        i_vec_cluster[digit_pos] += 1;
        g.fusion_node[i_node] = true;
        g.fusion_partner[i_node] = i_nb_cluster * num_mu + 1 + digit_pos*2 + 1;
      } else if (periodic){ // start from other side --> periodic boundary
        int64_t save = i_vec_cluster[digit_pos];
        i_vec_cluster[digit_pos] = lsize-1; // start from other side
        i_nb_cluster = digital_vec_to_index(i_vec_cluster, n_dim, lsize);
        i_vec_cluster[digit_pos] = save; // reset
        g.fusion_node[i_node] = true;
        g.fusion_partner[i_node] = i_nb_cluster * num_mu + 1 + digit_pos*2 + 1;
      } else {
        g.fusion_node[i_node] = false; // don't treat as fusion node when at edge
      }
      g.static_node[i_node] = false; // make all fusion nodes non-static
      g.nn[i_node * num_nb_max] = i_cluster * num_mu; // connect to central node
      g.len_nb[i_node] = 1; // only one neighbor (fusion node)
    }
  }
  if(!get_size){
    memset(g.start,0,g.nnode*sizeof(bool));
    memset(g.stop,0,g.nnode*sizeof(bool));
    for (int64_t i_cluster=0; i_cluster<num_cluster/lsize; i_cluster++) g.start[i_cluster * num_mu] = 1;  // only set central node
    for (int64_t i_cluster=num_cluster/lsize*(lsize-1); i_cluster<num_cluster; i_cluster++) g.stop[i_cluster * num_mu] = 1;
  }
  return g;
}

/* detached 5-qubit star custers to aim for construction of a 2d square lattice */
Graph get_5qubit_stars_for_2d_square(int64_t lsize, bool static_center, bool periodic, bool get_size){
  uint8_t num_nb_max = 4; // maximum number of neighbors
  int64_t num_mu = 5; // size of micro-cluster
  int64_t num_cluster = lsize*lsize;
  Graph g = new_graph(num_cluster * num_mu, num_nb_max, get_size, false, true); // 4 qubits for fusion per site, + 1 qubits in center
  for(int64_t i=0; i<num_cluster; i++){ // loop over all micro-clusters
    for(uint8_t k=0;k<num_nb_max;k++) g.nn[i*num_mu*num_nb_max + k] = i*num_mu + k + 1; // central node
    g.len_nb[i*num_mu + 0] = num_nb_max;
    g.fusion_node[i*num_mu + 0] = false;
    g.static_node[i*num_mu + 0] = (static_center ? true : false);
    for(uint8_t k=1;k<num_mu;k++){ // loop all fusion nodes
      g.nn[(i*num_mu + k) * num_nb_max + 0] = i*num_mu + 0; // fusion-nodes of star-cluster
      g.len_nb[i*num_mu + k] = 1; // only one neighbor (fusion node)
      g.fusion_node[i*num_mu + k] = true;
      g.static_node[i*num_mu + k] = false; // make all fusion nodes non-static
      // specify the fusion partner
      if(k==1){
        if(i>=lsize) g.fusion_partner[i*num_mu + k] = (i-lsize)*num_mu + k + 2;
        else if (periodic) g.fusion_partner[i*num_mu + k] = (i+lsize*(lsize-1))*num_mu + k + 2; // periodic boundaries for fusion (fusion qubits must have defined partner)
        else {g.fusion_node[i*num_mu + k] = false;} // don't treat as fusion node when at edge
      } else if (k==3){
        if(i<lsize*(lsize-1)) g.fusion_partner[i*num_mu + k] = (i+lsize)*num_mu + k - 2;
        else if (periodic) g.fusion_partner[i*num_mu + k] = (i-lsize*(lsize-1))*num_mu + k - 2; // periodic boundaries for fusion (fusion qubits must have defined partner)
        else {g.fusion_node[i*num_mu + k] = false;} // don't treat as fusion node when at edge
      } else if(k==2){
        if(i%lsize < lsize-1) g.fusion_partner[i*num_mu + k] = (i+1)*num_mu + k + 2;
        else if (periodic) g.fusion_partner[i*num_mu + k] = (i-lsize+1)*num_mu + k + 2; // periodic boundaries for fusion (fusion qubits must have defined partner)
        else {g.fusion_node[i*num_mu + k] = false;} // don't treat as fusion node when at edge
      } else if(k==4){
        if(i%lsize > 0) g.fusion_partner[i*num_mu + k] = (i-1)*num_mu + k - 2;
        else if (periodic) g.fusion_partner[i*num_mu + k] = (i+lsize-1)*num_mu + k - 2; // periodic boundaries for fusion (fusion qubits must have defined partner)
        else {g.fusion_node[i*num_mu + k] = false;} // don't treat as fusion node when at edge
      }
    }
  }
  if(!get_size){
    memset(g.start,0,g.nnode*sizeof(bool));
    memset(g.stop,0,g.nnode*sizeof(bool));
    for (int64_t i_cluster=0; i_cluster<num_cluster/lsize; i_cluster++) g.start[i_cluster * num_mu] = 1;  // only set central node
    for (int64_t i_cluster=num_cluster/lsize*(lsize-1); i_cluster<num_cluster; i_cluster++) g.stop[i_cluster * num_mu] = 1;
  }
  return g;
}

/* construction from 10.1038/s41467-019-08948-x to get 2d square lattice from 3photon GHZ */
Graph get_3qubit_GHZ_for_2d_square(int64_t lsize, bool static_center, bool periodic, bool get_size){
  uint8_t num_nb_max = 2; // maximum number of neighbors
  int64_t num_mu = 9; // size of micro-cluster unit cell
  int64_t num_cluster = lsize*lsize;
  Graph g = new_graph(num_cluster * num_mu, num_nb_max, get_size, false, true);
  int64_t x, y, k, idx_node;
  for(x=0;x<lsize;x++){
    for(y=0;y<lsize;y++){
      for(k=0;k<num_mu;k++){
        idx_node = (x*lsize + y)*num_mu + k; // index of current node
        g.static_node[idx_node] = false; // node = photon (can be lost) as default
        if(k==1){
          g.len_nb[idx_node] = 2; // central GHZ, central node
          g.fusion_node[idx_node] = false; // central node
          g.static_node[idx_node] = (static_center ? true : false);
          g.nn[idx_node*num_nb_max + 0] = idx_node + 4; g.nn[idx_node*num_nb_max + 1] = idx_node + 5;
        } else if (k<3){
          g.len_nb[idx_node] = 2; // two neighbors (center of GHZ)
          g.fusion_node[idx_node] = true; // used for fusion
          g.nn[idx_node*num_nb_max + 0] = idx_node + 3 + k; g.nn[idx_node*num_nb_max + 1] = idx_node + 4 + k;
          g.fusion_partner[idx_node] = (k==0) ? (idx_node + 5) : (idx_node + 4);
        } else {
          g.len_nb[idx_node] = 1; // edge of GHZ
          g.fusion_node[idx_node] = true; // used for fusion
          if(k==3 || k==4){g.nn[idx_node*num_nb_max + 0] = idx_node - k;}
          else if(k==5 || k==6){g.nn[idx_node*num_nb_max + 0] = idx_node - k + 1;}
          else{g.nn[idx_node*num_nb_max + 0] = idx_node - k + 2;}
          if(k==5){g.fusion_partner[idx_node] = idx_node-5;}
          else if(k==6){g.fusion_partner[idx_node] = idx_node-4;}
          else if(k==3){
            if(x>0) g.fusion_partner[idx_node] = idx_node - lsize*num_mu +1; // number 4 from neighbor unit cell
            else if (periodic) g.fusion_partner[idx_node] = idx_node + lsize*(lsize-1)*num_mu +1; // periodic boundaries for fusion (fusion qubits must have defined partner)
            else {g.fusion_node[idx_node] = false;} // don't treat as fusion node when at edge
          }
          else if(k==4){
            if(x<lsize-1) g.fusion_partner[idx_node] = idx_node + lsize*num_mu -1; // number 3 from neighbor unit cell
            else if (periodic) g.fusion_partner[idx_node] = idx_node - lsize*(lsize-1)*num_mu -1; // periodic boundaries for fusion (fusion qubits must have defined partner)
            else {g.fusion_node[idx_node] = false;} // don't treat as fusion node when at edge
          }
          else if(k==7){
            if(y>0) g.fusion_partner[idx_node] = idx_node - num_mu +1; // number 8 from neighbor unit cell
            else if (periodic) g.fusion_partner[idx_node] = idx_node + (lsize-1)*num_mu +1; // periodic boundaries for fusion (fusion qubits must have defined partner)
            else {g.fusion_node[idx_node] = false;} // don't treat as fusion node when at edge
          }
          else if(k==8){
            if(y<lsize-1) g.fusion_partner[idx_node] = idx_node + num_mu -1; // number 7 from neighbor unit cell
            else if (periodic) g.fusion_partner[idx_node] = idx_node - (lsize-1)*num_mu -1; // periodic boundaries for fusion (fusion qubits must have defined partner)
            else {g.fusion_node[idx_node] = false;} // don't treat as fusion node when at edge
          }
        }
      }
    }
  }
  if(!get_size){
    memset(g.start,0,g.nnode*sizeof(bool));
    memset(g.stop,0,g.nnode*sizeof(bool));
    for (int64_t i_cluster=0; i_cluster<num_cluster/lsize; i_cluster++) g.start[i_cluster * num_mu + 1] = 1;  // only set central node
    for (int64_t i_cluster=num_cluster/lsize*(lsize-1); i_cluster<num_cluster; i_cluster++) g.stop[i_cluster * num_mu + 1] = 1;
  }
  return g;
}
