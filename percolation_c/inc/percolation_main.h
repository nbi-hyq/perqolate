#ifndef PERCOLATION
#define PERCOLATION

int64_t findroot(Graph* g, int64_t i);
int64_t merge_root(Graph* g, int64_t r1, int64_t r2, int64_t* big);
int64_t get_edge_idx(Graph* g, int64_t s1, int64_t s2);
uint64_t rand_range(uint64_t n);
int64_t* permutation(int64_t len);
int64_t* permutation_static_first(int64_t len, bool staticFirst, bool* arryStatic, int64_t* numStatic);
int64_t bfs_full(Graph* g);
int64_t* percolate_site(Graph* g, float pBond, int64_t* idxLambda);
int64_t* percolate_bond(Graph* g, float pSite, int64_t* idxLambda);
int64_t fusion_and_loss_no_newman_ziff(Graph* g, float pFusion, float pNoLoss);
int64_t* fusion_and_loss_nz(Graph* g, float pFusion, bool considerStaticNode, int64_t* numQbts, int64_t* numStaticNode, int64_t* idxLambda);
int64_t* fusion_repeat_and_loss_nz(Graph* g, float pFusion, uint8_t n_fusion, int64_t* numQbts, int64_t* numStaticNode, int64_t* idxLambda);
int64_t analyse_fusion_net(Graph* g, float pSite, uint8_t nFusionPhotons, float pFail);
bool* apply_loss_to_fusion_net(Graph* g, float pSite, uint8_t nFusionPhotons, float pFail);
int64_t apply_loss(Graph* g, float pSite);
int64_t apply_loss_bfs(Graph* g, float pSite);
int64_t* apply_loss_nz(Graph* g, float pBond, bool considerStaticNode, int64_t* numStaticNode, int64_t* idxLambda);
float* get_expectation_value(float* sweep_prob, int sweep_len, int64_t* largest_component, int64_t N, int64_t nStatic, int64_t size_max, bool approx, bool mute);

#endif
