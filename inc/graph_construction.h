#ifndef GRAPH_CONSTRUCT
#define GRAPH_CONSTRUCT

Graph get_2d_square(int64_t lsize, bool periodic, bool get_size, bool edge_list);
Graph get_2d_triangular(int64_t lsize, bool periodic, bool get_size, bool edge_list);
Graph get_tree_lattice(uint8_t depth_max, uint8_t* num_child, bool static_center, bool get_size, bool edge_list);
Graph get_nd_simple_cubic(int64_t lsize, uint8_t n_dim, bool static_center, bool periodic, bool get_size, bool edge_list);
Graph get_nd_simple_cubic_block(int64_t lsize, int lblock, uint8_t n_dim, bool static_center, bool periodic, bool get_size, bool edge_list);
Graph get_raussendorf_lattice(int64_t lsize, bool static_center, bool periodic, bool get_size, bool edge_list);
Graph get_nd_simple_cubic_modified(int64_t lsize, uint8_t n_dim, bool static_center, bool vertical, bool diagonal, bool diamond, uint8_t num_diag_min, uint8_t num_diag_max, bool periodic, bool get_size, bool edge_list);
Graph get_GHZ_stars_for_nd_simple_cubic(int64_t lsize, uint8_t n_dim, bool static_center, bool periodic, bool get_size);
Graph get_5qubit_stars_for_2d_square(int64_t lsize, bool static_center, bool periodic, bool get_size);
Graph get_3qubit_GHZ_for_2d_square(int64_t lsize, bool static_center, bool periodic, bool get_size);
Graph get_lattice_from_nd_simple_cubic_and_vectors(int64_t lsize, uint8_t n_dim, int64_t* vecs, int64_t num_vecs, bool static_center, bool periodic, bool get_size, bool edge_list);
Graph get_lattice_from_unit_graph(UnitGraph* u, int64_t lsize, uint8_t n_dim, bool static_center, bool periodic, bool get_size, bool edge_list);
Graph renormalize_block_fusion_net(Graph* g, float pSite, uint8_t nFusionPhotons, float pFail, int lblock, int64_t lsize, uint8_t n_dim, bool periodic, bool get_size_rn, double* bondFraction, int64_t** bestNode);
int64_t* get_non_equivalent_vecs(uint8_t n_dim, uint8_t max_step);

#endif
