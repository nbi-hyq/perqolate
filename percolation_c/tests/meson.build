tests = [
  'fusion_nD_simple_cubic',
  'fusion_GHZ_for_2d_cluster_square',
  'fusion_nz_vs_apply_loss',
  'fusion_repeat_nD_simple_cubic',
  'site_2D_square_lattice',
  'block_renormalize',
  'bond_nd_simple_cubic',
  'bond_nd_simple_cubic_block',
  'bond_raussendorf',
  'bond_triangular_2D_lattice',
  'effective_bond_nD_simple_cubic',
  'loss_2d_triangular',
  'loss_nz_2d_triangular',
  'percol_fusion_nD_simple_cubic',
  'percol_fusion_repeat_nD_simple_cubic',
  'percol_fusion_GHZ_for_2d_cluster_square',
  'percol_site_2D_square_lattice',
  'percol_bond_nD_sc',
  'percol_bond_nD_bcc',
  'percol_bond_nD_fcc',
  'percol_bond_nD_diamond',
  'percol_bond_triangular_2D_lattice',
  'percol_loss_2d_triangular',
  'percol_loss_nz_2d_triangular',
  'tree_lattice_site',
  'complex_neighbor_reference_NN_4NN',
  'complex_neighbor_reference_2NN_4NN',
  'complex_neighbor_reference_3NN_4NN',
  'complex_neighbor_reference_NN_2NN_3NN',
  'check_permutation',
  'rand_range_a',
  'rand_range_b',
  'rand_lib',
]

foreach t : tests
  test_exe = executable(
    'test_' + t,
    [t + '.c', '../src/percolation_main.c', '../src/graph_type.c', '../src/graph_construction.c'],
    include_directories : ['../inc/'],
  )
  test(t, test_exe)
endforeach
