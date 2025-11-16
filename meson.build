project('percolation_c', 'c',
  version : '0.1',
  default_options : ['warning_level=3'])

m_dep = meson.get_compiler('c').find_library('m', required : true)
if m_dep.found()
  add_project_link_arguments('-lm', language : 'c')
endif

examples = [
  'example_site',
  'example_site_sweep',
  'example_bond',
  'example_bond_raussendorf',
  'example_bond_sweep',
  'example_loss_nz',
  'example_fusion_nz',
  'example_fusion_repeat_sweep_nz',
  'example_fusion_sweep_nz_and_no_nz',
  'example_fusion_sweep_nz',
  'lattice_optimization_add',
  'lattice_optimization_rm',
  'percol_example_fusion_cmp_nz',
  'percol_example_fusion_nz',
  'percol_example_loss_nz',
  'percol_example_bond',
  'percol_example_site',
  'complex_neighborhoods',
  'unit_graph',
  'unit_graph_load',
]

more_examples = [
  'percol_example_tree',
  'example_block_renormalize',
  'example_fusion',
  'example_fusion_sweep',
  'example_loss',
  'example_loss_sweep',
  'percol_example_fusion',
  'percol_example_loss',
]

foreach t : examples
  example_exe = executable(t, ['examples/' + t + '.c', 'src/percolation_main.c', 'src/graph_type.c', 'src/graph_construction.c'])
endforeach

foreach t : more_examples
  example_exe = executable(t, ['examples/more_examples/' + t + '.c', 'src/percolation_main.c', 'src/graph_type.c', 'src/graph_construction.c'])
endforeach

if get_option('tests') and not meson.is_cross_build()
  subdir('tests')
endif
