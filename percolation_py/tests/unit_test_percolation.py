import unittest
import sys
import numpy as np
import time
sys.path.insert(0, '../quantum_states')  # add directory above for import
from graph_states import GraphState
sys.path.insert(0, '../percolation/percolation_py')  # add directory above for import
from graph_construction import GraphPlusConstruction, index_to_digital_vec, digital_vec_to_index

np.random.seed(172473)  # fix the seed

class TestPercolationGraph(unittest.TestCase):
    def test_index_to_digital_vec_and_digital_vec_to_index(self):
        self.assertTrue((np.array(index_to_digital_vec(np.prod([20, 17, 3]) - 1, [20, 17, 3])) == np.array([19, 16, 2])).all())
        self.assertEqual(digital_vec_to_index([1, 2, 1, 0], [4, 4, 3, 2]), 38)
        box = [4, 2, 5, 3, 4, 9]
        for _ in range(100):
            idx = np.random.randint(0, np.prod(box))
            vec = index_to_digital_vec(idx, box)
            idx_new = digital_vec_to_index(vec, box)
            self.assertEqual(idx, idx_new)

    def test_classical_bond_site_percolation(self):
        time_start = time.time()

        lx, ly = 50, 50
        g = GraphPlusConstruction(np.product([lx, ly]))
        g.create_2d_cluster_square(lx, ly, interrupt=(1, 1), shift=(0, 0), static_node=False, get_gcc=False)
        g.set_random_edges_fail(0.4)  # threshold at 0.5 (bond percolation)
        self.assertTrue(g.check_start_target_percolation() == 0)
        g.reset_loss()
        g.set_random_edges_fail(0.6)  # threshold at 0.5 (bond percolation)
        self.assertTrue(g.check_start_target_percolation() == 1)
        g.reset_loss()
        g.set_random_nodes_lost(0.48)  # threshold at 0.5927 (site percolation)
        self.assertTrue(g.check_start_target_percolation() == 0)
        g.reset_loss()
        g.set_random_nodes_lost(0.72)  # threshold at 0.5927 (site percolation)
        self.assertTrue(g.check_start_target_percolation() == 1)
        g.reset_loss()
        g.set_random_nodes_lost(0.75)
        g.set_random_edges_fail(0.60)  # threshold at 0.7262 (site-bond percolation)
        self.assertTrue(g.check_start_target_percolation() == 0)
        g.reset_loss()
        g.set_random_nodes_lost(0.75)
        g.set_random_edges_fail(0.84)  # threshold at 0.7262 (site-bond percolation)
        self.assertTrue(g.check_start_target_percolation() == 1)

        g = GraphPlusConstruction(np.product([lx, ly]))
        g.create_2d_cluster_square(lx, ly, interrupt=(1, 1), shift=(0, 0), static_node=False, get_gcc=True)
        g.set_random_edges_fail(0.4) # threshold at 0.5 (bond percolation)
        g.get_largest_connected_component_size()
        self.assertTrue(g.largest_component_size / g.n_nodes < 0.2)
        g.reset_loss()
        g.set_random_edges_fail(0.6) # threshold at 0.5 (bond percolation)
        g.get_largest_connected_component_size()
        self.assertTrue(g.largest_component_size / g.n_nodes > 0.8)
        g.reset_loss()
        g.set_random_nodes_lost(0.5)  # threshold at 0.5927 (site percolation)
        g.get_largest_connected_component_size()
        self.assertTrue(g.largest_component_size / g.n_nodes < 0.2)
        g.reset_loss()
        g.set_random_nodes_lost(0.7)  # threshold at 0.5927 (site percolation)
        g.get_largest_connected_component_size()
        self.assertTrue(g.largest_component_size / g.n_nodes > 0.5)
        g.reset_loss()
        g.set_random_nodes_lost(0.75)
        g.set_random_edges_fail(0.60)  # threshold at 0.7262 (site-bond percolation)
        g.get_largest_connected_component_size()
        self.assertTrue(g.largest_component_size / g.n_nodes < 0.3)
        g.reset_loss()
        g.set_random_nodes_lost(0.75)
        g.set_random_edges_fail(0.84)  # threshold at 0.7262 (site-bond percolation)
        g.get_largest_connected_component_size()
        self.assertTrue(g.largest_component_size / g.n_nodes > 0.3)

        lx, ly = 50, 50
        g = GraphPlusConstruction(np.product([lx, ly]))
        g.create_2d_cluster_triangular(lx, ly, static_node=False, get_gcc=False)
        g.set_random_edges_fail(0.24)  # threshold at 0.3473 (bond percolation)
        self.assertTrue(g.check_start_target_percolation() == 0)
        g.reset_loss()
        g.set_random_edges_fail(0.44)  # threshold at 0.3473 (bond percolation)
        self.assertTrue(g.check_start_target_percolation() == 1)
        g.reset_loss()
        g.set_random_nodes_lost(0.4)  # threshold at 0.5 (site percolation)
        self.assertTrue(g.check_start_target_percolation() == 0)
        g.reset_loss()
        g.set_random_nodes_lost(0.6)  # threshold at 0.5 (site percolation)
        self.assertTrue(g.check_start_target_percolation() == 1)

        lx, ly = 50, 50
        g = GraphPlusConstruction(np.product([lx, ly]))
        g.create_2d_cluster_triangular(lx, ly, static_node=False, get_gcc=True)
        g.set_random_edges_fail(0.24)  # threshold at 0.3473 (bond percolation)
        g.get_largest_connected_component_size()
        self.assertTrue(g.largest_component_size / g.n_nodes < 0.2)
        g.reset_loss()
        g.set_random_edges_fail(0.44)  # threshold at 0.3473 (bond percolation)
        g.get_largest_connected_component_size()
        self.assertTrue(g.largest_component_size / g.n_nodes > 0.8)
        g.reset_loss()
        g.set_random_nodes_lost(0.38)  # threshold at 0.5 (site percolation)
        g.get_largest_connected_component_size()
        self.assertTrue(g.largest_component_size / g.n_nodes < 0.2)
        g.reset_loss()
        g.set_random_nodes_lost(0.62)  # threshold at 0.5 (site percolation)
        g.get_largest_connected_component_size()
        self.assertTrue(g.largest_component_size / g.n_nodes > 0.45)

        self.assertTrue(time.time() - time_start < 1.0)  # check running time

    def test_quantum_fusion(self):
        time_start = time.time()

        lx, ly, lz = 50, 50, 50
        g = GraphPlusConstruction(np.product([lx, ly, lz]))
        g.create_3d_simple_cubic(lx, ly, lz, interrupt=(1, 1, 1), shift=(0, 0, 0), static_node=False, get_gcc=False)
        g.target_graph_by_bell_fusion(0.5, 0.95)
        self.assertTrue(g.check_start_target_percolation() == 0)
        g.reset_loss()
        g.target_graph_by_bell_fusion(0.5, 0.97)
        self.assertTrue(g.check_start_target_percolation() == 1)

        g = GraphPlusConstruction(np.product([lx, ly, lz]))
        g.create_3d_simple_cubic(lx, ly, lz, interrupt=(1, 1, 1), shift=(0, 0, 0), static_node=False, get_gcc=True)  # 1655117247
        g.target_graph_by_bell_fusion(0.5, 0.95)
        g.get_largest_connected_component_size()
        self.assertTrue(g.largest_component_size / g.n_nodes < 0.1)
        g.reset_loss()
        g.target_graph_by_bell_fusion(0.5, 0.97)
        g.get_largest_connected_component_size()
        self.assertTrue(g.largest_component_size / g.n_nodes > 0.2)

        g = GraphPlusConstruction(np.product([lx, ly, lz]))
        g.create_3d_simple_cubic(lx, ly, lz, interrupt=(1, 1, 1), shift=(0, 0, 0), static_node=True, get_gcc=False)
        g.target_graph_by_bell_fusion(0.5, 0.93)
        self.assertTrue(g.check_start_target_percolation() == 0)
        g.reset_loss()
        g.target_graph_by_bell_fusion(0.5, 0.96)
        self.assertTrue(g.check_start_target_percolation() == 1)

        g = GraphPlusConstruction(np.product([lx, ly, lz]))
        g.create_3d_simple_cubic(lx, ly, lz, interrupt=(1, 1, 1), shift=(0, 0, 0), static_node=True, get_gcc=True)  # 1655131436
        g.target_graph_by_bell_fusion(0.5, 0.93)
        g.get_largest_connected_component_size()
        self.assertTrue(g.largest_component_size / g.n_nodes < 0.1)
        g.reset_loss()
        g.target_graph_by_bell_fusion(0.5, 0.96)
        g.get_largest_connected_component_size()
        self.assertTrue(g.largest_component_size / g.n_nodes > 0.2)

        g = GraphPlusConstruction(np.product([lx, ly, lz]))
        g.create_3d_simple_cubic(lx, ly, lz, interrupt=(3, 1, 1), shift=(0, 0, 0), static_node=False, get_gcc=False)
        g.target_graph_by_bell_fusion(0.5, 0.92)
        self.assertTrue(g.check_start_target_percolation() == 0)
        g.reset_loss()
        g.target_graph_by_bell_fusion(0.5, 0.95)
        self.assertTrue(g.check_start_target_percolation() == 1)

        g = GraphPlusConstruction(np.product([lx, ly, lz]))
        g.create_3d_simple_cubic(lx, ly, lz, interrupt=(3, 1, 1), shift=(0, 0, 0), static_node=False, get_gcc=True)  # 1655278659
        g.target_graph_by_bell_fusion(0.5, 0.92)
        g.get_largest_connected_component_size()
        self.assertTrue(g.largest_component_size / g.n_nodes < 0.1)
        g.reset_loss()
        g.target_graph_by_bell_fusion(0.5, 0.95)
        g.get_largest_connected_component_size()
        self.assertTrue(g.largest_component_size / g.n_nodes > 0.15)

        g = GraphPlusConstruction(np.product([lx, ly, lz]))
        g.create_3d_raussendorf(lx, ly, lz, static_node=False, get_gcc=False)
        g.apply_loss(0.75)
        self.assertTrue(g.check_start_target_percolation() == 0)
        g.reset_loss()
        g.apply_loss(0.80)
        self.assertTrue(g.check_start_target_percolation() == 1)
        g.reset_loss()

        g = GraphPlusConstruction(np.product([lx, ly, lz]))
        g.create_3d_raussendorf(lx, ly, lz, static_node=False, get_gcc=True)  # 1655715987
        g.apply_loss(0.76)
        g.get_largest_connected_component_size()
        self.assertTrue(g.largest_component_size / g.n_nodes < 0.05)
        g.reset_loss()
        g.apply_loss(0.80)
        g.get_largest_connected_component_size()
        self.assertTrue(g.largest_component_size / g.n_nodes > 0.1)

        self.assertTrue(time.time() - time_start < 40)  # check running time


if __name__ == '__main__':
    unittest.main()
