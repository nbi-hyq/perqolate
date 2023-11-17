from copy import copy
import numpy as np
import matplotlib.pyplot as plt
import unittest
import sys
sys.path.insert(0, '../quantum_states')  # add directory above for import
from graph_states import GraphState
from general_states import GeneralState
from stabilizer_states import StabilizerState
from transform_states import compute_edges, compute_basis, compute_stabilizers_from_hypergraph

plot = False
np.random.seed(527829)  # fix the seed

class TestQuantumStates(unittest.TestCase):
    def test_get_generator_matrix(self):
        g = GraphState(3)
        s = StabilizerState(np.zeros((4, 3)))
        g.set_edges(set([(0, 1), (1, 2)]))
        compute_stabilizers_from_hypergraph(s, g)
        h = (s.generator_matrix == np.array([[1, 0, 0, 0, 1, 0], [0, 1, 0, 1, 0, 1], [0, 0, 1, 0, 1, 0]]))
        self.assertTrue(h.all())
        g.set_edges(set([(0, 1), (1, 2), (1,), (1, 2, 3)]))
        smallest_edge, largest_edge = g.get_edge_boundaries()
        self.assertEqual(largest_edge, 3)
        self.assertEqual(smallest_edge, 1)

    def test_compute_basis(self):
        g = GraphState(3)
        q = GeneralState(3)
        g.set_edges(set([(0, 1), (1, 2)]))
        compute_basis(q, g)
        self.assertTrue((q.prefactor == np.array([1, 1, 1, -1, 1, 1, -1, 1])).all)

    def test_compute_edges(self):
        g = GraphState(3)
        q = GeneralState(3)
        q.prefactor = np.array([1, 1, 1, 1, 1, 1, 1, -1])
        compute_edges(g, q, keep_prefactor=False)
        self.assertEqual(g.edges, {(0, 1, 2)})
        g = GraphState(4)
        q = GeneralState(4)
        q.prefactor = np.array([1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, -1, 1, -1, -1])
        compute_edges(g, q, keep_prefactor=False)
        self.assertEqual(g.edges, {(0, 1), (0, 1, 2, 3), (3,)})

    def test_is_hypergraph_state(self):
        g = GeneralState(3)
        g.prefactor = (-0.7*np.exp(-1j)) * np.array([1, 1, 1, 1, 1, 1, 1, -1])
        g.get_is_hypergraph_state()
        self.assertTrue(g.is_hypergraph_state)
        g = GeneralState(3)
        g.prefactor = np.array([1, 1, 1, 1, 1, 1, 0.9, -1])  # no hyper-graph state
        g.get_is_hypergraph_state()
        self.assertFalse(g.is_hypergraph_state)
        g = GeneralState(3)
        g.prefactor = np.array([0, 1, 1, 1, 1, 1, 1, -1])  # no hyper-graph state
        g.get_is_hypergraph_state()
        self.assertFalse(g.is_hypergraph_state)
        g = GeneralState(3)
        g.prefactor = np.array([np.exp(-1j), 1, 1, np.exp(-15j), 1, 1, 1, -1])  # no hyper-graph state (just weighted hyper-graph state)
        g.get_is_hypergraph_state()
        self.assertFalse(g.is_hypergraph_state)
        self.assertTrue(g.is_weighted_hypergraph_state)

    def test_back_and_forth(self):
        n_qubit = 5
        g = GraphState(n_qubit)
        q = GeneralState(n_qubit)
        q.prefactor = np.rint(np.random.rand(2**n_qubit)) * 2 - 1  # random signs (+-1)
        prefactor_cpy = copy(q.prefactor)  # save signs in REW representation
        compute_edges(g, q, keep_prefactor=True)  # compute hyper-edges from REW representation (keep REW representation)
        self.assertTrue((prefactor_cpy == q.prefactor).all)  # check that REW stayed identical (keep_prefactor=True)
        g2 = GraphState(n_qubit)
        compute_edges(g2, q, keep_prefactor=False)  # compute hyper-edges again from REW representation
        self.assertTrue(np.all(q.prefactor == 1))  # after compute_edges with keep_pre-factor=False: empty graph in REW
        compute_basis(q, g)  # recompute the signs in REW representation
        self.assertTrue((prefactor_cpy == q.prefactor).all)  # check that REW is identical to original
        compute_basis(q, g, reset=False)  # apply all CZ gates a second time (no reset)
        self.assertTrue(np.all(q.prefactor == 1))  # all CZ gates 2x must result in empty graph (all REW signs +)
        if plot:
            fig, axs = plt.subplots(1, 2)
            fig.suptitle('test_back_and_forth')
            g.plt_graph(axs[0])
            g2.plt_graph(axs[1])
            plt.show()

    def test_adj(self):
        n_qubit = 4
        adj = np.array([[1, 1, 1, 0], [1, 1, 0, 1], [1, 0, 1, 1], [0, 1, 1, 1]], dtype='uint32')
        g = GraphState(n_qubit)
        g.set_adj(adj)
        g.get_edges_from_adj()
        self.assertEqual(g.edges, {(0, 1), (1, 3), (0, 2), (2, 3)})
        g.set_adj(np.eye(g.n_nodes, dtype='uint32'))
        g.get_adj_from_edges()
        self.assertTrue(np.sum(np.sum(g.adj - adj)) == 0)

    # test effect of projective measurement on graph state, TBD: check measurement effect on hypergraph state
    def test_measure_z(self):
        n_qubit = 5
        g = GraphState(n_qubit)
        q = GeneralState(n_qubit)
        g.set_edges(set([(0, 1), (1, 2), (2, 0), (2, 3), (3, 4)]))
        if plot:
            fig, axs = plt.subplots(1, 2)
            fig.suptitle('test_measure_z: measure_z(2, 0)')
            g.plt_graph(axs[0])
        compute_basis(q, g, reset=True)
        q.measure_z(2, 0)  # measure 0 on qubit 2
        q.apply_1qubit_operation(2, 1/np.sqrt(2)*np.array([[1,1],[1,-1]]))  # bring measured qubit to |+> state (such that everything is a graph state)
        compute_edges(g, q, keep_prefactor=False)  # compute hyper-edges from REW representation
        self.assertEqual(g.edges, {(0, 1), (3, 4)})
        if plot:
            g.plt_graph(axs[1])
            plt.show()
        g.set_edges(set([(0, 1), (1, 2), (2, 0), (2, 3), (3, 4)]))
        compute_basis(q, g, reset=True)
        q.measure_z(2, 1)  # measure 1 on qubit 2
        q.apply_1qubit_operation(2, 1 / np.sqrt(2) * np.array([[-1, 1], [1, 1]]))  # bring measured qubit to |+> state
        compute_edges(g, q, keep_prefactor=False)
        self.assertEqual(g.edges, {(0,), (0, 1), (1,), (3,), (3, 4)})
        g.set_edges(set([(1, 2), (2, 0, 1), (2, 3), (3, 4, 0)]))  # measurement effect on hyper-graph?
        if plot:
            fig, axs = plt.subplots(1, 2)
            fig.suptitle('test_measure_z: measure_z(2, 1)')
            g.plt_graph(axs[0])
        compute_basis(q, g, reset=True)
        q.measure_z(2, 1)  # measure 1 on qubit 2
        q.apply_1qubit_operation(2, 1 / np.sqrt(2) * np.array([[-1, 1], [1, 1]]))  # bring measured qubit to |+> state
        compute_edges(g, q, keep_prefactor=False)
        if plot:
            g.plt_graph(axs[1])
            plt.show()

    # test effect on Bell measurement on graph state (fusion), TBD: effect on hyper-graph state
    def test_measure_bell(self):
        n_qubit = 4
        g = GraphState(n_qubit)
        q = GeneralState(n_qubit)
        g.set_edges(set([(0, 1), (2, 3)]))
        if plot:
            fig, axs = plt.subplots(1, 2)
            fig.suptitle('test_measure_bell: measure_bell(1, 2, 0, 0)')
            g.plt_graph(axs[0])
        compute_basis(q, g, reset=True)
        q.apply_1qubit_operation(2, 1 / np.sqrt(2) * np.array([[1, 1], [1, -1]]))  # apply Hadamard gate
        q.measure_bell(1, 2, 0, 0)  # perform Bell measurement on qubits (fusion)
        q.apply_1qubit_operation(1, 1 / np.sqrt(2) * np.array([[1, 1], [1, -1]]))  # reset measured qubit to |+>
        q.apply_1qubit_operation(2, 1 / np.sqrt(2) * np.array([[1, 1], [1, -1]]))  # reset measured qubit to |+>
        self.assertTrue(np.all(np.round(2 * q.prefactor) == np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, -1, 1, -1, 1, -1, 1, -1])))
        compute_edges(g, q, keep_prefactor=False)  # rounding issue?
        if plot:
            g.plt_graph(axs[1])
            plt.show()

        # Bell measurement of B01 rather than B00
        n_qubit = 6
        g = GraphState(n_qubit)
        q = GeneralState(n_qubit)
        g.set_edges(set([(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3)]))
        if plot:
            fig, axs = plt.subplots(1, 2)
            fig.suptitle('test_measure_bell: measure_bell(2, 3, 0, 1)')
            g.plt_graph(axs[0])
        compute_basis(q, g, reset=True)
        q.apply_1qubit_operation(2, 1 / np.sqrt(2) * np.array([[1, 1], [1, -1]]))  # apply Hadamard gate
        q.measure_bell(2, 3, 0, 1)  # perform Bell measurement on qubits (fusion)
        q.apply_1qubit_operation(2, 1 / np.sqrt(2) * np.array([[1, 1], [1, -1]]))  # reset measured qubit to |+>
        q.apply_1qubit_operation(3, 1 / np.sqrt(2) * np.array([[1, 1], [1, -1]]))  # reset measured qubit to |+>
        compute_edges(g, q, keep_prefactor=False)
        if plot:
            g.plt_graph(axs[1])
            plt.show()

        # make edge disappear if it is supposed to be added as edge but is an edge already (like in graph complementation)
        n_qubit = 6
        g = GraphState(n_qubit)
        g.set_edges(set([(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 5)]))
        if plot:
            fig, axs = plt.subplots(1, 2)
            fig.suptitle('test_measure_bell: measure_bell(2, 3, 0, 0)')
            g.plt_graph(axs[0])
        compute_basis(q, g, reset=True)
        q.apply_1qubit_operation(2, 1 / np.sqrt(2) * np.array([[1, 1], [1, -1]]))  # apply Hadamard gate
        q.measure_bell(2, 3, 0, 0)  # perform Bell measurement on qubits (fusion)
        q.apply_1qubit_operation(2, 1 / np.sqrt(2) * np.array([[1, 1], [1, -1]]))  # reset measured qubit to |+>
        q.apply_1qubit_operation(3, 1 / np.sqrt(2) * np.array([[1, 1], [1, -1]]))  # reset measured qubit to |+>
        compute_edges(g, q, keep_prefactor=False)
        if plot:
            g.plt_graph(axs[1])
            plt.show()

        # TBD: simulate failed Bell measurement!


if __name__ == '__main__':
    unittest.main()
