from copy import copy
import numpy as np
import unittest
import sys
sys.path.insert(0, '../quantum_states')  # add directory above for import
from stabilizer_states import StabilizerState
from graph_states import GraphState
from transform_states import compute_stabilizers_from_hypergraph

plot_results = False

class TestStabilizer(unittest.TestCase):
    def test_stabilizer_init(self):
        s = StabilizerState(np.array([[4, 1], [0, 2]]))  # some element not in the Pauli group sigma_{0, 1, 2, 3}
        self.assertEqual(s.valid, False)
        s = StabilizerState(np.array([[1, 1, 0], [0, 2, 2], [0, 0, 0]]))  # 1st and 2nd element do not commute
        self.assertEqual(s.valid, False)
        s = StabilizerState(np.array([[1, 3, 3], [3, 1, 3], [3, 3, 1]]))  # stabilizer defining a triangle graph state
        self.assertEqual(s.valid, True)
        s = StabilizerState(np.array([[1, 3, 3], [3, 1, 3], [0, 2, 2]]))  # same as before but row 3 multiplied by row 2
        self.assertEqual(s.valid, True)

    def test_get_rref(self):
        s = StabilizerState(np.array([[1, 3, 3], [3, 1, 3], [2, 2, 0]])) # last row is linear dependent
        s.get_rref()
        self.assertEqual(s.rank, 2)
        s = StabilizerState(np.array([[1, 3, 3, 2], [2, 2, 0, 3], [3, 1, 3, 1], [2, 2, 0, 1]]))  # first three not independent
        s.get_rref()
        self.assertEqual(s.rank, 3)
        print(s.stabilizer_array)
        s.print_stabilizers()

    def test_get_height_function(self):
        import matplotlib.pyplot as plt
        s = StabilizerState(np.array([[1, 3, 3, 0], [3, 1, 3, 3], [3, 3, 1, 3], [0, 3, 3, 1]]))
        s.get_rref()
        s.get_height_function()
        self.assertTrue(abs(np.sum(s.height_function - np.array([0, 1, 2, 1, 0])))< 10**-10)

        # Fig. 2b in https://www.nature.com/articles/s41534-022-00522-6
        s = StabilizerState(np.array([[1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [3, 1, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3],
                                      [0, 0, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0], [0, 3, 3, 1, 0, 3, 0, 3, 0, 3, 0, 3],
                                      [0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 0], [0, 3, 0, 3, 3, 1, 0, 3, 0, 3, 0, 3],
                                      [0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0], [0, 3, 0, 3, 0, 3, 3, 1, 0, 3, 0, 3],
                                      [0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0], [0, 3, 0, 3, 0, 3, 0, 3, 3, 1, 0, 3],
                                      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3], [0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 3, 1]]))
        s.get_rref()
        s.get_height_function()
        self.assertTrue(abs(np.sum(s.height_function - np.array([0, 1, 1, 2, 1, 2, 1, 2, 1, 2, 1, 1, 0]))) < 10 ** -10)
        if plot_results:
            plt.plot(s.height_function, '-ob')
            plt.show()

        # Fig. 2c in https://www.nature.com/articles/s41534-022-00522-6
        s = StabilizerState(np.array([[1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [3, 1, 0, 0, 0, 3, 0, 3, 0, 3, 0, 3],
                                      [0, 0, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 3, 1, 0, 0, 0, 3, 0, 3, 0, 3],
                                      [0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 0], [0, 3, 0, 0, 3, 1, 0, 3, 0, 3, 0, 3],
                                      [0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0], [0, 3, 0, 3, 0, 3, 3, 1, 0, 0, 0, 3],
                                      [0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0], [0, 3, 0, 3, 0, 3, 0, 0, 3, 1, 0, 0],
                                      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3], [0, 3, 0, 3, 0, 3, 0, 3, 0, 0, 3, 1]]))
        s.get_rref()
        s.get_height_function()
        self.assertTrue(abs(np.sum(s.height_function - np.array([0, 1, 1, 2, 2, 2, 1, 2, 2, 2, 1, 1, 0]))) < 10 ** -10)
        if plot_results:
            plt.plot(s.height_function, '-ob')
            plt.show()

        # Fig. 2a in https://www.nature.com/articles/s41534-022-00522-6
        n_qubit = 12
        adj = np.array([[0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
                        [1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1], [0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1],
                        [0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1], [0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1],
                        [0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1], [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0]])
        s = StabilizerState(np.zeros((n_qubit, n_qubit)))
        g = GraphState(n_qubit)
        g.set_adj(adj)
        g.get_edges_from_adj()
        compute_stabilizers_from_hypergraph(s, g)
        s.get_rref()
        s.get_height_function()
        self.assertTrue(abs(np.sum(s.height_function - np.array([0, 1, 2, 3, 4, 5, 6, 5, 4, 3, 2, 1, 0]))) < 10 ** -10)
        if plot_results:
            g.plt_graph()
            plt.show()
            plt.plot(s.height_function, '-ob')
            plt.show()

if __name__ == '__main__':
    unittest.main()
