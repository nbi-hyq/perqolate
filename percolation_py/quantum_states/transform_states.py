import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from copy import copy


# compute hyper-edges from signs of basis states (hy: hyper-graph state, gnrl: general state)
def compute_edges(hy, gnrl, keep_prefactor=True):
    if not hy.n_nodes == gnrl.n_qubit:
        print("hy and gnrl have different number of qubits")
        return 1
    else:
        gnrl.get_is_hypergraph_state()
        if not gnrl.is_hypergraph_state:
            print("edges cannot be computed as gnrl is not a hyper-graph state")
            return 1
    hy.edges = set()
    if keep_prefactor:
        prefactor_save = copy(gnrl.prefactor)
    for i in range(2 ** hy.n_nodes):
        if gnrl.prefactor[i] == -1:
            h = tuple(np.nonzero(gnrl.basis[i, :])[0])
            hy.edges.add(h)
            gnrl.apply_cz(h)
    if keep_prefactor:
        gnrl.prefactor = prefactor_save
    return 0


# compute signs of basis states from hyper-edges
# gnrl: general state, hy: hyper-graph state, reset: reset gnrl.prefactor to empty graph before apply CZ gates
def compute_basis(gnrl, hy, reset=True):
    if not hy.n_nodes == gnrl.n_qubit:
        print("hy and gnrl have different number of qubits")
        return 1
    if reset:
        gnrl.set_prefactor(np.ones(2 ** gnrl.n_qubit))
    for h in hy.edges:
        gnrl.apply_cz(h)
    return 0

# compute stabilizer generators and generator matrix (1|adj-matrix)
# st: stabilizer states, hy: hyper-graph state
# TBD: other direction, tests
def compute_stabilizers_from_hypergraph(st, hy):
    if not hy.n_nodes == st.n_qubit:
        print("hy and st have different number of qubits")
        return 1
    else:
        smallest_edge, largest_edge = hy.get_edge_boundaries()
        if not (largest_edge == 2 and smallest_edge == 2):
            print("hy is not a stabilizer state")
            return 1
    st.stabilizers = [[i] for i in range(st.n_qubit)]  # a stabilizer is a list starting with i for X_i
    for n in range(hy.n_nodes):
        for h in hy.edges:
            if n in h and len(h) > 1:
                neighbors = set(h)
                neighbors.remove(n)
                st.stabilizers[n].append(tuple(neighbors))  # every tuple represents a CZ operation
    for e in hy.edges:
        st.generator_matrix[e[0], hy.n_nodes + e[1]] = 1
        st.generator_matrix[e[1], hy.n_nodes + e[0]] = 1
    st.stabilizer_array = np.zeros((st.n_qubit, st.n_qubit), dtype='uint8')
    for i in range(hy.n_nodes):
        st.stabilizer_array[i, i] = 1  # X operator represented by 1
        for j in hy.neighbors[i]:
            st.stabilizer_array[i, j] = 3  # Z operator on neighbors represented by 3
    return 0
