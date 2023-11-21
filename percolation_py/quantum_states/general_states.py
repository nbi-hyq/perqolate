import numpy as np


# general quantum state
class GeneralState:
    def __init__(self, n_qubit):
        self.n_qubit = n_qubit
        self.prefactor = np.ones(2**n_qubit)
        self.basis = np.zeros((2**n_qubit, n_qubit), dtype=bool)
        self.is_hypergraph_state = False
        self.is_weighted_hypergraph_state = False
        for i in range(2**n_qubit):
            b = np.binary_repr(i, width=n_qubit)
            for j in range(n_qubit):
                self.basis[i, j] = bool(int(b[j]))

    # set the prefactor of all basis states
    def set_prefactor(self, prefactor):
        if np.all(np.abs(prefactor) > 0):
            self.prefactor = prefactor

    # verify that the state is a hypergraph state (erases also global phase factor)
    def get_is_hypergraph_state(self, precision=10**-15):
        if self.prefactor[0] == 0:
            self.is_hypergraph_state = False
        else:
            self.prefactor = self.prefactor / self.prefactor[0]  # erase global phase (not that state 00..00 cannot get a - sign by Z, CZ)
            if np.all(np.absolute(np.absolute(self.prefactor) - 1) < precision):
                self.is_hypergraph_state = np.all(np.isreal(self.prefactor))  # all prefactors are +-1
                self.is_weighted_hypergraph_state = True
            else:
                self.is_hypergraph_state = False
                self.is_weighted_hypergraph_state = False

    # apply multi-qubit controlled Z operation (h is a tuple of all involved qubit indices)
    def apply_cz(self, h):
        flip_sign = np.ones(2 ** self.n_qubit, dtype=bool)
        for v in h:
            flip_sign = np.logical_and(flip_sign, self.basis[:, v])
        self.prefactor = self.prefactor * (-1)**flip_sign  # flip sign if all qubits CZ acts on are in 1

    # apply an arbitrary operation on a single qubit (can be unitary gate, or projector), TBD: vectorize
    def apply_1qubit_operation(self, idx, gate):
        if 0 <= idx < self.n_qubit:
            prefactor = np.zeros(2 ** self.n_qubit)
            offset = 2 ** (self.n_qubit - idx - 1)  # distance to next state where just the idx qubit differs (basis states are ordered)
            for i in range(2 ** self.n_qubit):
                if self.basis[i, idx] == 0:
                    r = np.dot(gate, np.array([[self.prefactor[i]],[self.prefactor[i + offset]]]))
                    prefactor[i] += r[0, 0]
                    prefactor[i + offset] += r[1, 0]
            self.prefactor = prefactor
        else:
            print("operation on undefined qubit")

    # perform a measurement in z-direction on one qubit
    def measure_z(self, idx, result):
        if result == 0:
            projector = np.array([[1,0], [0,0]])
        else:
            projector = np.array([[0,0], [0,1]])
        p0 = np.sum(np.multiply(self.prefactor, self.prefactor))  # p0 = 1 if state is normalized
        self.apply_1qubit_operation(idx, projector)
        p1 = np.sum(np.multiply(self.prefactor, self.prefactor))
        print("projection probability: ", str(p1/p0))

    # perform a Bell measurement
    # idx1, idx2: index of the two qubits that are measured
    # result1, result2: result of the Bell measurement (projector that is applied)
    def measure_bell(self, idx1, idx2, result1, result2):
        H = 1 / np.sqrt(2) * np.array([[1, 1], [1, -1]])
        self.apply_1qubit_operation(idx2, H)
        self.apply_cz((idx1, idx2))
        self.apply_1qubit_operation(idx1, H)
        self.apply_1qubit_operation(idx2, H)
        self.measure_z(idx1, result1)
        self.measure_z(idx2, result2)

    # print the state as sum of basis states
    def print_basis(self):
        for i in range(2 ** self.n_qubit):
            if self.prefactor[i] == +1:
                sign = '+'
            elif self.prefactor[i] == -1:
                sign = '-'
            else:
                sign = str(self.prefactor[i]) + " * "
            print(sign + np.binary_repr(i, width=self.n_qubit))
