import numpy as np
import copy

# TBD: implement algorithms from https://iopscience.iop.org/article/10.1088/1367-2630/7/1/170/meta
# TBD: keep track of stabilizer prefactors, compare two stabilizers - identical, difference, get state/code space
class StabilizerState:
    # stabilizer_array: matrix of Pauli-matrices with every row representing a stabilizer
    # it has elements 0-3 with 0: identity, 1: Pauli X, 2: Pauli Y, 3: Pauli Z
    def __init__(self, stabilizer_array):
        self.n_qubit = stabilizer_array.shape[1]
        self.stabilizer_array = np.array(stabilizer_array, dtype='uint8')  # make sure the data-type is correct
        self.generator_matrix = np.hstack((np.eye(self.n_qubit), np.zeros((self.n_qubit, self.n_qubit))))  # different representation of stabilizers
        self.stabilizers = []  # represent stabilizers as list
        self.height_function = np.zeros(self.n_qubit + 1)
        self.rank = None
        a = np.multiply(np.multiply(np.multiply(stabilizer_array, stabilizer_array - 1), stabilizer_array - 2), stabilizer_array - 3)
        if np.any(a):
            self.valid = False
            print('No valid stabilizer. Some elements are not in the Pauli group.')
        else:
            z = np.zeros((self.n_qubit, self.n_qubit))
            i = np.eye(self.n_qubit)
            p = np.hstack((np.vstack((z, i)), np.vstack((i, z))))

            x_matrix = (self.stabilizer_array == 1)
            y_matrix = (self.stabilizer_array == 2)
            z_matrix = (self.stabilizer_array == 3)
            self.xz_matrix = np.hstack((x_matrix + y_matrix, z_matrix + y_matrix))  # note: + acts like logical or

            if np.mod(np.dot(np.dot(self.xz_matrix, p), np.transpose(self.xz_matrix)), 2).any():
                self.valid = False
                print('No valid stabilizer. Some operators do not commute.')
            else:
                self.valid = True

    # transpose (swap) two rows in the stabilizer array
    def transpose_rows(self, row_a, row_b):
        if row_a != row_b:
            row_a_save = copy.copy(self.stabilizer_array[row_a,:])
            self.stabilizer_array[row_a,:] = copy.copy(self.stabilizer_array[row_b,:])
            self.stabilizer_array[row_b,:] = row_a_save

    # multiply row a on row b (following rules for Pauli matrices)
    def multiply_a_on_b(self, row_a, row_b):
        for n in range(self.n_qubit):
            if self.stabilizer_array[row_a, n] > 0:
                if self.stabilizer_array[row_a, n] == self.stabilizer_array[row_b, n]:
                    self.stabilizer_array[row_b, n] = 0
                elif self.stabilizer_array[row_a, n] == 1:
                    if self.stabilizer_array[row_b, n] == 2:
                        self.stabilizer_array[row_b, n] = 3
                    elif self.stabilizer_array[row_b, n] == 3:
                        self.stabilizer_array[row_b, n] = 2
                    else:
                        self.stabilizer_array[row_b, n] = 1
                elif self.stabilizer_array[row_a, n] == 2:
                    if self.stabilizer_array[row_b, n] == 1:
                        self.stabilizer_array[row_b, n] = 3
                    elif self.stabilizer_array[row_b, n] == 3:
                        self.stabilizer_array[row_b, n] = 1
                    else:
                        self.stabilizer_array[row_b, n] = 2
                elif self.stabilizer_array[row_a, n] == 3:
                    if self.stabilizer_array[row_b, n] == 2:
                        self.stabilizer_array[row_b, n] = 1
                    elif self.stabilizer_array[row_b, n] == 1:
                        self.stabilizer_array[row_b, n] = 2
                    else:
                        self.stabilizer_array[row_b, n] = 3

    # get row-reduced echelon form (RREF) representing the stabilizer generators + rank of the stabilizers
    # (algorithm RREF from https://iopscience.iop.org/article/10.1088/1367-2630/7/1/170/meta)
    def get_rref(self):
        K = self.stabilizer_array.shape[0]
        N = self.stabilizer_array.shape[1]
        k_u = 0  # defines active region
        n_l = 0  # defines active region
        
        while k_u < K and n_l < N:  # step 3): check that active region is non-empty
            # 1) get number of different Pauli operators in column
            pauli_in_column = np.array([0,0,0])
            k_first_pauli = -1  # first row with non identity operator
            k_second_pauli = -1  # second row with different non identity operator
            for k in range(k_u, K):
                if self.stabilizer_array[k, n_l] > 0:
                    pauli_in_column[self.stabilizer_array[k, n_l] - 1] = 1
                    if k_first_pauli < 0:
                        k_first_pauli = k
                    elif (k_second_pauli < 0) and (self.stabilizer_array[k_first_pauli, n_l] != self.stabilizer_array[k, n_l]):
                        k_second_pauli = k
            
            # 2) row multiplication and swapping (transposing) rows
            if np.sum(pauli_in_column) == 0:  # case a)
                n_l += 1
            elif np.sum(pauli_in_column) == 1:  # case b)
                self.transpose_rows(k_first_pauli, k_u)
                for k in range(k_u+1, K):
                    if self.stabilizer_array[k_u, n_l] == self.stabilizer_array[k, n_l]:
                        self.multiply_a_on_b(k_u, k)
                k_u += 1
                n_l += 1
            else:  # case c)
                self.transpose_rows(k_first_pauli, k_u)
                self.transpose_rows(k_second_pauli, k_u + 1)
                for k in range(k_u+2, K):
                    if self.stabilizer_array[k, n_l] > 0:
                        if self.stabilizer_array[k_u, n_l] == self.stabilizer_array[k, n_l]:
                            self.multiply_a_on_b(k_u, k)
                        elif self.stabilizer_array[k_u+1, n_l] == self.stabilizer_array[k, n_l]:
                            self.multiply_a_on_b(k_u+1, k)
                        else:
                            self.multiply_a_on_b(k_u, k)
                            self.multiply_a_on_b(k_u+1, k)
                k_u += 2
                n_l += 1
        
        # compute the rank (count all rows of stabilizer array that are not only identities)
        self.rank = 0
        for k in range(K):
            if not np.all(np.zeros(N) == self.stabilizer_array[k,:]):
                self.rank += 1

    # get height function as in https://www.nature.com/articles/s41534-022-00522-6
    # (so far only for self.n_qubit == self.rank meaning that the stabilizers define a state, not a state space)
    def get_height_function(self):
        first_nontrivial = np.zeros(self.n_qubit)
        for y in range(self.stabilizer_array.shape[0]):
            for x in range(self.n_qubit):
                if self.stabilizer_array[y, x] > 0:
                    first_nontrivial[y] = x
                    break
        self.height_function = np.zeros(self.n_qubit + 1)
        for x in range(self.n_qubit + 1):
            number_left = np.sum(first_nontrivial >= x)
            self.height_function[x] = self.n_qubit - x - number_left

    # print the stabilizers
    def print_stabilizers(self):
        for k in range(self.stabilizer_array.shape[0]):
            stabilizer_str = ""
            for n in range(self.stabilizer_array.shape[1]):
                if self.stabilizer_array[k, n] == 1:
                    stabilizer_str += ("X" + str(n))
                elif self.stabilizer_array[k, n] == 2:
                    stabilizer_str += ("Y" + str(n))
                elif self.stabilizer_array[k, n] == 3:
                    stabilizer_str += ("Z" + str(n))
            print(stabilizer_str)
