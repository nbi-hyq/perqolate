import sys
sys.path.insert(0, '../../quantum_states')  # add directory above for import
from graph_states import GraphState
from scipy.spatial import ConvexHull
from scipy.optimize import linprog
import matplotlib.pyplot as plt
from copy import copy
import numpy as np
from itertools import permutations


# convert decimal number to digital vector in system defined by box
def index_to_digital_vec(i, box):
    len_vec = len(box)
    factor = 1
    vec = [0 for _ in range(len_vec)]
    for pos in range(len_vec - 1, -1, -1):
        vec[pos] = (i % (factor * box[pos]) / factor)
        i -= vec[pos] * factor
        factor *= box[pos]
    return vec


# convert vector representing number in some digital system back to decimal index
def digital_vec_to_index(vec, box):
    i = 0
    factor = 1
    for pos in range(len(box) - 1, -1, -1):
        i += vec[pos] * factor
        factor *= box[pos]
    return int(i)


# Add methods for construction in this child class
class GraphPlusConstruction(GraphState):
    def __init__(self, n_nodes):
        GraphState.__init__(self, n_nodes)

    # example graph (rectangular 2d lattice)
    def create_2d_cluster_square(self, lx, ly, interrupt=(1, 1), shift=(0, 0), static_node=False, get_gcc=True):
        # loop over rectangular lattice, x has higher value in digit system for counting the nodes
        # interrupt[0] == 2: fusion every 2nd edge in x-direction (default: star cluster -- fusion for every edge)
        # shift[0] (shift[1]): shift line by line where the interruption in x(y)-direction is
        # static: make all nodes static (not removed in classical bond-site percolation, no photon loss in quantum case)
        if static_node:
            self.static_nodes = np.ones(self.n_nodes, dtype='bool')
        for x in range(lx):
            for y in range(ly):
                if x < lx-1:  # add the edge towards the right
                    self.add_edge(x * ly + y, (x + 1) * ly + y, fusion=((x + y*shift[0]) % interrupt[0] == 0))
                if y < ly-1:  # add the edge towards the bottom
                    self.add_edge(x * ly + y, x * ly + y + 1, fusion=((y + x*shift[1]) % interrupt[1] == 0))
        # define start and target nodes to decide if percolation happened
        if not get_gcc:
            for y in range(ly):
                self.set_start.add(y)
                self.set_target.add((lx - 1) * ly + y)

    # example graph (triangular 2d lattice) == rectangular with one diagonal
    def create_2d_cluster_triangular(self, lx, ly, static_node=False, get_gcc=True):
        if static_node:
            self.static_nodes = np.ones(self.n_nodes, dtype='bool')
        for x in range(lx):
            for y in range(ly):
                if x < lx-1:  # add the edge towards the right
                    self.add_edge(x * ly + y, (x + 1) * ly + y)
                if y < ly-1:  # add the edge towards the bottom
                    self.add_edge(x * ly + y, x * ly + y + 1)
                if x < lx-1 and y < ly-1:  # add diagonal edge
                    self.add_edge(x * ly + y, (x + 1) * ly + y + 1)
        # define start and target nodes to decide if percolation happened
        if not get_gcc:
            for y in range(ly):
                self.set_start.add(y)
                self.set_target.add((lx - 1) * ly + y)

    # example graph (honeycomb 2d lattice)
    def create_2d_cluster_hexagonal(self, lx, ly, static_node=False, get_gcc=True):
        if static_node:
            self.static_nodes = np.ones(self.n_nodes, dtype='bool')
        for x in range(lx):
            for y in range(ly):
                if x < lx-1:  # add the edge towards the right
                    self.add_edge(x * ly + y, (x + 1) * ly + y)
                if y < ly-1 and (y + x) % 2 == 0:  # add the edge towards the bottom
                    self.add_edge(x * ly + y, x * ly + y + 1)
        # define start and target nodes to decide if percolation happened
        if not get_gcc:
            for y in range(ly):
                self.set_start.add(y)
                self.set_target.add((lx - 1) * ly + y)

    # example graph (simple cubic lattice)
    def create_3d_simple_cubic(self, lx, ly, lz, interrupt=(1, 1, 1), shift=(0, 0, 0), static_node=False, get_gcc=True):
        if static_node:
            self.static_nodes = np.ones(self.n_nodes, dtype='bool')
        # loop over rectangular lattice, x, y have higher value in digit system for counting the nodes
        for x in range(lx):
            for y in range(ly):
                for z in range(lz):
                    if x < lx-1:  # add the edge in x-direction
                        self.add_edge(x * ly * lz + y * lz + z, (x+1) * ly * lz + y * lz + z, fusion=((x + (y+z)*shift[0]) % interrupt[0] == 0))
                    if y < ly-1:  # add the edge in y-direction
                        self.add_edge(x * ly * lz + y * lz + z, x * ly * lz + (y+1) * lz + z, fusion=((y + (x+z)*shift[1]) % interrupt[1] == 0))
                    if z < lz-1:  # add the edge in z-direction
                        self.add_edge(x * ly * lz + y * lz + z, x * ly * lz + y * lz + (z+1), fusion=((z + (x+y)*shift[2]) % interrupt[2] == 0))
        # define start and target nodes to decide if percolation happened
        if not get_gcc:
            for y in range(ly):
                for z in range(lz):
                    self.set_start.add(y * lz + z)
                    self.set_target.add((lx - 1) * ly * lz + y * lz + z)

    # example graph (Raussendorf lattice): create from 3d simple cubic by removing edges
    # !note: that means that the relative size of the largest connected component is underestimated
    def create_3d_raussendorf(self, lx, ly, lz, static_node=False, get_gcc=True):
        self.create_3d_simple_cubic(lx, ly, lz, static_node=static_node)
        for x in range(0, lx, 2):
            for y in range(0, ly, 2):
                for z in range(0, lz, 2):
                    self.rm_node(x * ly * lz + y * lz + z)
                    if x < lx-1 and y < ly-1 and z < lz-1:
                        self.rm_node((x + 1) * ly * lz + (y + 1) * lz + z + 1)
        # define start and target nodes to decide if percolation happened (note some start/target nodes are always detached)
        if not get_gcc:
            for y in range(ly):
                for z in range(lz):
                    self.set_start.add(y * lz + z)
                    self.set_target.add((lx - 1) * ly * lz + y * lz + z)

    # example graph (hypercubic lattice, Z^d), bcc: option adding diagonals {-1, +1}**dimension
    def create_nd_simple_cubic(self, box, static_node=False, bcc=False, get_gcc=True):
        if static_node:
            self.static_nodes = np.ones(self.n_nodes, dtype='bool')
        for i in range(self.n_nodes):
            i_vec = index_to_digital_vec(i, box)
            for digit_pos in range(len(box)):
                if i_vec[digit_pos] < box[digit_pos] - 1:
                    i_vec[digit_pos] += 1
                    self.add_edge(i, digital_vec_to_index(i_vec, box))
                    i_vec[digit_pos] -= 1
            if bcc:
                num_bcc = 2 ** (len(box) - 1)  # -1 avoids adding edges twice, len(box) is dimension
                for i_d in range(num_bcc):
                    i_d_vec = np.array(index_to_digital_vec(i_d, [2 for _ in range(len(box))])) * 2 - 1
                    i_vec += i_d_vec
                    if np.all(i_vec >= 0) and np.all(i_vec < np.array(box)):
                        self.add_edge(i, digital_vec_to_index(i_vec, box))  # add diagonal connection
                    i_vec -= i_d_vec
            if not get_gcc:
                if i_vec[0] == 0:
                    self.set_start.add(i)
                elif i_vec[0] == box[0] - 1:
                    self.set_target.add(i)

    # get lattice from vectors representing the bonds between nodes on hypercubic lattice
    def get_lattice_from_vectors(self, box, vecs, static_node=False, get_gcc=True, periodic=False):
        if static_node:
            self.static_nodes = np.ones(self.n_nodes, dtype='bool')
        num_vecs = np.shape(vecs)[0]  # number of vectors specifying the lattice
        dimension = np.shape(vecs)[1]  # dimension of the lattice
        for i in range(self.n_nodes):
            i_vec = np.array(index_to_digital_vec(i, box))
            i_vec_save = copy(i_vec)
            for k in range(num_vecs):
                i_vec += vecs[k, :]  # step in opposite direction will be done when the neighbor is called
                if np.all(i_vec >= 0) and np.all([i_vec[d]<box[d] for d in range(dimension)]):
                    i_nb = digital_vec_to_index(i_vec, box)
                    self.add_edge(i, i_nb)
                elif periodic:
                    for d in range(dimension):
                        i_vec[d] = i_vec[d] % box[d]
                    i_nb = digital_vec_to_index(i_vec, box)
                    self.add_edge(i, i_nb)
                i_vec = i_vec_save
            if not get_gcc:
                if i_vec[0] == 0:
                    self.set_start.add(i)
                elif i_vec[0] == box[0] - 1:
                    self.set_target.add(i)

    # 3-photon fusion square lattice from https://doi.org/10.1038/s41467-019-08948-x
    # (note: the construction here makes only sense to verify the bond-percolation threshold without losses)
    def create_cluster_square_from_2d_3photon(self, lx, ly, get_gcc=True):
        # loop over rectangular lattice, x, y have higher value in digit system for counting the nodes
        lz = 3
        for x in range(lx):
            for y in range(ly):
                self.add_edge(x * ly * lz + y * lz, x * ly * lz + y * lz + 1)  # 1st edge in z-direction
                self.add_edge(x * ly * lz + y * lz + 1, x * ly * lz + y * lz + 2)  # 2nd edge in z-direction
                if x < lx-1:
                    self.add_edge(x * ly * lz + y * lz, (x+1) * ly * lz + y * lz)  # z == 0 level (only in x-direction)
                if y < ly-1:
                    self.add_edge(x * ly * lz + y * lz + 2, x * ly * lz + (y+1) * lz + 2)  # z == 2 level (only in y direction)
        # define start and target nodes to decide if percolation happened
        if not get_gcc:
            for y in range(ly):
                for z in range(lz):
                    self.set_start.add(y * lz + z)
                    self.set_target.add((lx - 1) * ly * lz + y * lz + z)

