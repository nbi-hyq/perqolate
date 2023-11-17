import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import copy
from functools import reduce
from operator import mul, and_, add
from itertools import chain, combinations


# Undirected graph / graph states / can be used as hyper-graph (state), can be used for percolation simulations
class GraphState:
    def __init__(self, n_nodes=None, graph=None):
        self.n_nodes = n_nodes
        self.graph = None
        self.stab_gens = None
        self.n_edges = 0
        self.neighbors = [set() for _ in range(self.n_nodes)]  # neighbors of each node
        self.edges = set()  # set with edges (do not add edges for large graphs)
        self.adj = np.array([], dtype='uint32')  # adjacency matrix
        self.no_fusion_edges = set()  # graph state only: edges that are already present in an initial micro-cluster
        self.static_nodes = np.zeros(self.n_nodes, dtype='bool')  # graph state only: nodes that are not subject to photon loss (e.g. QD)
        self.edges_fail = set()  # edges that are temporarily not usable (due to e.g. fusion fail)
        self.node_lost = np.zeros(self.n_nodes, dtype='bool')  # node/photon that is lost
        self.node_z = np.zeros(self.n_nodes, dtype='bool')  # node where Z-measurement must be done
        self.largest_component_size = 0  # size of the largest connected component
        self.visited = np.zeros(self.n_nodes, dtype='bool')  # visited nodes (needed for crawling of graph)
        self.node_queue = [0]  # order of node traversal for breadth first crawling (avoids too deep recursion in DFS)
        self.node_queue_pos = 0  # position in the queue for breath first traversal
        self.set_start = set()  # where graph traversal starts (if check_start_target_percolation is used)
        self.set_target = set()  # where graph traversal terminates (if check_start_target_percolation is used)
        self.num_success = np.zeros(self.n_nodes, dtype=int)  # number of successful fusions for each node
        self.num_photon = np.zeros(self.n_nodes, dtype=int)  # number of photons left
        for n in range(self.n_nodes):
            self.num_photon[n] = len(self.neighbors[n])

    # add an undirected new edge to the graph
    def add_edge(self, i, j, fusion=True):
        if not i in self.neighbors[j]:  # if not edge yet
            self.n_edges += 1
            self.neighbors[i].add(j)
            self.neighbors[j].add(i)
            if not fusion:
                if i < j:  # convention: smaller index first (undirected graph)
                    self.no_fusion_edges.add((i, j))
                elif i > j:
                    self.no_fusion_edges.add((j, i))

    # set hyper-edges (TBD: check if valid edges)
    def set_edges(self, edges):
        self.edges = edges

    # remove node by disconnecting it (can be used to construct graph from other graph)
    def rm_node(self, n):
        self.n_edges -= len(self.neighbors[n])
        for i_nb in self.neighbors[n]:
            self.neighbors[i_nb].remove(n)
        self.neighbors[n] = set()

    # remove an undirected edge from the graph (can be used to construct graph from other graph)
    def rm_edge(self, i, j):
        if i in self.neighbors[j]:
            self.n_edges -= 1
            self.neighbors[i].remove(j)
            self.neighbors[j].remove(i)

    # reset everything that is related to crawling the graph
    def reset_crawling(self):
        self.largest_component_size = 0
        self.visited = np.zeros(self.n_nodes)
        self.node_queue = [0]
        self.node_queue_pos = 0

    # reset losses applied to edges, nodes
    def reset_loss(self):
        self.edges_fail = set()
        self.node_lost = np.zeros(self.n_nodes, dtype='bool')
        self.node_z = np.zeros(self.n_nodes, dtype='bool')

    # graph states only: given a graph state, photons are lost + Z is applied to neighbors
    def apply_loss(self, q):  # q: probability that a photon is not lost
        loss = (np.random.random(self.n_nodes) > q)
        for i_node in range(self.n_nodes):
            if loss[i_node] and self.static_nodes[i_node] == False:  # static node stays (e.g. quantum emitter)
                self.node_lost[i_node] = 1
                for i_neighbor in self.neighbors[i_node]:
                    if not ((i_node, i_neighbor) in self.edges_fail or (i_neighbor, i_node) in self.edges_fail):
                        self.node_z[i_neighbor] = 1 # mark neighbors with Z

    # graph state only: target at current graph by fusing micro-clusters
    # -> output is the obtained graph state under loss and finite fusion success
    def target_graph_by_bell_fusion(self, p, q, n_rep=0):
        # p: probability that fusion works (if both fusion qubits exist)
        # q: probability that photon is not lost (efficiency)
        # n_rep: number of repetitions when fusion fails (locally adaptive: create extra fusion arms via emitter)
        # note: n_rep > 0 makes only sense when the central node is static_node (a quantum emitter)
        rnd = np.random.random(self.n_edges)
        p_fail = ((1-p) * q**2)**(n_rep + 1)  # probability that all attempted fusions fail (but no loss)
        p_success = q**2 * p * (1 - p_fail) / (1 - (1-p) * q**2)  # probability of successful fusion
        fail = (rnd < p_fail)  # fail of fusion (but no loss after n_rep)
        loss = np.logical_and(p_fail <= rnd, rnd < 1 - p_success)  # fusion attempts terminated by loss of fusion qubit
        loss_edge = set()
        i = 0
        for n in range(self.n_nodes):  # this, next for loop effectively loops over all undirected edges
            for nb in self.neighbors[n]:
                if nb > n and (not (n, nb) in self.no_fusion_edges):  # 1) avoids looking at the same edge twice, 2) not consider at no_fusion edges
                    if loss[i]:  # note: rm edge is not necessary here (will automatically happen at the end)
                        self.node_z[n] = 1  # photon loss in fusion: mark both qubits with Z
                        self.node_z[nb] = 1 # photon loss in fusion : mark both qubits with Z
                        loss_edge.add((n, nb))
                    elif fail[i]:  # fusion failed (otherwise, do nothing to the graph if fusion succeeds)
                        self.edges_fail.add((n, nb))
                    i += 1
        self.apply_loss(q)  # apply loss to all nodes (qubits) on target graph, takes care of deleting at the end
        return loss_edge

    # measure some node in Z, update self.sum_success, and remove successful edges from self.edges
    def mz_and_update_num_success_after_mz(self, n):
        self.node_z[n] = True
        for nb in self.neighbors[n]:
            e = (min(n, nb), max(n, nb))
            if e in self.edges:
                self.edges.remove(e)
                self.edges_fail.add(e)
                self.num_success[n] -= 1
                self.num_success[nb] -= 1

    # make a fusion or just measure both qubits in Z depending on some policy
    def adaptive_fusion_action(self, p, q, n, nb, fusion_policy, hyperparam, adaptive_level):
        n_edges_alive_before = len(self.edges)  # number of edges before fusion action
        edge = (min(n, nb), max(n, nb))
        self.num_photon[n] -= 1  # fusion photon will be consumed
        self.num_photon[nb] -= 1  # fusion photon will be consumed
        if self.node_z[n] and self.node_z[nb]:
            self.edges_fail.add(edge)  # only effect: label the edge to not consider it 2nd time
        elif self.node_z[nb] and adaptive_level > 0:  # not make fusion attempt, measure fusion qubit in Z
            if np.random.random(1) > q:  # if fusion qubit lost when measure it
                self.mz_and_update_num_success_after_mz(n)
            self.edges_fail.add(edge)
        elif self.node_z[n] and adaptive_level > 0:  # not make fusion attempt, measure fusion qubit in Z
            if np.random.random(1) > q:  # if fusion qubit lost when measure it
                self.mz_and_update_num_success_after_mz(nb)
            self.edges_fail.add(edge)
        else:
            data_in = np.array([self.num_success[n], self.num_photon[n], self.num_success[nb], self.num_photon[nb], len(self.edges), self.n_nodes], dtype=int)
            if adaptive_level < 2 or fusion_policy(data_in, hyperparam[0]):  # make a fusion attempt
                if np.random.random(1) > q ** 2:
                    self.mz_and_update_num_success_after_mz(n)  # photon loss in fusion: mark both qubits with Z
                    self.mz_and_update_num_success_after_mz(nb)  # photon loss in fusion: mark both qubits with Z
                    self.edges_fail.add(edge)
                elif np.random.random(1) > p:  # fusion failed (otherwise, do nothing to the graph if fusion succeeds)
                    self.edges_fail.add(edge)
                else:
                    self.num_success[n] += 1  # count successful fusion
                    self.num_success[nb] += 1  # count successful fusion
                    self.edges.add(edge)
            else:  # do not make fusion attempt, measure fusion qubits in Z
                if np.random.random(1) > q:  # if fusion qubit lost when measuring it
                    self.mz_and_update_num_success_after_mz(n)
                if np.random.random(1) > q:  # if fusion qubit of neighbor lost when measuring it
                    self.mz_and_update_num_success_after_mz(nb)
                self.edges_fail.add(edge)
        reward = len(self.edges) - n_edges_alive_before  # no. edges after fusion action - no. edges before
        return reward

    # graph state only: target at current graph by fusing micro-clusters (adaptive)
    def target_graph_by_bell_fusion_adaptive(self, p, q, adaptive_level=2, hyperparam=[1.0], fusion_policy=lambda x, y:True):
        # p: probability that fusion works (if both fusion qubits exist)
        # q: probability that photon is not lost (efficiency)
        # adaptive_level: 0 not adaptive, 1: not connect to node_z[i]==True, 2 more clever policy
        self.reset_loss()
        self.edges = set()  # reset, used as list for all edges that are alive
        self.num_success = np.zeros(self.n_nodes, dtype=int)  # reset: number of successful fusions for each node
        for n in range(self.n_nodes):  # reset: number of photons left
            self.num_photon[n] = len(self.neighbors[n])
        order = np.arange(self.n_nodes)
        np.random.shuffle(order)  # go through nodes in random order
        # loop over all edges by looping over nodes and their neighbors in given order
        for n in order:
            for nb in self.neighbors[n]:
                edge = (min(n, nb), max(n, nb))
                if not edge in self.no_fusion_edges and not edge in self.edges_fail and not edge in self.edges:  # not consider no_fusion edges, nor edges considered already
                    rwd = self.adaptive_fusion_action(p, q, n, nb, fusion_policy, hyperparam, adaptive_level)
        self.apply_loss(q)  # apply loss to all nodes (qubits) on target graph, takes care of deleting at the end

    # set edges as non-usable (for classical bond-percolation)
    def set_random_edges_fail(self, p):  # p: probability that an edge stays
        rand = (np.random.random(self.n_edges) > p)
        i = 0
        for n in range(self.n_nodes):  # this, next for loop effectively loops over all undirected egdes
            for nb in self.neighbors[n]:
                if nb > n:  # avoids looking at the same edge twice
                    if rand[i]:
                        self.edges_fail.add((n, nb))
                    i += 1

    # temporarily remove node by setting them to "lost" (for classical site-percolation)
    def set_random_nodes_lost(self, q):  # q: probability that a node stays
        rand = (np.random.random(self.n_nodes) > q)
        self.node_lost = np.logical_and(rand, np.logical_not(self.static_nodes))

    # crawl the graph following the queue for breadth-first traversal
    def crawl_graph(self, no_node):
        while self.node_queue_pos < len(self.node_queue):
            idx = self.node_queue[self.node_queue_pos]
            for idx_neighbor in self.neighbors[idx]:
                if (not self.visited[idx_neighbor]) and (no_node[idx_neighbor] == False) and (not (idx, idx_neighbor) in self.edges_fail) and (not (idx_neighbor, idx) in self.edges_fail):
                    self.node_queue.append(idx_neighbor)  # append neighbors to queue if not visited
                    self.visited[idx_neighbor] = True  # mark node as visited
            self.node_queue_pos += 1

    # get the size of the largest connected component of the graph
    def get_largest_connected_component_size(self, get_gcc=False):
        self.reset_crawling()
        gcc = [0]
        no_node = np.logical_or(self.node_z, self.node_lost)
        for i_node in range(self.n_nodes):
            if not self.visited[i_node] and no_node[i_node] == False:  # not visited in a previous crawling, nor loss/Z
                self.visited[i_node] = True  # mark node as visited
                self.node_queue = [i_node]  # initialize the queue for breadth first traversal
                self.node_queue_pos = 0  # start at queue beginning
                self.crawl_graph(no_node)  # crawl graph with breath first traversal
                if get_gcc and self.node_queue_pos > self.largest_component_size:
                    gcc = copy.copy(self.node_queue)
                self.largest_component_size = max(self.node_queue_pos, self.largest_component_size)
        return gcc

    # crawl the graph following the queue for breadth-first traversal (stop when target is reached)
    def crawl_graph_to_target(self, no_node):
        while self.node_queue_pos < len(self.node_queue):
            idx = self.node_queue[self.node_queue_pos]
            for idx_neighbor in self.neighbors[idx]:
                if (not self.visited[idx_neighbor]) and (no_node[idx_neighbor] == False) and (not (idx, idx_neighbor) in self.edges_fail) and (not (idx_neighbor, idx) in self.edges_fail):
                    if idx_neighbor in self.set_target:
                        return 1  # target is reached with at least one path
                    self.node_queue.append(idx_neighbor)  # append neighbors to queue if not visited
                    self.visited[idx_neighbor] = True  # mark node as visited
            self.node_queue_pos += 1
        return 0

    # check if there is a connection from a set of start-nodes to some target-nodes
    def check_start_target_percolation(self):
        self.reset_crawling()
        no_node = np.logical_or(self.node_z, self.node_lost)
        for i_node in self.set_start:
            if not self.visited[i_node] and no_node[i_node] == False:  # not visited in a previous crawling, nor loss/Z
                self.visited[i_node] = True  # mark node as visited
                self.node_queue = [i_node]  # initialize the queue for breadth first traversal
                self.node_queue_pos = 0  # start at queue beginning
                if self.crawl_graph_to_target(no_node):
                    return 1  # target is reached with at least one path
        return 0

    # get edges from node neighbors representation
    def get_edges(self):
        self.edges = set()
        for i in range(self.n_nodes):
            for nb in self.neighbors[i]:
                if nb > i:
                    self.edges.add((i, nb))

    # set adjacency matrix
    def set_adj(self, adj):
        if np.shape(adj)[0] == np.shape(adj)[1] == self.n_nodes:
            self.adj = np.array(adj, dtype='uint32')

    # get edges and node neighbor representation from adjacency matrix
    def get_edges_from_adj(self):
        self.neighbors = [set() for _ in range(self.n_nodes)]  # reset node neighbor representation
        for i in range(np.shape(self.adj)[0]):
            for j in range(np.shape(self.adj)[0]):
                if self.adj[i, j] == 1 and not i==j:
                    self.neighbors[i].add(j)
        self.get_edges()  # get set with all edges

    # apply multi-qubit controlled Z operation (corresponds to adding corresponding hyper-edge)
    # h: hyper-edge e.g. (1,) or (1,4,6)
    def apply_cz(self, h):
        hy = tuple(sorted(h))  # make sure that there are not several equivalent hyper-edges
        if hy in self.edges:
            self.edges.remove(hy)  # remove hyper-edge
        else:
            self.edges.add(hy)

    # get the size of the largest and the smallest hyper-edge
    def get_edge_boundaries(self):
        largest_edge = 1
        smallest_edge = self.n_nodes
        for e in self.edges:
            largest_edge = max(largest_edge, len(e))
            smallest_edge = min(smallest_edge, len(e))
        return smallest_edge, largest_edge

    # True: if normal graph, False: if there are hyper-edges
    def is_normal_graph(self):
        _, largest_edge = self.get_edge_boundaries()
        if largest_edge <= 2:
            return True
        else:
            return False

    # get adjacency matrix from edges
    def get_adj_from_edges(self):
        edge_min, edge_max = self.get_edge_boundaries()
        if edge_min == edge_max == 2:  # check that there are only normal edges
            self.adj = np.eye(self.n_nodes, dtype='uint32')
            for e in self.edges:
                self.adj[e[0], e[1]] = 1
                self.adj[e[1], e[0]] = 1

    # plot the hyper-graph state as hyper-graph (qubits with just Z, no other connections are not plotted)
    def plt_graph(self, ax=None):
        graph = nx.Graph()
        num_hyper_edges = 0
        set_z = set()  # store single qubit Z-operations here
        node_labels = {}
        for edge in self.edges:
            if len(edge) == 1:
                set_z.add(edge[0])
            elif len(edge) == 2:
                graph.add_edge(*edge)
                node_labels[edge[0]] = str(edge[0])
                node_labels[edge[1]] = str(edge[1])
            else:
                num_hyper_edges -= 1
                for v in edge:
                    graph.add_edge(num_hyper_edges, v)
                    node_labels[num_hyper_edges] = ''  # don't label the nodes representing hyper-edges
                    node_labels[v] = str(v)
        for v in set_z:
            if v in graph.nodes.keys():
                node_labels[v] += 'Z'

        node_list = list(graph.nodes)  # can be unordered: unclear how to assign color
        node_color = ['gray' if node_list[i] >= 0 else 'red' for i in range(len(graph.nodes))]
        layout = nx.spring_layout(graph)
        nx.draw_networkx_nodes(graph, layout, node_color=node_color, ax=ax)
        nx.draw_networkx_edges(graph, layout, ax=ax)
        nx.draw_networkx_labels(graph, layout, labels=node_labels, ax=ax)
