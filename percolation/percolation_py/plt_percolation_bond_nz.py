from graph_construction import GraphPlusConstruction
import matplotlib.pyplot as plt
import random

random.seed(235679)  # fix the seed

# plot illustration of Newman-Ziff algorithm for bond-perolation 
# (just for plotting illustrations, no efficient version of the Newman-Ziff algorithm)
lsize = 7
l_size = []
cnt = 0

g = GraphPlusConstruction(lsize*lsize)
g.create_2d_cluster_square(lsize, lsize, static_node=True, get_gcc=True)
g.get_edges() # create set with all edges
g.edges_fail = g.edges  # set all lost in the beginning
l_edges = list(g.edges)
random.shuffle(l_edges)  # list with randomly shuffled edges to loop through

for e in l_edges:
    g.edges_fail.remove(e) # add a new random edge
    gcc = g.get_largest_connected_component_size(get_gcc=True)
    l_size.append(g.largest_component_size / g.n_nodes)

    if len(l_edges)/1.80 > cnt > len(l_edges)/2.20:  # divmod(cnt, 10)[1] == 0:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for x in range(lsize):
            for y in range(lsize):
                i_node = x * lsize + y
                # plot edges
                for nb in g.neighbors[i_node]:
                    x_nb, y_nb = divmod(nb, lsize)
                    if not ((i_node, nb) in g.edges_fail or (nb, i_node) in g.edges_fail):
                        plt.plot([x, x_nb],[y, y_nb], '-k', zorder=-1)
                    elif i_node < nb:
                        plt.plot([x, x_nb], [y, y_nb], ':', color='#BEBEBE', zorder=-1)
                # plot node
                plt.plot([x], [y], 'o', color='#000000')
                if i_node in gcc:
                    plt.plot([x],[y],'o', color='#64C864', markeredgecolor='k')
        plt.title(str(cnt+1) + ' / ' + str(len(l_edges)))
        ax.set_aspect('equal')
        ax.axis('off')
        fig.set_size_inches(2.5, 2.5)
        plt.savefig('percolation_bond_nz' + str(cnt) + '.pdf', format='pdf')
        plt.show()
    cnt += 1

fig = plt.figure()
plt.plot([i for i in range(len(l_size))],l_size)
plt.xlabel('#bonds')
plt.ylabel('component size')
fig.set_size_inches(2.3, 2.0)
plt.gca().set_ylim([0,1])
plt.savefig('percolation_bond_nz_result.pdf', format='pdf')
plt.show()
