from graph_construction import GraphPlusConstruction
import matplotlib.pyplot as plt
import numpy as np

np.random.seed(267239)  # fix the seed

# plot site percolation curve of classical graph (just for illustration)

lsize = 7
l_p = [0.01*i for i in range(100)]
l_size = []
cnt = 0

g = GraphPlusConstruction(lsize*lsize)
g.create_2d_cluster_square(lsize, lsize, static_node=False, get_gcc=True)

for p in l_p:
    g.reset_loss()  # reset previous losses (edges, nodes)
    g.set_random_nodes_lost(p)  # for classical bond-percolation
    gcc = g.get_largest_connected_component_size(get_gcc=True)
    l_size.append(g.largest_component_size / g.n_nodes)

    if divmod(cnt, 10)[1] == 0:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for x in range(lsize):
            for y in range(lsize):
                i_node = x * lsize + y
                # plot the node
                if g.node_lost[i_node]:
                    plt.plot([x], [y], 'o', color='#C86464')
                else:
                    plt.plot([x], [y], 'o', color='#000000')
                if i_node in gcc:
                    plt.plot([x], [y], 'o', color='#64C864', markeredgecolor='k')
                # plot the edges
                for nb in g.neighbors[i_node]:
                    x_nb, y_nb = divmod(nb, lsize)
                    if not (g.node_lost[i_node] or g.node_lost[nb] or (i_node, nb) in g.edges_fail or (nb, i_node) in g.edges_fail):
                        plt.plot([x, x_nb],[y, y_nb], '-k', zorder=-1)
                    elif i_node < nb:
                        plt.plot([x, x_nb], [y, y_nb], ':', color='#BEBEBE', zorder=-1)
        plt.title(str(p))
        ax.set_aspect('equal')
        ax.axis('off')
        fig.set_size_inches(2.5, 2.5)
        plt.savefig('percolation_site' + str(cnt) + '.pdf')
        plt.show()
    cnt += 1

fig = plt.figure()
plt.plot(l_p,l_size)

# repeat the simulation for a larger graph and without plots
lsize = 200
l_size = []

g = GraphPlusConstruction(lsize*lsize)
g.create_2d_cluster_square(lsize, lsize, static_node=False, get_gcc=True)

for p in l_p:
    g.reset_loss()  # reset previous losses (edges, nodes)
    g.set_random_nodes_lost(p)  # for classical bond-percolation
    g.get_largest_connected_component_size()
    l_size.append(g.largest_component_size / g.n_nodes)

plt.plot(l_p,l_size)
plt.xlabel('$p_{site}$')
plt.ylabel('component size')
fig.set_size_inches(2.3, 2.0)
plt.gca().set_ylim([0,1])
plt.savefig('percolation_site_result.pdf')
plt.show()
