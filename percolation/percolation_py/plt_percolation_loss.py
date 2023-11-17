from graph_construction import GraphPlusConstruction
import matplotlib.pyplot as plt
import numpy as np

np.random.seed(246789)  # fix the seed

plt_lcc = False  # highlight largest connected component
lsize_x = 10
lsize_y = 5
l_p = [0.5 + 0.01*i for i in range(51)]
l_size = []
cnt = 0

g = GraphPlusConstruction(lsize_x * lsize_y)
g.create_2d_cluster_square(lsize_x, lsize_y, static_node=False, get_gcc=True)

for p in l_p:
    g.reset_loss()  # reset previous losses (edges, nodes)
    g.apply_loss(p)  # for classical bond-percolation
    gcc = g.get_largest_connected_component_size(get_gcc=True)
    l_size.append(g.largest_component_size / g.n_nodes)

    if divmod(cnt, 10)[1] == 0:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for x in range(lsize_x):
            for y in range(lsize_y):
                i_node = x * lsize_y + y
                # plot the node
                if g.node_lost[i_node]:
                    plt.plot([x], [y], 'o', color='#C86464')
                elif g.node_z[i_node]:
                    plt.plot([x], [y], 'o', color='#BEBEBE')
                else:
                    plt.plot([x], [y], 'o', color='#000000')
                if plt_lcc and i_node in gcc:
                    plt.plot([x], [y], 'o', color='#64C864', markeredgecolor='k')
                # plot th edges
                for nb in g.neighbors[i_node]:
                    x_nb, y_nb = divmod(nb, lsize_y)
                    if not (g.node_lost[i_node] or g.node_lost[nb] or g.node_z[i_node] or g.node_z[nb] or (i_node, nb) in g.edges_fail or (nb, i_node) in g.edges_fail):
                        plt.plot([x, x_nb],[y, y_nb], '-k', zorder=-1)
                    elif i_node < nb:
                        plt.plot([x, x_nb], [y, y_nb], ':', color='#BEBEBE', zorder=-1)
        plt.title(str(p))
        ax.set_aspect('equal')
        ax.axis('off')
        fig.set_size_inches(5, 2.5)
        plt.savefig('percolation_loss' + str(cnt) + '.pdf')
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
    g.apply_loss(p)  # for classical bond-percolation
    g.get_largest_connected_component_size()
    l_size.append(g.largest_component_size / g.n_nodes)

plt.plot(l_p,l_size)
plt.xlabel('1 - $p_{loss}$')
plt.ylabel('component size')
fig.set_size_inches(2.3, 2.0)
plt.gca().set_ylim([0,1])
plt.savefig('percolation_loss_result.pdf')
plt.show()
