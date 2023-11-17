from graph_construction import GraphPlusConstruction
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse

np.random.seed(91045)  # fix the seed

# plot percolation simulation for constructing graph state from small pieces (in presence of photon loss)

plt_lcc = False  # highlight largest connected component
lsize_x = 10
lsize_y = 5
l_p = [0.8 + 0.005*i for i in range(41)]
l_size = []
cnt = 0

g = GraphPlusConstruction(lsize_x * lsize_y)
g.create_2d_cluster_square(lsize_x, lsize_y, static_node=False, get_gcc=True)

for p in l_p:
    g.reset_loss()  # reset previous losses (edges, nodes)
    lossy_edge = g.target_graph_by_bell_fusion(0.5, p)
    gcc = g.get_largest_connected_component_size(get_gcc=True)
    l_size.append(g.largest_component_size / g.n_nodes)

    if divmod(cnt, 5)[1] == 0:
        fig = plt.figure("before")
        ax = fig.add_subplot(111)
        fig2 = plt.figure("after")
        ax2 = fig2.add_subplot(111)
        for x in range(lsize_x):
            for y in range(lsize_y):
                i_node = x * lsize_y + y
                # plot the node
                if g.node_lost[i_node]:
                    plt.figure("before")
                    plt.plot([x], [y], 'o', color='#C86464')  # lost node
                    plt.figure("after")
                    plt.plot([x], [y], 'o', color='#C86464')  # lost node
                else:
                    plt.figure("before")
                    plt.plot([x], [y], 'o', color='#000000')  # all other nodes
                    plt.figure("after")
                    if g.node_z[i_node]:
                        plt.plot([x], [y], 'o', color='#BEBEBE')  # measure in Z node
                    else:
                        plt.plot([x], [y], 'o', color='#000000')  # all other nodes
                if plt_lcc and i_node in gcc:
                    plt.figure("after")
                    plt.plot([x],[y],'o', color='#64C864', markeredgecolor='k') # nodes of the largest connected component in green
                plt.figure("before")
                if x > 0:
                    plt.plot([x, x - 0.3], [y, y] ,'-k', zorder=-1)
                    plt.plot([x - 0.3], [y], '.k', zorder=-1)
                if y > 0:
                    plt.plot([x, x], [y, y - 0.3], '-k', zorder=-1)
                    plt.plot([x], [y - 0.3], '.k', zorder=-1)
                if x < lsize_x - 1:
                    plt.plot([x, x + 0.3], [y, y] ,'-k', zorder=-1)
                    plt.plot([x + 0.3], [y], '.k', zorder=-1)
                if y < lsize_y - 1:
                    plt.plot([x, x], [y, y + 0.3], '-k', zorder=-1)
                    plt.plot([x], [y + 0.3], '.k', zorder=-1)
                plt.figure("after")
                # draw edges
                for nb in g.neighbors[i_node]:
                    if i_node < nb:
                        x_nb, y_nb = divmod(nb, lsize_y)
                        if (i_node, nb) in lossy_edge or (nb, i_node) in lossy_edge: # draw all edges with a lossy qubit on it
                            plt.plot([x, x_nb], [y, y_nb], ':', color='#C86464', zorder=-1)
                            el = Ellipse(((x+x_nb)/2, (y+y_nb)/2), 0.3*(1+abs(x-x_nb)), 0.3*(1+abs(y-y_nb)), angle=0, color='#C86464', zorder=-1, alpha=100/255)
                            ax.add_artist(el)
                        elif (i_node, nb) in g.edges_fail or (nb, i_node) in g.edges_fail: # draw all other edges with failed fusion
                            plt.plot([x, x_nb], [y, y_nb], ':', color='#6464C8', zorder=-1)
                            el = Ellipse(((x+x_nb)/2, (y+y_nb)/2), 0.3*(1+abs(x-x_nb)), 0.3*(1+abs(y-y_nb)), angle=0, color='#6464C8', zorder=-1, alpha=100/255)
                            ax.add_artist(el)
                        elif g.node_lost[i_node] or g.node_lost[nb] or g.node_z[i_node] or g.node_z[nb]: # draw successful fusion between nodes that are measured in Z or lost as gray
                            plt.plot([x, x_nb], [y, y_nb], ':', color='#BEBEBE', zorder=-1)
                            el = Ellipse(((x+x_nb)/2, (y+y_nb)/2), 0.3*(1+abs(x-x_nb)), 0.3*(1+abs(y-y_nb)), angle=0, color='#64C864', zorder=-1, alpha=100/255)
                            ax.add_artist(el)
                        else: # all successful edges between existing nodes
                            plt.plot([x, x_nb], [y, y_nb], '-k', zorder=-1)
                            el = Ellipse(((x+x_nb)/2, (y+y_nb)/2), 0.3*(1+abs(x-x_nb)), 0.3*(1+abs(y-y_nb)), angle=0, color='#64C864', zorder=-1, alpha=100/255)
                            ax.add_artist(el)
        plt.figure("before")
        plt.title(str(p))
        ax.set_aspect('equal')
        ax.axis('off')
        fig.set_size_inches(5, 2.5)
        plt.savefig('percolation_fusion_with_loss_before' + str(cnt) + '.pdf')
        plt.figure("after")
        plt.title(str(p))
        ax2.set_aspect('equal')
        ax2.axis('off')
        fig2.set_size_inches(5, 2.5)
        plt.savefig('percolation_fusion_with_loss_after' + str(cnt) + '.pdf')
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
    g.target_graph_by_bell_fusion(0.5,p)  # for classical bond-percolation
    g.get_largest_connected_component_size()
    l_size.append(g.largest_component_size / g.n_nodes)

plt.plot(l_p,l_size)
plt.xlabel('1 - $p_{loss}$')
plt.ylabel('component size')
fig.set_size_inches(2.3, 2.0)
plt.gca().set_ylim([0,1])
plt.savefig('percolation_fusion_with_loss_result.pdf')
plt.show()
