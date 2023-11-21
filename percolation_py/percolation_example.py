import numpy as np
from graph_construction import GraphPlusConstruction
import matplotlib.pyplot as plt
import dill
import time
import shutil
from itertools import product
from policy import do_fusion_policy


time_stamp = time.time()
load_session = False
save_session = True
path_name_load = '1670314534.pkl'
folder_save = ''

if load_session:
    dill.load_session(path_name_load)
else:
    if save_session:
        shutil.copy2('percolation_example.py', folder_save + str(round(time_stamp)) + '_' + 'percolation_example.py')

    # parameters that are not a property of the graph
    p = 0.5  # probability that edge stays/fusion works
    l_q = [0.8+0.005*i for i in range(40)]  # probability that node/photon stays

    # parameters that determine the graph
    box = [20, 20, 20]  # size, dimension of the simulation
    static = True  # central nodes cannot be lost
    gcc = True  # look for percolation (gcc==True looks for largest component size)

    # construct graph for given parameters
    g = GraphPlusConstruction(np.product(box))
    g.create_nd_simple_cubic(box, static_node=static, bcc=True, get_gcc=gcc)

    # perform percolation simulation
    result = []
    for q in l_q:
        g.reset_loss()  # reset previous losses (edges, nodes)
        # g.target_graph_by_bell_fusion(p, q, n_rep=2)  # aim at graph with fusions
        g.target_graph_by_bell_fusion_adaptive(p, q, adaptive_level=2, hyperparam=[6], fusion_policy=do_fusion_policy)  # not ballistic, but adaptive
        # g.apply_loss(q)  # apply loss to existing graph
        # g.set_random_edges_fail(p)  # for classical bond-percolation
        # g.set_random_nodes_lost(q)  # for classical site-percolation
        if gcc:
            g.get_largest_connected_component_size()
            result.append(g.largest_component_size / g.n_nodes)
        else:
            result.append(g.check_start_target_percolation())

    time_stop = time.time()
    print('time passed: ' + str(time_stop - time_stamp))
    if save_session:
        dill.dump_session(folder_save + str(round(time_stamp)) + '.pkl')

plt.plot(l_q, result)
plt.savefig(folder_save + str(round(time_stamp)) + '_' + '.pdf')
plt.show()
