import numpy as np
from graph_construction import GraphPlusConstruction
import matplotlib.pyplot as plt
import dill
import time
import shutil
import multiprocessing
from itertools import product
from policy import do_fusion_policy


# set up a certain graph, apply loss etc., return relative size of the largest connected component
def percolation_process(p, q, n_rep, graph, dummy):
    graph.reset_loss()  # reset previous losses (edges, nodes)
    graph.target_graph_by_bell_fusion(p, q, n_rep=n_rep)  # aim at graph with fusions
    # graph.target_graph_by_bell_fusion_adaptive(p, q, adaptive_level=2, hyperparam=[n_rep], fusion_policy=do_fusion_policy)  # more advanced policy
    # graph.apply_loss(q)  # apply loss to existing graph
    # graph.set_random_edges_fail(p)  # for classical bond-percolation
    # graph.set_random_nodes_lost(q)  # for classical site-percolation
    if len(graph.set_start) == 0:  # if there is no start set for percolation defined
        graph.get_largest_connected_component_size()
        return graph.largest_component_size / graph.n_nodes
    else:
        return graph.check_start_target_percolation()


if __name__ == '__main__':
    pool = multiprocessing.Pool()
    pool = multiprocessing.Pool(processes=4)
    time_stamp = time.time()
    load_session = False
    save_session = True
    path_name_load = '1655131436.pkl'
    folder_save = ''

    if load_session:
        dill.load_session(path_name_load)
    else:
        if save_session:
            shutil.copy2('percolation_parallel.py', folder_save + str(round(time_stamp)) + '_' + 'percolation_parallel.py')

        l_param = ['probability: edge stays/fusion works', 'probability: node/photon stays', 'simulation box size',
                   'fusion repetition', 'fusion frequency', 'shift', 'static node', 'percolation method']  # parameters that can be tuned

        # parameters that are no property of the graph
        l_p = [0.5]  # probability that edge stays/fusion works
        l_q = [0.8+0.005*i for i in range(40)] # probability that node/photon stays
        n_repeat = 1  # number of repetitions for averaging

        # parameters that determine the graph
        l_box = [[30, 30, 30]]  # size, dimension of the simulation
        l_n_rep = [0]  # no. of times a fusion is repeated if it fails (repeating makes only sense if l_static[]==True)
        l_inter = [(1, 1, 1)]  # if > 1, only edges at interruption point can be lost
        l_shift = [(0, 0, 0)]  # in addition to l_inter
        l_static = [False]  # define nodes that cannot be lost
        l_gcc = [True]  # decide between methods to determine percolation

        simulation_shape = (len(l_p), len(l_q), len(l_box), len(l_n_rep), len(l_inter), len(l_shift), len(l_static), len(l_gcc))
        simulation_result = np.zeros(simulation_shape)
        parallelize_shape = (len(l_p), len(l_q), n_repeat)

        for i_box in range(len(l_box)):
            for i_rep in range(len(l_n_rep)):
                for i_inter in range(len(l_inter)):
                    for i_shift in range(len(l_shift)):
                        for i_static in range(len(l_static)):
                            for i_gcc in range(len(l_gcc)):
                                # construct graph for given parameters
                                g = GraphPlusConstruction(np.product(l_box[i_box]))
                                #g.create_nd_simple_cubic(l_box[i_box], static_node=l_static[i_static], get_gcc=l_gcc[i_gcc])
                                g.create_3d_simple_cubic(l_box[i_box][0], l_box[i_box][1], l_box[i_box][2], interrupt=l_inter[i_inter], shift=l_shift[i_shift], static_node=l_static[i_static], get_gcc=l_gcc[i_gcc])
                                #g.create_2d_cluster_square(l_box[i_box][0], l_box[i_box][1], static_node=l_static[i_static], get_gcc=l_gcc[i_gcc])
                                #g.create_3d_raussendorf(l_box[i_box][0], l_box[i_box][1], l_box[i_box][2], static_node=l_static[i_static], get_gcc=l_gcc[i_gcc])
                                # perform simulation for different parameters in parallel
                                results = pool.starmap(percolation_process, product(l_p, l_q, [l_n_rep[i_rep]], [g], [i for i in range(n_repeat)]))
                                simulation_result[:,:,i_box,i_rep,i_inter,i_shift,i_static,i_gcc] = np.mean(np.array(results).reshape(parallelize_shape), axis=len(parallelize_shape) - 1)
        time_stop = time.time()
        print('time passed: ' + str(time_stop - time_stamp))
        if save_session:
            pool.close()
            del pool  # terminate, delete pool object as it cannot be saved with dill
            dill.dump_session(folder_save + str(round(time_stamp)) + '.pkl')

    # plot simulation results
    if len(l_p) == 1:
        sweep_axis = np.argmax(simulation_shape[2:]) + 2
        idx = [0 for _ in range(len(simulation_shape) - 2)]
        for i_sweep in range(simulation_shape[sweep_axis]):
            idx[sweep_axis - 2] = i_sweep
            plt.plot(l_q, simulation_result[0, :, idx[0], idx[1], idx[2], idx[3], idx[4]])
        plt.xlabel('probability that node/photon stays')
        plt.ylabel('largest connected component 0-1')
        plt.title('sweep: ' + l_param[sweep_axis])
        plt.savefig(folder_save + str(round(time_stamp)) + '.pdf')
        plt.show()
    elif len(l_q) == 1:
        sweep_axis = np.argmax(simulation_shape[2:]) + 2
        idx = [0 for _ in range(len(simulation_shape) - 2)]
        for i_sweep in range(simulation_shape[sweep_axis]):
            idx[sweep_axis - 2] = i_sweep
            plt.plot(l_p, simulation_result[:, 0, idx[0], idx[1], idx[2], idx[3], idx[4]])
        plt.xlabel('probability that edge stays/fusion works')
        plt.ylabel('largest connected component 0-1')
        plt.title('sweep: ' + l_param[sweep_axis])
        plt.savefig(folder_save + str(round(time_stamp)) + '.pdf')
        plt.show()
    else:
        idx = [0 for _ in range(len(simulation_shape) -1 - 2)]
        fig, ax = plt.subplots()
        c = ax.imshow(simulation_result[:, :, idx[0], idx[1], idx[2], idx[3], idx[4]], cmap='viridis', vmin=0, vmax=1, extent=[min(l_q), max(l_q), min(l_p), max(l_p)], interpolation='nearest', origin='lower')
        fig.colorbar(c, ax=ax)
        ax.set_title('largest connected component 0-1, param = ' + str(idx))
        plt.xlabel('q: probability that node/photon stays')
        plt.ylabel('p: probability that edge stays/fusion works')
        plt.savefig(folder_save + str(round(time_stamp)) + '.pdf')
        plt.show()

        # plot line cut
        p_cut = 0.5
        p_cut_idx = np.where(np.array(l_p) == p_cut)[0][0]
        plt.plot(l_q, simulation_result[p_cut_idx, :, idx[0], idx[1], idx[2], idx[3], idx[4]])
        plt.title('p_fusion = ' + str(p_cut))
        plt.xlabel('probability that node/photon stays')
        plt.ylabel('largest connected component 0-1')
        plt.savefig(folder_save + str(round(time_stamp)) + '_line_cut.pdf')
        plt.show()
