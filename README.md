# perqolate
- ```C```-code for percolation simulations of graphs, graph states, fusion networks. The main part of the code is in ```percolation_c/```.
- Some simulation results can be found here: https://arxiv.org/abs/2304.03796

## percolation simulations: features
- Percolation simulations on undirected graphs/lattices. Percolation models: site-, bond-, and site-bond, bond-site
    - fast speed due to the use of the union-find algorithm by Newman and Ziff (https://doi.org/10.1103/PhysRevE.64.016706)
- Determine cluster spanning (percolation) or largest connected component size of graphs
- construct cluster states via rotated type-II fusion of star-shaped resource states. Check the resulting graph for percolation or the largest connected component size as a function of photon loss. Loss of photons used for the fusions is treated by removing/measuring the neighbors of both fusion photons. (see https://arxiv.org/abs/2304.03796 (Fig. A5(d)), https://doi.org/10.1103/PhysRevLett.115.020502 (Fig. 4) for the used linear optics fusion.)
- Simulation of fusion network using repeat-until-success method to boost the fusion probability. Photon loss can be simulated. (speed-up with modified Newman-Ziff algorithm)
- Simutate loss on graph states where loss is simulated by removing neighbors of a lost qubit from the graph by measuring them in the Z-basis. (speed-up with modified Newman-Ziff algorithm)
- Tests and continuous integration
- Examples and easy compilation via meson

## percolation simulations: overview

In the folder ```percolation_c``` you find a ```C``` implementation that forms the main part of the code. The folder ```percolation_py``` contains a ```Python```-only implementation that serves mostly as a reference for consistency checks and for plotting illustrations (and thus should not be used for extensive simulations). The following explanation follows the ```C``` implementation.

Generally, there are two types of simulations that can be done. The first type of simulation determines the size of the largest connected component of a graph after some random process such as removing edges or loss has been applied. The second type of simulation determines whether or not two sides of a graph lattice are connected (cluster spanning). In the given examples, the parameter ```get_size``` determines which simulation is used (```get_size==true``` for method1 and ```get_size==false``` for method2).

Given some initial graph, there are several types of random processes that can be applied to the graph and decreasing/increasing its connectivity. These processes are applied by the following functions that can be found in the file ```src/percolation_main.c```:
- ```percolate_site```: site-percolation simulation on a classical graph using the Newman-Ziff method to increase the speed.
- ```percolate_bond```: bond-percolation simulation on a classical graph using the Newman-Ziff method to increase the speed.
- ```apply_loss_nz```: apply loss to a graph representing a graph state. Similar to ```apply_loss```, ```apply_loss_bfs``` but more efficient.
- ```fusion_and_loss_nz```: the input graph represents a fusion network of star-shaped resource states (edges represent fusions).
  - considers qubit/photon loss, and uses a modified version of the Newman-Ziff algorithm to speed up the simulation.
- ```fusion_repeat_and_loss_nz```: similar to ```fusion_and_loss_nz``` but here fusions can be repeated until success. Assumes that two quantum emitters can generate two new fusion photons when the fusion fails. (uses modified Newman-Ziff method to speed up the simulation)
- ```get_expectation_value```: The input is a simulation of the largest connected component or the cluster spanning as a function of the number of lost photons, existing bond, or existing sites (result of a simulation with the Newman-Ziff method). Output: expectation value of the largest connected component or the cluster spanning as a function of the loss probability, bond occupation rate, or site occupation rate.

The following four functions implement similar things with different methods. They should serve only as references for consistency checks:
- ```apply_loss```: apply loss to a graph representing a quantum graph state, then determine the largest connected component or cluster spanning. (no Newman-Ziff method)
- ```apply_loss_bfs```: same as ```apply_loss``` but it uses breadth-first graph traversal (```apply_loss``` uses the tree-data structure that is also used by the Newman-Ziff method).
- ```fusion_and_loss_no_newman_ziff```: apply fusions to a graph state (consisting of many unconnected subgraphs).
  - considers qubit/photon loss, and eventually determines the largest connected component or cluster spanning.
  - could be used for arbitrary resource states but the physics is only correct for star-shaped resource states.
- ```analyse_fusion_net```: takes a connected graph/lattice as input and assumes that it represents a fusion network of star-shaped resource states (as in ```fusion_and_loss_nz```, every edge represents a fusion, every node the central qubit of a star-shaped resource state). No Newman-Ziff method.

Some functions for constructing a graph can be found in ```src/graph construction.c```:
- ```get_2d_square```: square lattice graph in two dimensions.
- ```get_2d_triangular```: triangular lattice graph in two dimensions.
- ```get_tree_lattice```: graph that is a tree (e.g. the Bethe lattice)
- ```get_nd_simple_cubic```: graph that is an n-dimensional simple cubic lattice.
- ```get_nd_simple_cubic_block```: many blocks with a hypercubic lattice. Between the blocks, there is one layer of nodes which are only connected in the direction perpendicular to the two blocks next to them.
- ```get_nd_simple_cubic_modified```: graph that is an n-dimensional lattice including simple cubic, diamond, fcc, bcc
- ```get_raussendorf_lattice```: lattice from https://doi.org/10.1016/j.aop.2006.01.012
- ```get_5qubit_stars_for_2d_square```: 5-qubit star-shaped (mutually unconnected) graph states on a square lattice
- ```get_3qubit_GHZ_for_2d_square```: use to determine bond-percolation threshold from https://doi.org/10.1038/s41467-019-08948-x (Fig. 4)
- ```get_lattice_from_nd_simple_cubic_and_vectors```: create a lattice by specifying connection vectors between points on a hypercubic lattice.
- ```get_GHZ_stars_for_nd_simple_cubic```: graph consisting of many micro-graph-states. Those can be fused to a large graph that will form an n-dimensional simple cubic lattice (if all fusions work and there is no photon loss).

## percolation simulations: getting started

To get started, we suggest having a look at the various examples in the folder ```percolation_c```. The ```C``` implementation uses the build system ```meson``` and you can compile all examples with the following commands:

```
$ meson setup build --buildtype=release
$ cd build/
$ ninja
```

By selecting the option ```--buildtype=release```, the compiler will do optimizations which will reduce the running time. For debugging with tools such as ```cgdb``` or ```valgrind```, please do not select this option. Note that the outlined compilation procedure requires ```meson``` the compiler ```ninja``` to be present on your system. If you haven't installed it yet, you can do so by running (on a ubuntu distribution):

```
$ sudo apt install meson
$ sudo apt install ninja
```

If you make any changes/improvements, we suggest validating everything with the provided tests in the folder ```tests/``` e.g. this way:
```
$ cd build/
$ meson test
```
In addition, the continous integration will run the tests with ```valgrind``` as a wrapper.

When adding new features, you can also add new tests by simply adding the name of the tests to the file ```tests/meson.build```.

## percolation simulations: Python vs. C
The implementations in ```percolation_py``` and ```percolation_C``` are different in some aspects. In contrast to the ```Python``` implementation, the ```C``` implementation has the advantage that a (modified) Newman-Ziff algorithm can be used which makes the simulations much faster.
