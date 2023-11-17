# graphHyQ
- ```C```-code for percolation simulations of graphs, graph states, fusion networks. The main part of the code is in ```percolation/percolation_c/```.
- Some simulation results can be found here: https://arxiv.org/abs/2304.03796
- Add-on: display and transform quantum states such as graph states and hypergraph states (folder ```quantum_states/```).

## Percolation simulations: features
- Percolation simulations on undirected graphs/lattices. Percolation models: site-, bond-, and site-bond, bond-site
    - fast speed due to the use of the union-find algorithm by Newman and Ziff (https://doi.org/10.1103/PhysRevE.64.016706)
- Determine cluster spanning (percolation) or largest connected component size of graphs
- Percolation simulation of graph states subject to loss. (speed-up with modified Newman-Ziff algorithm)
- Simulation of fusion networks where resource states are fused by probabilistic Bell measurements. Photon loss is simulated with the model from https://arxiv.org/abs/2304.03796 (speed-up with modified Newman-Ziff algorithm)
- Simulation of fusion network using repeat-until-success method to boost the fusion probability. Photon loss can be simulated. (speed-up with modified Newman-Ziff algorithm)
- Tests and continuous integration
- Examples and easy compilation via meson
