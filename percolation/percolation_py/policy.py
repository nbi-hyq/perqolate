import numpy as np

# policy deciding if a fusion is done or not (just some examples)
def do_fusion_policy(data_in, policy_hyperparam):
    num_success_n = data_in[0]
    num_photon_n = data_in[1]
    num_success_nb = data_in[2]
    num_photon_nb = data_in[3]
    edge_node_ratio = data_in[4] / data_in[5]

    if policy_hyperparam <= 1:
        return np.random.random(1) < policy_hyperparam
    if policy_hyperparam == 2:
        if edge_node_ratio > 0.5 or (num_success_n >= 2 and num_success_nb <= 1 and num_photon_nb < 5) or (num_success_nb >= 2 and num_success_n <= 1 and num_photon_n < 5):
            return False
    if policy_hyperparam == 3:
        if num_success_n >= 1 and num_success_nb >= 1:
            return False
    if policy_hyperparam == 4:
        if num_success_n >= 2 or num_success_nb >= 2:
            return False
    if policy_hyperparam == 5:
        if (num_success_n >= 4 and num_success_nb < 1 and num_photon_nb <= 2) or (num_success_nb >= 4 and num_success_n < 1 and num_photon_n <= 2):
            return False
    if policy_hyperparam == 6:
        if num_success_n >= 3 and num_success_nb <= 1 and num_photon_nb < 2 or num_success_nb >= 3 and num_success_n <= 1 and num_photon_n < 2:
            return False
    return True
