import numpy as np
from scipy import optimize as opt
import matplotlib.pyplot as plt


# for fitting percolation threshold towards infinity
def pc_scaling(size, f, pc_inf, d):  # size: lattice size, f: pre-factor, pc_inf: threshold at inf size, d: dimension
    l_v = [None, None, 4/3, 0.88, 0.68, 0.57, 0.5]  # D. Stauffer p.52 (resp. p. 60)
    v = l_v[round(d)]
    pc = pc_inf + f * size**(-1/v)
    return pc


# alternative for fitting percolation threshold towards infinity
def pc_scaling_alt(pc_std, f, pc_inf):  # pc_std: standard deviation of pc(L), f: pre-factor, pc_inf: threshold at inf
    pc = pc_inf + f * pc_std
    return pc


# reshape and average C simulations and get the corresponding percolation thresholds
def get_thresholds(path_name, file_name, l_size, idx_fusion, l_fusion, n_rep, sweep_len, dimension, plot=False):
    n_size = len(l_size)  # parameter sweep in outer loop
    x_list = []  # x-axis of averaged percolation simulation
    y_list = []  # y-axis of averaged percolation simulation
    thr_list = []  # thresholds from individual runs
    with open(path_name + file_name, "r") as file:
        for line in file:
            h = line.split(" ")
            if len(h) == 2:
                x_list.append(float(h[0]))
                y_list.append(float(h[1]))
            elif len(h) == 1:
                thr_list.append(float(h[0]))
    x_array = np.array(x_list[0:sweep_len])
    y_array = np.array(y_list)
    y_array = np.reshape(y_array, (n_size, len(l_fusion), sweep_len))  # order given by loop nesting in C
    y_array = y_array[:, idx_fusion, :]  # take only part of the full simulation
    if len(thr_list) > 0:  # there are no thresholds when e.g. just plotting largest component size
        thr_array = np.reshape(thr_list, (n_size, len(l_fusion), n_rep))  # order given by loop nesting in C
        thr_array = thr_array[:, idx_fusion, :]  # take only part of the full simulation
    thr_mean = []  # mean value of the threshold
    thr_error = []  # error on the threshold
    n_cutoff_below = np.zeros(n_size)  # keep track if percolation threshold is too close to cutoff
    n_cutoff_above = np.zeros(n_size)  # keep track if percolation threshold is too close to cutoff
    for i in range(n_size):
        xrange = np.max(x_array) - np.min(x_array)
        if len(thr_list) > 0:
            n_cutoff_below[i] = np.sum(thr_array[i, :] < np.min(x_array) + xrange * 0.05)
            n_cutoff_above[i] = np.sum(thr_array[i, :] > np.max(x_array) - xrange * 0.05)
            thr_mean.append(np.mean(thr_array[i, :]))
            thr_error.append(np.std(thr_array[i, :]) / np.sqrt(n_rep - 1))
    if plot:
        fig1, ax1 = plt.subplots()
        for i in range(n_size):
            ax1.plot(x_array, y_array[i, :])
        ax1.set_xlabel('probability')
        ax1.set_ylabel('order parameter')
        ax1.legend(l_size)
        fig1.suptitle(file_name + ', # fusion: ' + str(l_fusion[idx_fusion]))
        plt.show()
        fig2, ax2 = plt.subplots()
        fig2.suptitle(file_name + ', # fusion: ' + str(l_fusion[idx_fusion]))
        if len(thr_list) > 0:  # only the case if the percolation transition is recorded
            for i in range(n_size):
                ax2.hist(thr_array[i, :])
            ax2.set_xlabel('probability')
            ax2.set_ylabel('frequency')
            plt.show()
        else:  # makes only sense when the largest connected component is simulated
            l_v = [None, None, 4 / 3, 0.88, 0.68, 0.57, 0.5]  # D. Stauffer p.52 (resp. p. 60)
            l_b = [None, None, 5 / 36, 0.41, 0.64, 0.84]  # D. Stauffer p.52 (resp. p. 60)
            for i in range(n_size):
                ax2.plot(x_array, y_array[i, :] * l_size[i]**(l_b[dimension]/l_v[dimension]))
            ax2.set_xlabel('probability')
            ax2.set_ylabel('$P_{max}*L^{β/ν}$')
            ax2.legend(l_size)
            plt.show()

    return np.array(thr_mean), np.array(thr_error), n_cutoff_below, n_cutoff_above


# fitting to estimate the percolation threshold at infinite size (see: https://doi.org/10.1142/S0129183198000431)
def get_infinity_threshold(list_l, avg_thr, err_thr, d, title='get infinity threshold', plot=True):
    p_start = np.array([1, avg_thr[np.argmax(list_l)]])
    p_fit, p_cov = opt.curve_fit(lambda l, f, pc_inf: pc_scaling(l, f, pc_inf, d), list_l, avg_thr, p0=p_start, sigma=err_thr)
    p_fit_b, p_cov_b = opt.curve_fit(lambda l, f, pc_inf: pc_scaling(l, f, pc_inf, d), list_l[-2:], avg_thr[-2:], p0=p_start, sigma=err_thr[-2:])  # fit only to last two data points
    p_fit_alt, p_cov_alt = opt.curve_fit(pc_scaling_alt, err_thr, avg_thr, p0=p_start)  # alternative (see eq. 6), (arbitrary errors for now TBD)
    if plot:
        plt.errorbar(list_l, avg_thr, err_thr)
        l_fit_plt = np.array(list(range(min(list_l), max(list_l) + 1)))
        plt.plot(l_fit_plt, [pc_scaling(size, p_fit[0], p_fit[1], d) for size in l_fit_plt])
        plt.plot(l_fit_plt, [pc_scaling(size, p_fit_b[0], p_fit_b[1], d) for size in l_fit_plt])
        plt.xlabel('L (box size)')
        plt.ylabel('p_c')
        plt.title(title)
        plt.show()
        plt.plot(err_thr, avg_thr, 'ob')
        err_fit_plt = np.array([min(err_thr), max(err_thr)])
        plt.plot(err_fit_plt, pc_scaling_alt(err_fit_plt, p_fit_alt[0], p_fit_alt[1]), '-r')
        plt.xlabel('std(L)')
        plt.ylabel('p_c')
        plt.title(title)
        plt.show()
    return p_fit, p_cov, p_fit_b, p_cov_b, p_fit_alt, p_cov_alt


def run_analysis(path, l_name, l_dim, l_l_size, l_num_fusion, l_num_rep, l_num_sweep):
    l_thresh = [[] for _ in range(len(l_num_fusion))]  # percolation thresholds
    l_thresh_alt = [[] for _ in range(len(l_num_fusion))]  # thresholds obtained with alternative method
    l_err_stat = [[] for _ in range(len(l_num_fusion))]  # statistical error
    l_err_syst = [[] for _ in range(len(l_num_fusion))]  # estimate for systematic error
    for idx in range(len(l_num_fusion)):
        for k in range(len(l_name)):
            avg_threshold, err_threshold, cutoff_below, cutoff_above = get_thresholds(path, l_name[k], l_l_size[k], idx, l_num_fusion, l_num_rep[k], l_num_sweep[k], l_dim[k], plot=plot_get_threshold)
            param_fit, param_cov, param_fit_b, param_cov_b, param_fit_alt, param_cov_alt = get_infinity_threshold(l_l_size[k], avg_threshold, err_threshold, l_dim[k], title=l_name[k]+', fusion # = '+str(l_num_fusion[idx]), plot=plot_inf_fit)
            param_err = np.sqrt(np.diag(param_cov))  # just look at diagonal of covariance matrix
            param_err_b = np.sqrt(np.diag(param_cov_b))  # just look at diagonal of covariance matrix
            l_thresh[idx].append(param_fit[1])
            l_thresh_alt[idx].append(param_fit_alt[1])
            l_err_stat[idx].append(param_err[1])
            l_err_syst[idx].append(abs(param_fit[1] - param_fit_b[1]))
            print(l_name[k], ' np. fusion attempts', l_num_fusion[idx], param_fit[1], ' +-(syst) ', l_err_syst[idx][-1],' +-(stat) ', l_err_stat[idx][-1])
            print('check: ', l_name[k], ' np. fusion attempts', l_num_fusion[idx], ' ', l_thresh_alt[idx][-1])
            if plot_cutoff:
                plt.plot(l_l_size[k], cutoff_above)
                plt.plot(l_l_size[k], cutoff_below)
                plt.xlabel('L (size)')
                plt.ylabel('# cutoff of ' + str(l_num_rep[k]))
                plt.title(l_name[k] + ', # fusion attempts: ' + str(l_num_fusion[idx]))
                plt.legend(['cutoff above', 'cutoff below'])
                plt.show()
    # quadratically add systematic and statistical error
    return l_thresh, l_thresh_alt, np.sqrt(np.array(l_err_stat) ** 2 + np.array(l_err_syst) ** 2)


if __name__ == '__main__':
    plot_get_threshold = False  # plot fits for getting percolation thresholds for finite size simulations
    plot_inf_fit = False  # plot extrapolation for getting percolation threshold at infinite size
    plot_cutoff = False  # make plot for checking if simulation ranges were cut off
    plot_largest_component = False

    # plot size of largest component
    if plot_largest_component:
        lname = ["size_nz_2d.txt", "size_nz_3d.txt", "size_nz_4d.txt", "size_nz_5d.txt", "size_nz_2d_static.txt", "size_nz_3d_static.txt", "size_nz_4d_static.txt", "size_nz_5d_static.txt"]  # file names
        llsize = [[i for i in range(500, 3001, 500)], [i for i in range(50, 201, 50)], [i for i in range(20, 57, 12)], [i for i in range(10, 26, 5)]]
        llsize.extend(llsize)
        l_dimension = [2, 3, 4, 5]  # dimension of lattices
        l_dimension.extend(l_dimension)
        for i in range(len(lname)):
            _, _, _, _ = get_thresholds("G:\\My Drive\\percolation_project\\simulations\\2023_01_04\\main_lattices\\sc\\", lname[i], llsize[i], 0, [1], 100, 1000, l_dimension[i], plot=True)

    # analyse percolation
    figs, ax = plt.subplots()
    spath = "G:\\My Drive\\percolation_project\\simulations\\2023_01_04\\main_lattices\\sc\\"  # path
    lname = ["percol_nz_2d.txt", "percol_nz_3d.txt", "percol_nz_4d.txt", "percol_nz_5d.txt"]  # file names
    ldimension = [2, 3, 4, 5]  # dimension of lattices
    llsize = [[i for i in range(500, 3001, 500)], [i for i in range(50, 201, 50)], [i for i in range(20, 57, 12)], [i for i in range(10, 26, 5)]]
    lnum_fusion = [1]  # number of fusion attempts
    lnum_rep = [100, 100, 100, 100]  # number of repetitions for averaging
    lnum_sweep = [1000, 1000, 1000, 1000]  # sweep of site/bond probability (most inner loop)
    ll_threshold, _, arry_err = run_analysis(spath, lname, ldimension, llsize, lnum_fusion, lnum_rep, lnum_sweep)
    for fusion_idx in range(len(lnum_fusion)):
        ax.errorbar(ldimension, ll_threshold[fusion_idx], arry_err[fusion_idx, :], label='$\lambda_c$, # fusion attempts: ' + str(lnum_fusion[fusion_idx]))

    ax.set_xlabel('dimension')
    ax.set_ylabel('percolation threshold $\lambda_c$')
    ax.set_ylim(top=1.0)
    ax.legend()
    plt.show()
