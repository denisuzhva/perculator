import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import math

def plotter(x_arr, y_arr, y_error, y_max, title):
    plt.figure(figsize=(9,6))
    plt.errorbar(x_arr, y_arr, y_error, marker='.', color='k', capsize=2, elinewidth=0.5, markeredgewidth=0.5, linewidth=0)
    plt.xlabel(r'$\eta$', fontsize=14)
    plt.ylabel(r'$<n_{fusion}^F> / <n_{nofusion}^F>$', fontsize=14)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.grid(b=True, which='minor', color='#999999', linestyle='--', alpha=0.2)
    plt.axis([0, 12.0, 0, y_max])
    plt.title(title)
    patch = plt.plot([],[], marker=".", ms=10, ls="", mec=None, color='k', label=r'$\frac{<n_{fusion}^F>}{<n_{nofusion}^F>}(\eta)$')[0]
    plt.legend(handles=[patch], fontsize=18)
    plt.show()


if __name__ == "__main__":

    R = 7.5
    rs = 0.225
    S_0 = math.pi*R*R
    stringSigma = math.pi*rs*rs
    gamma_sq = 0.5

    nn = np.array([100.0, 250.0, 550.0, 1100.0, 3300.0, 5600.0, 12200.0])
    eta = nn * stringSigma / S_0


    nF_fusion = np.loadtxt('./data_nF_fusion.txt')
    nF_nofusion = np.loadtxt('./data_nF_nofusion.txt')
    nF_nofusion = nF_nofusion[:, 1:]
    mulrat = nF_fusion / nF_nofusion
    mulrat_mean = np.mean(mulrat, axis=1)
    mulrat_error = np.std(mulrat, axis=1)
    mulrat_mean = mulrat[:, 0]
    

    plotter(eta, mulrat_mean, mulrat_error, 1.1, 'decrease of multiplicity due to the string fusion')