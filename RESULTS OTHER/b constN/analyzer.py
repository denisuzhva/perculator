import numpy as np
import math
import matplotlib.pyplot as plt


def plotter(x, y, y_arr, x_arr, title, patches):
    plt.figure(figsize=(9,6))
    plt.errorbar(x, y, y_arr, x_arr, marker='.', color='k', capsize=2, elinewidth=0.5, markeredgewidth=0.5, linewidth=0)
    plt.xlabel(r'$<\eta>$', fontsize=14)
    plt.ylabel(r'$b_{nn}$', fontsize=14)
    #plt.yscale('log')
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.grid(b=True, which='minor', color='#999999', linestyle='--', alpha=0.2)
    plt.axis([0, 12.0, -0.02, 0.1])
    plt.title(title)
    plt.legend(handles=patches, fontsize=20, loc=(0.61, 0.71))
    plt.show()


if __name__ == "__main__":
    eta = np.loadtxt('./eta.txt')
    b = np.loadtxt('./b.txt')
    eta_err = np.zeros(eta.shape[0])
    b_err = eta_err
    patches = [plt.plot([],[], marker=".", ms=10, ls="", mec=None, color='k', label=r'$b_{nn}(<\eta>)$')[0]]
    plotter(eta, b, b_err, eta_err, r'correlations at $N = <N>$', patches)