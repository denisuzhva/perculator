import numpy as np
import math
import matplotlib.pyplot as plt


def plotter(x, y, y_arr, x_arr, title, patches):
    plt.figure(figsize=(9,6))
    plt.errorbar(x, y, y_arr, x_arr, marker='.', color='k', capsize=2, elinewidth=0.5, markeredgewidth=0.5, linewidth=0)
    plt.xlabel(r'$\eta$', fontsize=14)
    plt.ylabel(r'$N_{cl}/N$', fontsize=14)
    #plt.yscale('log')
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.grid(b=True, which='minor', color='#999999', linestyle='--', alpha=0.2)
    plt.axis([0, 12.0, 0, 1.1])
    plt.title(title)
    plt.legend(handles=patches, fontsize=20, loc=(0.61, 0.71))
    plt.show()

def calc_NclN_matlab():
    NN = np.loadtxt('./matlab/NN_matlab.txt')
    eta_arr = np.loadtxt('./matlab/eta_NN_matlab.txt')
    eta_error = np.zeros(eta_arr.shape[0])
    NN_mean = NN[:, 0]
    NN_error = NN[:, 1]
    patches = [plt.plot([],[], marker=".", ms=10, ls="", mec=None, color='k', label=r'$\frac{N_{cl}}{N}(\eta)$')[0]]
    plotter(eta_arr, NN_mean, NN_error, eta_error, 'uniform distribution, 5 simulations', patches)


def calc_SclS_matlab():
    SS = np.loadtxt('./matlab/SS_matlab.txt')
    eta_arr = np.loadtxt('./matlab/eta_SS_matlab.txt')
    eta_error = np.zeros(eta_arr.shape[0])
    SS_mean = SS[0, :]
    SS_error = SS[1, :]
    patches = [plt.plot([],[], marker=".", ms=10, ls="", mec=None, color='k', label=r'$\frac{S_{cl}}{S}(\eta)$')[0]]
    plotter(eta_arr, SS_mean, SS_error, eta_error, 'uniform distribution, 60 simulations', patches)


def calc_NclN_c():
    NN = np.loadtxt('./c/NclN.txt')
    eta_arr = np.loadtxt('./c/eta.txt')
    NN_mean = NN[0, :]
    NN_error = NN[1, :]
    eta_mean = eta_arr[0, :]
    eta_error = eta_arr[1, :]
    patches = [plt.plot([],[], marker=".", ms=10, ls="", mec=None, color='k', label=r'$\frac{N_{cl}}{N}(\eta)$')[0]]
    plotter(eta_mean, NN_mean, NN_error, eta_error, 'parabolic distribution, 3638 simulations', patches)


def calc_SclS_c():
    SS = np.loadtxt('./c/SclS.txt')
    eta_arr = np.loadtxt('./c/eta.txt')
    SS_mean = SS[0, :]
    SS_error = SS[1, :]
    eta_mean = eta_arr[0, :]
    eta_error = eta_arr[1, :]
    patches = [plt.plot([],[], marker=".", ms=10, ls="", mec=None, color='k', label=r'$\frac{S_{cl}}{S}(\eta)$')[0]]
    plotter(eta_mean, SS_mean, SS_error, eta_error, 'parabolic distribution, 3638 simulations', patches)


if __name__ == "__main__":

    #nn = np.array([100.0, 250.0, 550.0, 1100.0, 3300.0, 5600.0, 12200.0])
    #eta = nn * stringSigma / S_0

    calc_SclS_c()

