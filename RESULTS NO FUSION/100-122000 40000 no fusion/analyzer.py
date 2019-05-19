import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.patches as mpatches
import math



def b_nn_fun(x_arr):
    y_arr = 1/(1+4*np.sqrt(x_arr))
    return y_arr

def b_pp_fun(x_arr, gamma_sq):
    y_arr = 1/(1+16*gamma_sq*np.sqrt(x_arr))
    return y_arr

def plot_b(eta, b, x_fun, y_fun, y_max, title, formula):
    plt.plot(eta, b, 'o', color='orangered')
    plt.plot(x_fun, y_fun, color='royalblue')
    plt.xlabel('eta')
    plt.ylabel('b')
    plt.grid(True)
    plt.axis([0, 12.0, 0, y_max])
    plt.title(title)
    gen_patch = mpatches.Patch(color='orangered', label='generated')
    cal_patch = mpatches.Patch(color='royalblue', label=formula)
    plt.legend(handles=[gen_patch, cal_patch])
    plt.show()


if __name__ == "__main__":
    nF = np.loadtxt('./1_run/data_nF_i.txt')
    nB = np.loadtxt('./1_run/data_nB_i.txt')
    pF = np.loadtxt('./1_run/data_pF_i.txt')
    pB = np.loadtxt('./1_run/data_pB_i.txt')

    nF = nF[:, 1:]
    nB = nB[:, 1:]
    pF = pF[:, 1:]
    pB = pB[:, 1:]

    nFnB_av = np.average((nF * nB), axis=1)
    nF_av = np.average(nF, axis=1)
    nB_av = np.average(nB, axis=1)
    nF_av2 = nF_av*nF_av
    nF2_av = np.average((nF * nF), axis=1)

    pFpB_av = np.average((pF * pB), axis=1)
    pF_av = np.average(pF, axis=1)
    pB_av = np.average(pB, axis=1)
    pF_av2 = pF_av*pF_av
    pF2_av = np.average((pF * pF), axis=1)

    b_nn = (nFnB_av - nF_av*nB_av) / (nF2_av - nF_av2)
    b_pp = (pFpB_av - pF_av*pB_av) / (pF2_av - pF_av2)

    R = 7.5
    rs = 0.225
    S_0 = math.pi*R*R
    stringSigma = math.pi*rs*rs
    gamma_sq = 0.5

    nn = np.array([100.0, 250.0, 550.0, 1100.0, 3300.0, 5600.0, 12200.0])
    eta = nn * stringSigma / S_0
    x_fun = np.arange(0, np.max(eta), 0.01)
    y_fun_nn = b_nn_fun(x_fun)
    y_fun_pp = b_pp_fun(x_fun, gamma_sq)

    nn_formula = '1/(1+4*sqrt(eta))'
    pp_formula = '1/(1+16*(gamma^2)*sqrt(eta)) (gamma^2 = 0.5)'
    plot_b(eta, b_nn, x_fun, y_fun_nn, 1.0, 'b_nn(eta)', nn_formula)



