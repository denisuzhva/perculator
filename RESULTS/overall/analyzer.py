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

def plot_b(eta, b, x_fun, y_fun, title, formula):
    plt.plot(eta, b, 'o', color='orangered')
    plt.plot(x_fun, y_fun, color='royalblue')
    plt.xlabel('eta')
    plt.ylabel('b')
    plt.grid(True)
    plt.axis([0, 12.0, 0, 0.2])
    plt.title(title)
    gen_patch = mpatches.Patch(color='orangered', label='generated')
    cal_patch = mpatches.Patch(color='royalblue', label=formula)
    plt.legend(handles=[gen_patch, cal_patch])
    plt.show()


if __name__ == "__main__":
    b_nn = np.array([0.46393217, 0.4445169,  0.42124118, 0.37571875, 0.18053061914954044, 0.111852, 0.074625])
    b_pp = np.array([0.00563547, -0.00808659,  0.0012177,   0.02890763, 0.10029719786728043, 0.091229, 0.071035])
    b_pp_test = np.array([0.00325804, 0.00312066, 0.00944548, 0.00705924, 0.10236134212705346])

    R = 7.5
    rs = 0.225
    S_0 = math.pi*R*R
    stringSigma = math.pi*rs*rs
    nn_arr = np.array([100.0, 250.0, 550.0, 1100.0, 3300.0, 5600.0, 12200.0])
    eta = nn_arr * stringSigma / S_0
    gamma_sq = 0.5
    x_fun = np.arange(0, np.max(eta), 0.01)
    y_fun_nn = b_nn_fun(x_fun)
    y_fun_pp = b_pp_fun(x_fun, gamma_sq)


    #plot_b(eta, b_nn, x_fun, y_fun_nn, 'b_nn(eta)', '1/(1+4*sqrt(eta))')
    plot_b(eta, b_pp, x_fun, y_fun_pp, 'b_pp(eta)', '1/(1+16*gamma_squared*sqrt(eta))')
    
