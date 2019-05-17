import numpy as np
import math
import matplotlib.pyplot as plt


def b_nn_fun(x_arr):
    y_arr = 1/(1+4*np.sqrt(x_arr))
    return y_arr

def b_pp_fun(x_arr, gamma_sq):
    y_arr = 1/(1+16*gamma_sq*np.sqrt(x_arr))
    return y_arr

def plotter(eta, b_f, b_nf, x_fun, y_fun, y_max, title, formula, ylabel):
    plt.figure(figsize=(9,6))
    plt.plot(eta, b_f, 'v', color='k')
    plt.plot(eta, b_nf, '^', color='k')
    plt.plot(x_fun, y_fun, color='k', linewidth=0.5)
    plt.xlabel(r'$<\eta>$', fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    plt.yscale('log')
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.grid(b=True, which='minor', color='#999999', linestyle='--', alpha=0.2)
    plt.axis([0, 12.0, 0, y_max])
    plt.title(title)
    #f_patch = mpatches.Patch(marker='v', label='fusion')
    #nf_patch = mpatches.Patch(marker='^', label='nofusion')
    #calc_patch = mpatches.Patch(marker='-', label='asymptotics')
    patches = [plt.plot([],[], marker="v", ms=10, ls="", mec=None, color='k', label='fusion')[0],
               plt.plot([],[], marker="^", ms=10, ls="", mec=None, color='k', label='no-fusion')[0],
               plt.plot([],[], marker="_", ms=10, ls="", mec=None, color='k', label=formula)[0],]
    plt.legend(handles=patches, fontsize=18, loc=(0.41, 0.05))
    #plt.legend(handles=patches, fontsize=18, loc='best')
    plt.show()


if __name__ == "__main__":
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

    nF_f = np.loadtxt('./fusion/data_nF_i_f.txt')
    nB_f = np.loadtxt('./fusion/data_nB_i_f.txt')
    pF_f = np.loadtxt('./fusion/data_pF_i_f.txt')
    pB_f = np.loadtxt('./fusion/data_pB_i_f.txt')

    nF_nf = np.loadtxt('./nofusion/data_nF_i_nf.txt')
    nB_nf = np.loadtxt('./nofusion/data_nB_i_nf.txt')
    pF_nf = np.loadtxt('./nofusion/data_pF_i_nf.txt')
    pB_nf = np.loadtxt('./nofusion/data_pB_i_nf.txt')

    nFnB_av_f = np.average((nF_f * nB_f), axis=1)
    nF_av_f = np.average(nF_f, axis=1)
    nB_av_f = np.average(nB_f, axis=1)
    nF_av2_f = nF_av_f*nF_av_f
    nF2_av_f = np.average((nF_f * nF_f), axis=1)

    pFpB_av_f = np.average((pF_f * pB_f), axis=1)
    pF_av_f = np.average(pF_f, axis=1)
    pB_av_f = np.average(pB_f, axis=1)
    pF_av2_f = pF_av_f*pF_av_f
    pF2_av_f = np.average((pF_f * pF_f), axis=1)

    nFnB_av_nf = np.average((nF_nf * nB_nf), axis=1)
    nF_av_nf = np.average(nF_nf, axis=1)
    nB_av_nf = np.average(nB_nf, axis=1)
    nF_av2_nf = nF_av_nf*nF_av_nf
    nF2_av_nf = np.average((nF_nf * nF_nf), axis=1)

    pFpB_av_nf = np.average((pF_nf * pB_nf), axis=1)
    pF_av_nf = np.average(pF_nf, axis=1)
    pB_av_nf = np.average(pB_nf, axis=1)
    pF_av2_nf = pF_av_nf*pF_av_nf
    pF2_av_nf = np.average((pF_nf * pF_nf), axis=1)

    b_nn_f = (nFnB_av_f - nF_av_f*nB_av_f) / (nF2_av_f - nF_av2_f)
    b_nn_nf = (nFnB_av_nf - nF_av_nf*nB_av_nf) / (nF2_av_nf - nF_av2_nf)
    b_pp_f = (pFpB_av_f - pF_av_f*pB_av_f) / (pF2_av_f - pF_av2_f)
    b_pp_nf = (pFpB_av_nf - pF_av_nf*pB_av_nf) / (pF2_av_nf - pF_av2_nf)



    nn_formula = r'$\frac{1}{(1+4 \cdot \sqrt{\eta)}}$'
    pp_formula = r'$\frac{1}{(1+16 \cdot \gamma^2 \cdot \sqrt{\eta})} \quad \gamma^2 = 0.5)$'

    #plotter(eta, b_pp_f, b_pp_nf, x_fun, y_fun_pp, 1.1, r'$b_{p_tp_t}$ — fusion, no-fusion and asymptotics', pp_formula, r'$b_{p_tp_t}$')

    



    