import numpy as np
import math

########
# 5600 #
########


def load_data(title):
    data = np.zeros((4, 5000), dtype=float)
    data_temp = np.loadtxt('./1_run/1_verse/data_%s_i.txt' % title)
    data[0] = np.delete(data_temp, 0)
    data_temp = np.loadtxt('./1_run/2_verse/data_%s_i.txt' % title)
    data[1] = np.delete(data_temp, 0)
    data_temp = np.loadtxt('./1_run/3_verse/data_%s_i.txt' % title)
    data[2] = np.delete(data_temp, 0)
    data_temp = np.loadtxt('./1_run/4_verse/data_%s_i.txt' % title)
    data[3] = np.delete(data_temp, 0)

    #data_temp = np.loadtxt('./2_run/1_verse/data_%s_i.txt' % title)
    #data[1, 0] = np.delete(data_temp, 0)
    #data_temp = np.loadtxt('./2_run/2_verse/data_%s_i.txt' % title)
    #data[1, 1] = np.delete(data_temp, 0)
    #data_temp = np.loadtxt('./2_run/3_verse/data_%s_i.txt' % title)
    #data[1, 2] = np.delete(data_temp, 0)
    #data_temp = np.loadtxt('./2_run/4_verse/data_%s_i.txt' % title)
    #data[1, 3] = np.delete(data_temp, 0)

    data_tot = data.reshape(-1)
    return data_tot

def b_nn_fun(x_arr):
    y_arr = 1/(1+4*np.sqrt(x_arr))
    return y_arr

def b_pp_fun(x_arr, gamma_sq):
    y_arr = 1/(1+16*gamma_sq*np.sqrt(x_arr))
    return y_arr

if __name__ == "__main__":
    titles = ['nF', 'nB', 'pF', 'pB']
    nF = load_data(titles[0])
    nB = load_data(titles[1])
    pF = load_data(titles[2])
    pB = load_data(titles[3])

    nFnB_av = np.average(nF * nB)
    nF_av = np.average(nF)
    nB_av = np.average(nB)
    nF_av2 = nF_av*nF_av
    nF2_av = np.average(nF * nF)

    pFpB_av = np.average(pF * pB)
    pF_av = np.average(pF)
    pB_av = np.average(pB)
    pF_av2 = pF_av*pF_av
    pF2_av = np.average(pF * pF)

    b_nn = (nFnB_av - nF_av*nB_av) / (nF2_av - nF_av2)
    b_pp = (pFpB_av - pF_av*pB_av) / (pF2_av - pF_av2)

    print('b_nn (meas): %f' % b_nn)
    print('b_pp (meas): %f' % b_pp)

    R = 7.5
    rs = 0.225
    S_0 = math.pi*R*R
    stringSigma = math.pi*rs*rs
    gamma_squared = 0.5

    nn = 5600.0
    eta = nn * stringSigma / S_0

    print('b_nn (calc): %f' % b_nn_fun(eta))
    print('b_pp (calc): %f' % b_pp_fun(eta, gamma_squared))