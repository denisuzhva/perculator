import numpy as np
import os
import math


def load_data(title):
    data = np.loadtxt('./1_run/total/data_%s_i.txt' % title)
    print(data.shape)
    return data

def b_nn_fun(x_arr):
    y_arr = 1/(1+4*np.sqrt(x_arr))
    return y_arr

def b_pp_fun(x_arr, gamma_sq):
    y_arr = 1/(1+16*gamma_sq*np.sqrt(x_arr))
    return y_arr


if __name__ == "__main__":
    titles = ['nF', 'nB', 'pF', 'pB', 'pF_test', 'pB_test']
    nF = load_data(titles[0])
    nB = load_data(titles[1])
    pF = load_data(titles[2])
    pB = load_data(titles[3])
    pF_test = load_data(titles[4])
    pB_test = load_data(titles[5])


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

    pFpB_av_test = np.average(pF_test * pB_test)
    pF_av_test = np.average(pF_test)
    pB_av_test = np.average(pB_test)
    pF_av2_test = pF_av_test*pF_av_test
    pF2_av_test = np.average(pF_test * pF_test)


    b_nn = (nFnB_av - nF_av*nB_av) / (nF2_av - nF_av2)
    b_pp = (pFpB_av - pF_av*pB_av) / (pF2_av - pF_av2)
    b_pp_test = (pFpB_av_test - pF_av_test*pB_av_test) / (pF2_av_test - pF_av2_test)



    R = 7.5
    rs = 0.225
    S_0 = math.pi*R*R
    stringSigma = math.pi*rs*rs
    gamma_squared = 0.5

    nn = 3300
    eta = nn * stringSigma / S_0

    print('b_nn (calc):')
    print(b_nn_fun(eta))
    print('b_nn (meas):')
    print(b_nn)
    print('b_pp (calc): %f')
    print(b_pp_fun(eta, gamma_squared))
    print('b_pp (meas):')
    print(b_pp)
    print('b_pp (meas TEST):')
    print(b_pp_test)

    