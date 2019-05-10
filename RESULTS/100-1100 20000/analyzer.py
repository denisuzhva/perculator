import numpy as np
import os
import math


def clear_from_3300():
    dir_path = './1_run/1_verse/'
    for filename in os.listdir(dir_path):
        if filename.endswith('i.txt'):
            filename = dir_path + filename
            f = open(filename, 'r')
            lines = f.readlines()
            f.close()
            if len(lines) == 6:
                del lines[5]
            f = open(filename, 'w')
            f.write('\n')
            for line in lines:
                f.write(line)
            f.close()

def load_data(title):
    data = np.zeros((4, 20000), dtype=float)
    data_temp = np.loadtxt('./1_run/1_verse/data_%s_i.txt' % title)
    data = data_temp[:,1:]
    #data_temp = np.loadtxt('./1_run/2_verse/data_%s_i.txt' % title)
    #data[1] = data_temp[:, 1:]
    #data_temp = np.loadtxt('./1_run/3_verse/data_%s_i.txt' % title)
    #data[2] = data_temp[:, 1:]
    #data_temp = np.loadtxt('./1_run/4_verse/data_%s_i.txt' % title)
    #data[3] = data_temp[:, 1:]

    #data_temp = np.loadtxt('./2_run/1_verse/data_%s_i.txt' % title)
    #data[1, 0] = np.delete(data_temp, 0)
    #data_temp = np.loadtxt('./2_run/2_verse/data_%s_i.txt' % title)
    #data[1, 1] = np.delete(data_temp, 0)
    #data_temp = np.loadtxt('./2_run/3_verse/data_%s_i.txt' % title)
    #data[1, 2] = np.delete(data_temp, 0)
    #data_temp = np.loadtxt('./2_run/4_verse/data_%s_i.txt' % title)
    #data[1, 3] = np.delete(data_temp, 0)

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

    pFpB_av_test = np.average((pF_test * pB_test), axis=1)
    pF_av_test = np.average(pF_test, axis=1)
    pB_av_test = np.average(pB_test, axis=1)
    pF_av2_test = pF_av_test*pF_av_test
    pF2_av_test = np.average((pF_test * pF_test), axis=1)


    b_nn = (nFnB_av - nF_av*nB_av) / (nF2_av - nF_av2)
    b_pp = (pFpB_av - pF_av*pB_av) / (pF2_av - pF_av2)
    b_pp_test = (pFpB_av_test - pF_av_test*pB_av_test) / (pF2_av_test - pF_av2_test)



    R = 7.5
    rs = 0.225
    S_0 = math.pi*R*R
    stringSigma = math.pi*rs*rs
    gamma_squared = 0.5

    nn = np.array([100, 250, 550, 1100])
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

    