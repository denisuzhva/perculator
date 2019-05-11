import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def processDataset(title):
    nn = 7
    data_tens = np.zeros((4, nn, 5001))
    for verse_iter in list(range(4)):
        #print('verse_iter: %i' % verse_iter)
        filename_read = './1_run/{0}_verse/data_{1}_i.txt'.format(verse_iter + 1, title)
        dataset_file = open(filename_read, 'r')
        dataset_lines = dataset_file.readlines()
        dataset_file.close()
        if nn == 6:
            del dataset_lines[-1]
        del dataset_lines[0]
        for i in list(range(len(dataset_lines))):
            data_tens[verse_iter, i, :] = np.asarray([float(x) for x in dataset_lines[i].split()])
    data_tens = data_tens[:, :, 1:]
    new_data_tens = np.zeros([nn, 20000])
    for i in list(range(nn)):
        new_data_tens[i, :] = np.reshape(data_tens[:, i, :], -1)
    
    filename_save = './1_run/total/data_{}_i.txt'.format(title)
    np.savetxt(filename_save, new_data_tens)


def makeTotal():
    dataset_title = ['pF', 'pB']
    for title in dataset_title:
        processDataset(title)


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
    makeTotal() # EXECUTE ONLY ONCE

    pF = np.loadtxt('./1_run/total/data_pF_i.txt')
    pB = np.loadtxt('./1_run/total/data_pB_i.txt')
    pFpB_av = np.average((pF * pB), axis=1)
    pF_av = np.average((pF), axis=1)
    pB_av = np.average((pB), axis=1)
    pF_av2 = pF_av*pF_av
    pF2_av = np.average((pF * pF), axis=1)
    b_pp = (pFpB_av - pF_av*pB_av) / (pF2_av - pF_av2)

    R = 7.5
    rs = 0.225
    S_0 = math.pi*R*R
    stringSigma = math.pi*rs*rs
    gamma_squared = 0.5
    nn_arr = np.array([100.0, 250.0, 550.0, 1100.0, 3300.0, 5600.0, 12200.0])
    eta = nn_arr * stringSigma / S_0
    x_fun = np.arange(0, np.max(eta), 0.01)
    y_fun_pp = b_pp_fun(x_fun, gamma_squared)

    plot_b(eta, b_pp, x_fun, y_fun_pp, 'b_pp(eta)', '1/(1+16*gamma_squared*sqrt(eta))')