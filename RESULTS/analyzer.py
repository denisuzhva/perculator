import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import math


def concatenator(dataset_name):
    dataset_5_25 = np.loadtxt("./5-25/{}.txt".format(dataset_name))
    dataset_50_12000 = np.loadtxt("./50-12000/{}.txt".format(dataset_name))
    dataset_50_12000_1 = np.loadtxt("./50-12000/{}1.txt".format(dataset_name))
    
    dataset_50_12000_concat = np.concatenate((dataset_50_12000, dataset_50_12000_1), axis=0)
    dataset_50_12000_concat = dataset_50_12000_concat[:, 1:]
    dataset_5_25_new = dataset_5_25[:, 1:] 
    dataset = np.concatenate((dataset_5_25_new, dataset_50_12000_concat), axis=1)

    return dataset


def plotter(x_arr, y_arr, x_error, y_error, y_max, title):
    plt.errorbar(x_arr, y_arr, x_error, y_error, '.', color='orangered')
    plt.xlabel(r'$\eta$')
    plt.ylabel(r'$<n_{fusion} / n_{nofusion}>$')
    plt.grid(True)
    plt.axis([0, 12.0, 0, y_max])
    plt.title(title)
    #gen_patch = mpatches.Patch(color='orangered', label=r'$<n_{fusion} / n_{nofusion}>(\eta)$')
    #plt.legend(gen_patch)
    plt.show()
    

if __name__ == "__main__":
    names = ['data_MulRatio', 'data_Eta']
    mulrat = concatenator(names[0])
    eta = concatenator(names[1])

    mulrat_mean = np.mean(mulrat, axis=0)
    mulrat_error = np.std(mulrat, axis=0)
    eta_mean = np.mean(eta, axis=0)
    eta_error = np.std(eta, axis=0)

    
    plotter(eta_mean, mulrat_mean, eta_error, mulrat_error, 1.1, '')


