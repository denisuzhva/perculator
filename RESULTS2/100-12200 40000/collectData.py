import numpy as np


def collect_p():
    path_pF_1 = './1_run/data_pF_i.txt'
    path_pF_2 = './2_run/data_pF_i.txt'
    path_pB_1 = './1_run/data_pB_i.txt'
    path_pB_2 = './2_run/data_pB_i.txt'
    pF = np.zeros((2, 7, 20000))
    pF[0, :, :] = np.loadtxt(path_pF_1)
    pF[1, :, :] = np.loadtxt(path_pF_2)
    pF = np.concatenate((pF[0, :, :], pF[1, :, :]), axis=1)

    pB = np.zeros((2, 7, 20000))
    pB[0, :, :] = np.loadtxt(path_pB_1)
    pB[1, :, :] = np.loadtxt(path_pB_2)
    pB = np.concatenate((pB[0, :, :], pB[1, :, :]), axis=1)

    np.savetxt('./data_pF_i_40000.txt', pF)
    np.savetxt('./data_pB_i_40000.txt', pB)
    

def collect_100_1100(dataset_name):
    dataset = np.loadtxt('../../RESULTS/100-1100 20000/1_run/1_verse/data_{}_i.txt'.format(dataset_name))
    dataset = dataset[:, 1:]
    return dataset


def collect_3300(dataset_name):
    dataset = np.loadtxt('../../RESULTS/3300 20000/1_run/total/data_{}_i.txt'.format(dataset_name))
    return dataset


def collect_5600(dataset_name):
    dataset = np.zeros((4, 5001))
    for verse_iter in list(range(4)):
        dataset[verse_iter] = np.loadtxt('../../RESULTS/5600 20000/1_run/{0}_verse/data_{1}_i.txt'.format(verse_iter+1, dataset_name))
    dataset = dataset[:, 1:]
    dataset = np.reshape(dataset, -1)
    return dataset


def collect_12200(dataset_name):
    dataset = np.zeros((2, 4, 2501))
    for run_iter in list(range(2)):
        for verse_iter in list(range(4)):
            dataset[run_iter, verse_iter, :] = np.loadtxt('../../RESULTS/12200 20000/{0}_run/{1}_verse/data_{2}_i.txt'.format(run_iter+1, 
                                                                                                                              verse_iter+1, 
                                                                                                                              dataset_name))
    dataset = dataset[:, :, 1:]
    dataset = np.reshape(dataset, -1)
    return dataset


def make_1run():
    nF = np.zeros((7, 20000))
    nB = np.zeros((7, 20000)) 
    nF[0:4, :] = collect_100_1100('nF')
    nF[4, :] = collect_3300('nF')
    nF[5, :] = collect_5600('nF')
    nF[6, :] = collect_12200('nF')

    nB[0:4, :] = collect_100_1100('nB')
    nB[4, :] = collect_3300('nB')
    nB[5, :] = collect_5600('nB')
    nB[6, :] = collect_12200('nB')
    np.savetxt('./1_run/data_nF_i.txt', nF)
    np.savetxt('./1_run/data_nB_i.txt', nB)


def make_2run():
    nF_old = np.zeros((4, 7, 5001))
    nB_old = np.zeros((4, 7, 5001))
    for verse_iter in list(range(4)):
        nF_old[verse_iter, :, :] = np.loadtxt('../../RESULTS/100-12200 20000 2/2_run/{}_verse/data_nF_i.txt'.format(verse_iter + 1))
        nB_old[verse_iter, :, :] = np.loadtxt('../../RESULTS/100-12200 20000 2/2_run/{}_verse/data_nB_i.txt'.format(verse_iter + 1))
    nF_old = nF_old[:, :, 1:]
    nB_old = nB_old[:, :, 1:]
    nF = np.concatenate((nF_old[0, :, :], nF_old[1, :, :], nF_old[2, :, :], nF_old[3, :, :]), axis=1)
    nB = np.concatenate((nB_old[0, :, :], nB_old[1, :, :], nB_old[2, :, :], nB_old[3, :, :]), axis=1)
    np.savetxt('./2_run/data_nF_i.txt', nF)
    np.savetxt('./2_run/data_nB_i.txt', nB)


def collect_n():
    make_1run()
    make_2run()
    
    path_nF_1 = './1_run/data_nF_i.txt'
    path_nF_2 = './2_run/data_nF_i.txt'
    path_nB_1 = './1_run/data_nB_i.txt'
    path_nB_2 = './2_run/data_nB_i.txt'
    nF = np.zeros((2, 7, 20000))
    nF[0, :, :] = np.loadtxt(path_nF_1)
    nF[1, :, :] = np.loadtxt(path_nF_2)
    nF = np.concatenate((nF[0, :, :], nF[1, :, :]), axis=1)

    nB = np.zeros((2, 7, 20000))
    nB[0, :, :] = np.loadtxt(path_nB_1)
    nB[1, :, :] = np.loadtxt(path_nB_2)
    nB = np.concatenate((nB[0, :, :], nB[1, :, :]), axis=1)

    np.savetxt('./data_nF_i_40000.txt', nF)
    np.savetxt('./data_nB_i_40000.txt', nB)



if __name__ == "__main__":
    #collect_p()
    collect_n()

