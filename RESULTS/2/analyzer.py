import numpy as np

nF = np.loadtxt('./data_nF_i.txt')
nF = np.delete(nF, 0)
nB = np.loadtxt('./data_nB_i.txt')
nB = np.delete(nB, 0)
#print(nF)

nFnB_av = np.average(nF * nB)
nF_av = np.average(nF)
nB_av = np.average(nB)
nF_av2 = nF_av*nF_av
nF2_av = np.average(nF * nF)

b_nn = (nFnB_av - nF_av*nB_av) / (nF2_av - nF_av2)

print(b_nn)