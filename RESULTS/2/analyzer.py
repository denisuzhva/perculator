import numpy as np

nF = np.loadtxt('./data_nF_i.txt')
nF = np.delete(nF, 0)
nB = np.loadtxt('./data_nB_i.txt')
nB = np.delete(nB, 0)
pF = np.loadtxt('./data_pF_i.txt')
pF = np.delete(nF, 0)
pB = np.loadtxt('./data_pB_i.txt')
pB = np.delete(nB, 0)

#print(nF)

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

print(b_nn)
print(b_pp)