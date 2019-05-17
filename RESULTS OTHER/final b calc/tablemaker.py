import numpy as np

b_nn_f = np.array([0.46764241, 0.44827848, 0.42119959, 0.36760529, 0.1733475, 0.11454716, 0.07851881])
b_nn_nf = np.array([0.49288867, 0.49497225, 0.50489103, 0.50545169, 0.49806821, 0.50702784, 0.50170181])
b_pp_f = np.array([-0.00777478,  0.00426224,  0.01042918,  0.02343103,  0.05499484,  0.04512132, 0.03323535])
b_pp_nf = np.array([0.97951202, 0.99866255, 0.99987742, 0.99998517, 0.99999951, 0.99999988, 0.99999999])


D_nn = np.zeros(7, dtype=float)
D_pp = np.zeros(7, dtype=float)

for i in list(range(7)):
    D_nn[i] = 1 - (min(b_nn_f[i], b_nn_nf[i]))/(max(b_nn_f[i], b_nn_nf[i]))
    D_pp[i] = 1 - (min(b_pp_f[i], b_pp_nf[i]))/(max(b_pp_f[i], b_pp_nf[i]))

D_nn_mean = np.mean(D_nn)
D_pp_mean = np.mean(D_pp)

print(D_nn_mean)
print(D_pp_mean)

#for number in b_pp_nf:
#    print('%.3f' % number)