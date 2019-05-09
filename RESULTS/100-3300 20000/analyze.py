import numpy as np

title = 'nF'
data_temp = np.loadtxt('./1_run/1_verse/data_%s_i.txt' % title)
print(data_temp.shape)