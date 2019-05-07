import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('bin/Release/data_Pos.txt')

plt.plot(data[:, 0], data[:, 1], '.')
plt.axes().set_aspect('equal', 'datalim')
plt.show()