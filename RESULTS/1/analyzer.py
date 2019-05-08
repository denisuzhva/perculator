import numpy as np
import matplotlib.pyplot as plt
import math

def b_nn_fun(x_arr):
    y_arr = 1/(1+4*np.sqrt(x_arr))
    return y_arr

R = 7.5
rs = 0.225
S_0 = math.pi*R*R
stringSigma = math.pi*rs*rs
data = np.loadtxt('./data_b.txt')
b_nn = data[:, 1]

nn_arr = np.array([100.0, 250.0, 550.0, 1100.0, 3300.0, 5600.0, 12200.0])

eta = nn_arr * stringSigma / S_0
print(eta)

x_fun = np.arange(0, np.max(eta), 0.01)
y_fun = b_nn_fun(x_fun)

plt.plot(eta, b_nn, '.')
plt.plot(x_fun, y_fun)
#plt.axis([0, np.max(eta)+1, -1, 1])
plt.show()