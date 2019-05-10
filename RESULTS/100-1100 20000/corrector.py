import numpy as np
#import os

###########
# WARNING #
###########

# EXECUTE ONLY ONCE


def correct_pB():
    filename = './1_run/1_verse/data_pB_test_i.txt'
    pB_bad_np = np.loadtxt(filename)
    pB_bad_np = pB_bad_np[0:80000]
    pB_new = np.reshape(pB_bad_np, (4, 20000))
    nn_arr = np.array([100, 250, 550, 1100])
    pB_new = np.insert(pB_new, 0, nn_arr, axis=1)
    np.savetxt(filename, pB_new, delimiter='\t')

def correct_pF():
    filename = './1_run/1_verse/data_pF_test_i.txt'
    file = open(filename, 'r')
    file_lines = file.readlines()
    file.close()
    file_lines_new = [file_lines[2], file_lines[4], file_lines[6], file_lines[8]]
    file = open(filename, 'w')
    for line in file_lines_new:
        file.write(line)
    file.close()


if __name__ == "__main__":
    correct_pF()

    