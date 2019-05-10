import numpy as np
import os

###########
# WARNING #
###########

# EXECUTE ONLY ONCE


def correct_pB():
    pB_new = np.zeros((4, 5000))
    for verse_iter in list(range(4)):
        filename = './1_run/{0}_verse/data_pB_test_i.txt'.format(verse_iter + 1)
        pB_bad_np = np.loadtxt(filename)
        pB_new[verse_iter, :] = pB_bad_np[80000:85000]
    
    pB_new = pB_new.reshape(-1)
    pB_new_list = pB_new.tolist()

    savefilename = './1_run/total/data_pB_test_i.txt'
    f = open(savefilename, 'w')
    correct_line = str(pB_new_list).strip('[]')
    correct_line = correct_line.replace(',', ' ')
    correct_line = correct_line.replace('\'', ' ')
    correct_line = correct_line.replace('[', ' ')
    correct_line = correct_line.replace(']', ' ')
    f.write(correct_line)
    f.close()
    

def correct_pF():
    pF_new = np.zeros((4, 5000))
    for verse_iter in list(range(4)):
        filename = './1_run/{0}_verse/data_pF_test_i.txt'.format(verse_iter + 1)
        file = open(filename, 'r')
        file_lines = file.readlines()
        file.close()
        correct_line = file_lines[10]
        pF_new[verse_iter, :] = np.fromstring(correct_line, dtype=float, sep='\t')[1:5001]

    pF_new.reshape(-1)
    pF_new_list = pF_new.tolist()

    savefilename = './1_run/total/data_pF_test_i.txt'    
    f = open(savefilename, 'w')
    correct_line = str(pF_new_list).strip('[]')
    correct_line = correct_line.replace(',', ' ')
    correct_line = correct_line.replace('\'', ' ')
    correct_line = correct_line.replace('[', ' ')
    correct_line = correct_line.replace(']', ' ')
    f.write(correct_line)
    file.close()

def clear_from_100_1100():
    for verse_iter in list(range(4)):
        dir_path = './1_run/{0}_verse/'.format(verse_iter + 1)
        for filename in os.listdir(dir_path):
            if filename.endswith('F_i.txt') or filename.endswith('B_i.txt'):
                filename_wdir = dir_path + filename
                f = open(filename_wdir, 'r')
                lines = f.readlines()
                f.close()
                #print(filename)
                #print(len(lines))
                if len(lines) == 1:
                    correct_line = lines[0]
                else:
                    correct_line = lines[-1]
                
                correct_line_list = correct_line.split()
                del correct_line_list[0]
                correct_line_list = correct_line_list[0:5000]
                correct_line = str(correct_line_list).strip('[]')
                correct_line = correct_line.replace(',', ' ')
                correct_line = correct_line.replace('\'', ' ')
                correct_line = correct_line.replace('[', ' ')
                correct_line = correct_line.replace(']', ' ')
                print(filename)
                savefilename = './1_run/total/{0}'.format(filename)
                f_save = open(savefilename, 'a+')
                f_save.write(correct_line)
                


if __name__ == "__main__":
    correct_pB()
    correct_pF()
    clear_from_100_1100()

    