import numpy as np
import os



if __name__ == "__main__":
    dir_path = './1_run/1_verse/'
    for filename in os.listdir(dir_path):
        if filename.endswith('i.txt'):
            filename = dir_path + filename
            f = open(filename)
            lines = f.readlines()
            print(filename)
            print(len(lines))

    