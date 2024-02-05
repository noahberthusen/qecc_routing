import numpy as np
import os

full_path = os.path.realpath(__file__)
path, filename = os.path.split(full_path)

def read_code(file_name):
    file = open(file_name, 'r')
    tmp_r = file.readline().strip(',\n').split(',')
    tmp_s = file.readline().strip(',\n').split(',')
    tmp_n = file.readline().strip(',\n').split(',')
    tmp_m1 = file.readline().strip(',\n').split(',')
    tmp_m2 = file.readline().strip(',\n').split(',')
    if tmp_n[0] != 'n' or tmp_r[0] != 'r' or tmp_s[0] != 's' or tmp_m1[0] != 'm1' or tmp_m2[0] != 'm2':
        raise NameError('Bad file format1')
    r = int(tmp_r[1])
    s = int(tmp_s[1])
    n = int(tmp_n[1])
    m1 = int(tmp_m1[1])
    m2 = int(tmp_m2[1])

    H = np.zeros((m1+m2, n), dtype=int)

    for i in range(m1+m2):
        for j in file.readline().strip(',\n').split(','):
            H[i][int(j)] = 1

    params = { 'r': r, 's': s, 'n': n, 'm1': m1, 'm2': m2 }
    return (params, H)


if __name__ == "__main__":
    file = os.path.join(path, "./codes/4_5_60.txt")
    params, H = read_code(file)
    print(params)
    print(H.shape)
    for r in H:
        print(np.count_nonzero(r))