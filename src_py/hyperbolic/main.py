from read_code import read_code
from result import Result, save_new_res
from pymatching import Matching
import numpy as np
from code_utils import par2gen, SGSOP
import uuid
import os

full_path = os.path.realpath(__file__)
path, filename = os.path.split(full_path)
file_name = os.path.join(path, "../codes/5_4_160.txt")

p_masks = [0.2]
# Ts = np.linspace(0,100,11,dtype=np.uint8)
# Ts = [10, 20, 30, 40, 50]
Ts = [200]
p0 = 0.001
no_test = 1


def decode_iterative(masked_matching, unmasked_matching, H, logicals, p, p_mask, T):
    noise = np.zeros(H.shape[1], dtype=np.uint8)
    mask = np.where(np.random.random(H.shape[0]) < p_mask)[0]
    partial_mask = np.random.choice(mask, int(0.1*len(mask)))
    partial_mask2 = np.random.choice(mask, int(0.1*len(partial_mask)))
    print(mask)
    print()
    print(partial_mask)
    print()
    print(partial_mask2)

    for t in range(1, T+1):
        noise = noise ^ (np.random.random(H.shape[1]) < p).astype(np.uint8)
        print(f"({np.count_nonzero(noise)},",end='')

        shots = (noise @ H.T) % 2

        
        if (t % 20) == 0:
            masked_matching.set_boundary_nodes(set(partial_mask2))
            shots[partial_mask2] = 0
            print("**",end='')
        elif (t % 5) == 0:
            masked_matching.set_boundary_nodes(set(partial_mask))
            shots[partial_mask] = 0
            print("*",end='')
        else:
            masked_matching.set_boundary_nodes(set(mask))
            shots[mask] = 0

        predicted_error = masked_matching.decode(shots)
        noise = noise ^ predicted_error
        print(f"{np.count_nonzero(noise)}),",end='')


    noise = noise ^ (np.random.random(H.shape[1]) < p).astype(np.uint8)
    print(f"{np.count_nonzero(noise)},",end='')

    shots = (noise @ H.T) % 2
    actual_observables = (noise @ logicals.T) % 2
    predicted_observables = unmasked_matching.decode(shots)

    print(np.any(predicted_observables != actual_observables))
    return np.any(predicted_observables != actual_observables)        
          


if __name__ == "__main__":
    
    params, H = read_code(file_name)

    n = H.shape[1]
    m = H.shape[0]
    m1 = params['m1']
    m2 = params['m2']
    r = params['r']
    s = params['s']
    Hx = H[:m1]
    Hz = H[m1:]

    Gx, col_Gx = par2gen(Hx)
    Gz, col_Gz = par2gen(Hz)
    logicals, generators = SGSOP(Gx, Gz, n)

    logX = np.array([l[1][n:] for l in logicals])
    logZ = np.array([l[0][:n] for l in logicals])

    unmasked_matching_X = Matching.from_check_matrix(Hx, faults_matrix=logX)
    masked_matching_X = Matching.from_check_matrix(Hx)
    # unmasked_matching_Z = Matching.from_check_matrix(Hz, faults_matrix=logZ)
    # masked_matching_Z = Matching.from_check_matrix(Hz)

    id = str(uuid.uuid4())[:8]
    res_file_name = os.path.join(path, "../results/py_" + id + ".res")

    for _ in range(no_test):
        for p_mask in p_masks:
            rs = []
            for t in Ts:
                res_X = decode_iterative(masked_matching_X, unmasked_matching_X, Hx, logX, p0, p_mask, t)
                # res_Z = decode_iterative(masked_matching_Z, unmasked_matching_Z, Hz, logZ, p0, p_mask, t)
                
                res = Result(t, r, s, n, p0, p_mask, 1, int(not res_X))
                rs.append(res)
            # save_new_res(res_file_name, rs)
