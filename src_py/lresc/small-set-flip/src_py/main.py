import default_names
from print_erase import*
import decoder
from decoder import xerror_to_array, xerror_to_list
from read_result import Result, save_new_res
# import resource
import time
import numpy as np
import itertools


def find_best_1D_embedding(H, minimize_checks=False):
    # find best ordering of bit and checks to reduce lr edges
    # good thing the codes are small
    n = H.shape[1]
    m = H.shape[0]
    bits = np.zeros(n)
    checks = np.zeros(m)

    check_nbhd = []
    for check in range(m):
        check_nbhd.append(np.where(H[check])[0])

    best_per = None
    min_val = np.count_nonzero(H)
    best_lr_checks = None

    for per in itertools.permutations([f"b{i}" for i in range(n)] + [f"c{i}" for i in range(m)]):
        # probably no checks at the end (at least if all checks have at least two bits in them)
        # also add condition for two checks next to each other, shouldn't really happen probably
        if per <= per[::-1]:
            if (per[0][0] == "c") or (per[-1][0] == "c"): continue
            fail = 0
            for i in range(n+m-1):
                if (per[i][0] == per[(i+1)%(n+m)][0] == "c"):
                    fail = 1
                    break
            if fail: continue

            for i, node in enumerate(per):
                if node[0] == "b":
                    bits[int(node[1:])] = i
                else:
                    checks[int(node[1:])] = i

            num_lr = 0
            lr_checks = set()

            for check in range(m):
                for bit in check_nbhd[check]:
                    if abs(checks[check]-bits[bit]) != 1:
                        num_lr += 1
                        lr_checks |= {(check, bit)}
                        if minimize_checks:
                            break

            if num_lr < min_val:
                best_per = per
                min_val = num_lr
                best_lr_checks = lr_checks

    return best_per, min_val, best_lr_checks

def concatenate(H, rep_size):
    # rep = repetition_code(rep_size)
    n = H.shape[1]*rep_size
    k = H.shape[1]-np.linalg.matrix_rank(H)
    new_H = np.zeros((n-k,n), dtype=np.uint8)

    new_H[0:H.shape[0], 0:H.shape[1]] = H
    for i in range(H.shape[1]):
        for j in range(rep_size-1):
            new_H[H.shape[0]+(i*(rep_size-1)+j)][i+(j*H.shape[1])] = 1
            new_H[H.shape[0]+(i*(rep_size-1)+j)][i+((j+1)*H.shape[1])] = 1

    return new_H

def embed_1D(oH, rep_size, best_emb=None):
    H = concatenate(oH, rep_size)
    if best_emb == None: best_emb = find_best_1D_embedding(oH)
    lr_connections = best_emb[2]

    checks = np.zeros(H.shape[0], dtype=int)
    bits = np.zeros(H.shape[1], dtype=int)

    emb = np.empty(oH.shape[0]+oH.shape[1], dtype=object)
    for i in range(oH.shape[0]+oH.shape[1]):
        if best_emb[0][i][0] == "b":
            emb[i] = np.empty(2*rep_size-1, dtype=object)
            bits[int(best_emb[0][i][1:])] = i
        else:
            emb[i] = best_emb[0][i]

    for bit in range(oH.shape[1]):
        # place the repetition codes
        for i in range(2*rep_size-1):
            if (i % 2 == 0):
                emb[bits[bit]][i] = f"b{bit+((i//2)*oH.shape[1])}"
            else:
                emb[bits[bit]][i] = f"c{oH.shape[0]+((i-1)//2)+(bit*(rep_size-1))}"

    for i, node in enumerate(np.hstack(emb)):
        if node[0] == "b":
            bits[int(node[1:])] = i
        else:
            checks[int(node[1:])] = i

    emb = np.hstack(emb)
    # print(emb)
    # print(bits)
    # print(checks)
    # print()

    new_lr_connections = []
    for check in range(oH.shape[0]):
        # only need to weight balance the original checks
        print(f"c{check}: {np.where(oH[check])[0]}")

        for bit in np.where(oH[check])[0]:
            min_dist = len(emb)
            new_bit = bit
            for i in range(rep_size):
                tmp_dist = abs(checks[check] - bits[bit+(i*oH.shape[1])])
                if tmp_dist < min_dist:
                    min_dist = tmp_dist
                    new_bit = bit + (i*oH.shape[1])

            print(f"{bit}-->{new_bit}")
            if (check, bit) in lr_connections: new_lr_connections.append((check, new_bit))
            H[check][bit] = 0
            H[check][new_bit] = 1

    # print(H)
    return H, emb, new_lr_connections

rep_size = 3

# oH = np.array([[1,1,1]])
oH = np.array([ # [5,2,3]
    [1,0,0,1,0],
    [0,1,1,0,0],
    [0,1,0,1,1]
])
H, emb, lr_cons = embed_1D(oH, rep_size)
print(emb)
print(H)
print(lr_cons)

bit_nbhd = []
for bit in range(H.shape[1]):
    bit_nbhd.append(np.where(H[:,bit])[0])
check_nbhd = []
for check in range(H.shape[0]):
    check_nbhd.append(np.where(H[check])[0])

ccode = decoder.Classical_code(H.shape[1], H.shape[0], bit_nbhd, check_nbhd, "tmp")

no_mask = np.array([[True for c2 in range(H.shape[0])] for v1 in range(H.shape[1])])
mask = np.array([[True for c2 in range(H.shape[0])] for v1 in range(H.shape[1])])
for lr_con in lr_cons:
    mask[:, lr_con[0]] = False
    mask[lr_con[1]] = False
p_mask = 1 - (np.count_nonzero(mask) / (ccode.m * ccode.n))


###########################################################
# Simulation parameters
###########################################################
# Change this value for another algorithm (for example for the parallel version)
Ts = np.arange(10, 101, 10)
# Ts = [10]
no_runs = 1000000000
p = 0.001
###########################################################

res_file_name = f"../results/tmp.res"
start = time.time()
rs = []
logical2 = decoder.Logical2(ccode)

for i in range(no_runs):
    for T in Ts:
        vv_errors = np.zeros((ccode.n, ccode.n), dtype=np.uint8)
        cc_errors = np.zeros((ccode.m, ccode.m), dtype=np.uint8)

        for t in range(T):
            new_vv_error, new_cc_error = xerror_to_array(ccode, decoder.random_error(ccode, p))
            vv_errors ^= np.array(new_vv_error, dtype=np.uint8)
            cc_errors ^= np.array(new_cc_error, dtype=np.uint8)

            # if (t % 5 == 0):
            #     res, guessed_error = decoder.run_algo_qcode(ccode, xerror_to_list(ccode, (vv_errors, cc_errors)), no_mask, logical2)
            # else:
            #     res, guessed_error = decoder.run_algo_qcode(ccode, xerror_to_list(ccode, (vv_errors, cc_errors)), mask, logical2)
            res, guessed_error = decoder.run_algo_qcode(ccode, xerror_to_list(ccode, (vv_errors, cc_errors)), mask, logical2)

            new_vv_error, new_cc_error = xerror_to_array(ccode, guessed_error)
            vv_errors ^= np.array(new_vv_error, dtype=np.uint8)
            cc_errors ^= np.array(new_cc_error, dtype=np.uint8)

        new_vv_error, new_cc_error = xerror_to_array(ccode, decoder.random_error(ccode, p))
        vv_errors ^= np.array(new_vv_error, dtype=np.uint8)
        cc_errors ^= np.array(new_cc_error, dtype=np.uint8)

        res, guessed_error = decoder.run_algo_qcode(ccode, xerror_to_list(ccode, (vv_errors, cc_errors)), no_mask, logical2)

        if res == 2:
            r = Result(T,0,0,ccode.n,ccode.m,ccode.id,p,p_mask,1,0,1)
            # print(xerror_to_list(ccode, (vv_errors, cc_errors)))
            # print()
        else:
            r = Result(T,0,0,ccode.n,ccode.m,ccode.id,p,p_mask,1,res,0)
        rs.append(r)

    if (i % 100 == 0):
        save_new_res(res_file_name, rs)
        rs = []

stop = time.time()
print("Time taken: " + str(stop - start) + "\n")
