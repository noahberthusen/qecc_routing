import numpy as np
from numpy.linalg import matrix_power, matrix_rank
import matplotlib.pyplot as plt
import galois

import sys
sys.path.append("../../../../ldpc_masked2/src/")
from ldpc import bp_decoder
from ldpc.osd import bposd_decoder
from pymatching import Matching

from result import Result, save_new_res, res_to_line
import itertools
from mec import make_circle


GF = galois.GF(2)
def find_distance(H):
    # good thing these codes are small
    n = H.shape[1]
    min_weight = n
    for i in range(2**n):
        cw = bin(i)[2:].zfill(n)
        cw = [int(digit) for digit in cw]
        if not np.any((H @ cw) % 2):
            weight = np.count_nonzero(cw)
            if 0 < weight < min_weight:
                min_weight = weight
    return min_weight

def repetition_code(n):
    H = np.zeros((n-1,n))
    for i in range(n-1):
        H[i][i] = H[i][i+1] = 1
    return H

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
                        lr_checks |= {check}
                        if minimize_checks:
                            break

            if num_lr < min_val:
                best_per = per
                min_val = num_lr
                best_lr_checks = lr_checks

    return best_per, min_val, best_lr_checks

def embed_1D(oH, rep_size, best_emb=None):
    H = concatenate(oH, rep_size)
    if best_emb == None: best_emb = find_best_1D_embedding(oH)[0]

    checks = np.zeros(H.shape[0], dtype=int)
    bits = np.zeros(H.shape[1], dtype=int)

    emb = np.empty(oH.shape[0]+oH.shape[1], dtype=object)
    for i in range(oH.shape[0]+oH.shape[1]):
        if best_emb[i][0] == "b":
            emb[i] = np.empty(2*rep_size-1, dtype=object)
            bits[int(best_emb[i][1:])] = i
        else:
            emb[i] = best_emb[i]

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
            H[check][bit] = 0
            H[check][new_bit] = 1

    # print(H)
    return H, emb

def get_relative_positions(H, emb):
    n = H.shape[1]
    m = H.shape[0]

    bits = np.zeros(n, dtype=int)
    checks = np.zeros(m, dtype=int)

    bit_nbhd = []
    for bit in range(n):
        bit_nbhd.append(np.where(H[:,bit])[0])
    check_nbhd = []
    for check in range(m):
        check_nbhd.append(np.where(H[check])[0])

    for i, node in enumerate(emb):
        if node[0] == "b":
            bits[int(node[1:])] = i
        else:
            checks[int(node[1:])] = i

    rel_positions = np.empty(n+m, dtype=object)

    for i, bit in enumerate(bits):
        rel_positions[bit] = [checks[check] - bit for check in bit_nbhd[i]]
    for i, check in enumerate(checks):
        rel_positions[check] = [bits[bit] - check for bit in check_nbhd[i]]

    return rel_positions

def embed_hgp(H, emb):
    lattice = np.empty((len(emb), len(emb)), dtype=object)
    rel_pos = get_relative_positions(H, emb)

    n = H.shape[1]
    k = n-np.linalg.matrix_rank(GF(H))
    qbts = np.array([None for _ in range(n**2+(n-k)**2)])

    qb_ct = 0
    x_ct = 0
    z_ct = 0
    for i in range(len(emb)):
        for j in range(len(emb)):
            if emb[i][0] == emb[j][0]:
                lattice[i][j] = f"q{qb_ct}"
                qbts[qb_ct] = (i,j)
                qb_ct += 1
            elif emb[i][0] == "c" and emb[j][0] == "b":
                lattice[i][j] = f"x{x_ct}"
                x_ct += 1
            else:
                lattice[i][j] = f"z{z_ct}"
                z_ct += 1

    Hz = np.zeros((z_ct, qb_ct), dtype=np.uint8)
    Hx = np.zeros((x_ct, qb_ct), dtype=np.uint8)

    for i, qbt in enumerate(qbts):
        y, x = qbt
        hor_nbrs = rel_pos[x]
        ver_nbrs = rel_pos[y]

        for nbr in hor_nbrs:
            gen = lattice[y][x+nbr]
            gen_type = gen[0]
            gen_ind = int(gen[1:])
            if gen_type == "z":
                Hz[gen_ind][i] = 1
            else:
                Hx[gen_ind][i] = 1
        for nbr in ver_nbrs:
            gen = lattice[y+nbr][x]
            gen_type = gen[0]
            gen_ind = int(gen[1:])
            if gen_type == "z":
                Hz[gen_ind][i] = 1
            else:
                Hx[gen_ind][i] = 1

    return Hz, Hx, qbts

oH = np.array([[1,1,1]])
# oH = np.array([
#     [1,0,1,0,0,0],
#     [0,1,0,1,0,0],
#     [1,1,0,0,1,0],
#     [0,0,0,0,1,1]
# ])
oH = np.array([
    [1,0,0,1,0],
    [0,1,1,0,0],
    [0,1,0,1,1]
])
rep_size = 2

H, emb = embed_1D(oH, rep_size)
Hz, Hx, qbts = embed_hgp(H, emb)
print(f"[[{len(qbts)},{(oH.shape[1]-np.linalg.matrix_rank(GF(oH)))**2},{rep_size*find_distance(oH)}]]")

x_mask = np.zeros(Hx.shape[0], dtype=np.uint8)
z_mask = np.zeros(Hz.shape[0], dtype=np.uint8)

rs = []
for i in range(Hx.shape[0]):
    gen_qbts = qbts[np.where(Hx[i])[0]]
    if make_circle(gen_qbts)[2] > 1:
        x_mask[i] = 1
for i in range(Hz.shape[0]):
    gen_qbts = qbts[np.where(Hz[i])[0]]
    if make_circle(gen_qbts)[2] > 1:
        z_mask[i] = 1


Hxm = Hx[x_mask == 0]

p = 0.005
bp_dec = bp_decoder(
    Hx,
    Hxm,
    error_rate=p,
    max_iter=Hx.shape[1],
    bp_method="msl",
    ms_scaling_factor=0
)

bposd_dec_masked = bposd_decoder(
    Hx, # the parity check matrix
    Hxm, #the masked parity check matrix
    error_rate=p,
    channel_probs=[None], #assign error_rate to each qubit. This will override "error_rate" input variable
    max_iter=Hx.shape[1], #the maximum number of iterations for BP)
    bp_method="msl",
    ms_scaling_factor=0, #min sum scaling factor. If set to zero the variable scaling factor method is used
    osd_method="osd0", #the OSD method. Choose from:  1) "osd_e", "osd_cs", "osd0"
    # osd_order=min(Hx.shape[0]-1,60) #the osd search depth
)

bposd_dec_unmasked = bposd_decoder(
    Hx, # the parity check matrix
    Hx, #the masked parity check matrix
    error_rate=p,
    channel_probs=[None], #assign error_rate to each qubit. This will override "error_rate" input variable
    max_iter=Hx.shape[1], #the maximum number of iterations for BP)
    bp_method="msl",
    ms_scaling_factor=0, #min sum scaling factor. If set to zero the variable scaling factor method is used
    osd_method="osd_cs", #the OSD method. Choose from:  1) "osd_e", "osd_cs", "osd0"
    osd_order=min(Hx.shape[0]-1,60) #the osd search depth
)



num_iters = 3000000000
p_mask = np.count_nonzero(x_mask) / Hx.shape[0]
res_file_name = f"results/q{Hx.shape[1]}/bp_osd.res"
Ts = np.arange(10, 101, 10).astype(int)
rs = []

for i in range(num_iters):
    for T in Ts:
        error = np.zeros(Hx.shape[1]).astype(np.uint8)

        for j in range(T):
            error ^= (np.random.random(Hx.shape[1]) < p).astype(np.uint8)
            syndrome = (Hx @ error) % 2

            # bp_dec.decode(syndrome, x_mask)
            # bp_dec.decode(syndrome, np.zeros(Hx.shape[0]))
            # error ^= bp_dec.bp_decoding.astype(np.uint8)
            # if (j % 2 == 0):
            #     pass
            #     bposd_dec_masked.decode(syndrome, x_mask)
            #     error ^= bposd_dec_masked.osdw_decoding.astype(np.uint)
            # else:
            #     bposd_dec_unmasked.decode(syndrome, np.zeros(Hx.shape[0]))
            #     error ^= bposd_dec_unmasked.osdw_decoding.astype(np.uint8)

            # bposd_dec_unmasked.decode(syndrome, np.zeros(Hx.shape[0]))
            # error ^= bposd_dec_unmasked.osdw_decoding.astype(np.uint8)
            bposd_dec_masked.decode(syndrome, x_mask)
            error ^= bposd_dec_masked.osdw_decoding.astype(np.uint)

        error ^= (np.random.random(Hx.shape[1]) < p).astype(np.uint8)
        syndrome = (Hx @ error) % 2

        # bp_dec.decode(syndrome, np.zeros(Hx.shape[0]))
        # error ^= bp_dec.bp_decoding.astype(np.uint8)
        bposd_dec_unmasked.decode(syndrome, np.zeros(Hx.shape[0]))
        error ^= bposd_dec_unmasked.osdw_decoding.astype(np.uint8)

        res = Result(T, Hx.shape[1], Hx.shape[1]-2*Hx.shape[0],
                        0,0,0,0,0,0,
                        0,0,0,0,
                        p, p_mask, 1, int(not np.any(error)))
        rs.append(res)

    if (i % 100 == 0):
        save_new_res(res_file_name, rs)
        rs = []


# print(Hx.shape[1])
# for i in range(10000000000):
#     error = np.array([1]*2 + [0]*(Hx.shape[1]-2))
#     np.random.shuffle(error)

#     syndrome = (Hx @ error) % 2

#     bposd_dec_unmasked.decode(syndrome, np.zeros(Hx.shape[0]))

#     if np.any(error ^ bposd_dec_unmasked.osdw_decoding.astype(np.uint8)):
#         print(np.where(error))


