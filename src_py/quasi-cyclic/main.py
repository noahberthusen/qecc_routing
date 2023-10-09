import numpy as np
from numpy.linalg import matrix_power, matrix_rank
import matplotlib.pyplot as plt
import galois

import sys
sys.path.append("../../../../ldpc_masked2/src/")
from ldpc import bp_decoder
from ldpc.osd import bposd_decoder

from result import Result, save_new_res, res_to_line

GF = galois.GF(2)
codes = [
    [6,6,3,1,2,3,1,2],
    [15,3,9,1,2,0,2,7],
    [9,6,3,1,2,3,1,2],
    [12,6,3,1,2,3,1,2],
    [12,12,3,2,7,3,1,2],
    [30,6,9,1,2,3,25,26],
    [21,18,3,10,17,5,3,19]
]
code = codes[1]


# code = [30,6,19,2,3,1,24,11]
code = [15,3,11,1,0,0,14,13]


def cyclic_shift_matrix(l):
    arr = np.eye(l, dtype=int)
    return np.roll(arr, axis=1, shift=1)

ell = code[0]
m = code[1]

x = np.kron(cyclic_shift_matrix(ell), np.eye(m))
y = np.kron(np.eye(ell), cyclic_shift_matrix(m))

A1 = matrix_power(x, code[2])
A2 = matrix_power(y, code[3])
A3 = matrix_power(y, code[4])
A = ( A1 + A2 + A3 ) % 2

B1 = matrix_power(y, code[5])
B2 = matrix_power(x, code[6])
B3 = matrix_power(x, code[7])
B = ( B1 + B2 + B3 ) % 2

Hx = np.hstack([A, B]).astype(int)
Hz = np.hstack([B.T, A.T]).astype(int)

x_mask = np.zeros(Hx.shape[0]).astype(np.uint8)
# x_mask_inds = [0, 1, 2, 3, 4, 5, 18, 19, 20, 21, 22, 23, 36, 37, 38, 39, 40, 41, 66, 67, 68, 69, 70, 71, 84, 85, 86, 87, 88, 89, 132, 133, 134, 135, 136, 137, 150, 151, 152, 153, 154, 155]
x_mask_inds = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
x_mask[x_mask_inds] = 1

p = 0.001
Hxm = Hx[x_mask == 0]
bposd_dec_masked = bposd_decoder(
    Hx, # the parity check matrix
    Hxm, #the masked parity check matrix
    error_rate=p,
    channel_probs=[None], #assign error_rate to each qubit. This will override "error_rate" input variable
    max_iter=Hx.shape[1], #the maximum number of iterations for BP)
    bp_method="msl",
    ms_scaling_factor=0, #min sum scaling factor. If set to zero the variable scaling factor method is used
    osd_method="osd0", #the OSD method. Choose from:  1) "osd_e", "osd_cs", "osd0"
    #osd_order=min(Hx.shape[0],60) #the osd search depth
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
    osd_order=min(Hx.shape[0],60) #the osd search depth
)



num_iters = 3000000000
p_mask = np.count_nonzero(x_mask)/Hx.shape[0]
res_file_name = f"results/{ell}_{m}/bp_osd.res"
Ts = np.arange(10, 51, 10).astype(int)
rs = []

for i in range(num_iters):
    for T in Ts:
        # x_mask = (np.random.random(Hx.shape[0]) < p_mask).astype(np.uint8)
        error = np.zeros(Hx.shape[1]).astype(np.uint8)

        for j in range(T):
            error ^= (np.random.random(Hx.shape[1]) < p).astype(np.uint8)
            syndrome = (Hx @ error) % 2

            # bp_dec.decode(syndrome, x_mask)
            # bp_dec.decode(syndrome, np.zeros(Hx.shape[0]))
            # error ^= bp_dec.bp_decoding.astype(np.uint8)
            if (j % 5 == 0):
                # pass
                bposd_dec_masked.decode(syndrome, x_mask)
                error ^= bposd_dec_masked.osd0_decoding.astype(np.uint)
            else:
                bposd_dec_unmasked.decode(syndrome, np.zeros(Hx.shape[0]))
                error ^= bposd_dec_unmasked.osdw_decoding.astype(np.uint8)

            # bposd_dec_unmasked.decode(syndrome, np.zeros(Hx.shape[0]))
            # error ^= bposd_dec_unmasked.osdw_decoding.astype(np.uint8)
            # bposd_dec_masked.decode(syndrome, x_mask)
            # error ^= bposd_dec_masked.osd0_decoding.astype(np.uint)

        error ^= (np.random.random(Hx.shape[1]) < p).astype(np.uint8)
        syndrome = (Hx @ error) % 2

        # bp_dec.decode(syndrome, np.zeros(Hx.shape[0]))
        # error ^= bp_dec.bp_decoding.astype(np.uint8)
        bposd_dec_unmasked.decode(syndrome, np.zeros(Hx.shape[0]))
        error ^= bposd_dec_unmasked.osdw_decoding.astype(np.uint8)

        res = Result(T, m, ell,
                        code[2], code[3], code[4], code[5], code[6], code[7],
                        0,0,0,0,
                        p, p_mask, 1, int(not np.any(error)))
        rs.append(res)

    if (i % 100 == 0):
        save_new_res(res_file_name, rs)
        rs = []
