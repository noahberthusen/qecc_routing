import numpy as np
from numpy.linalg import matrix_power, matrix_rank
import matplotlib.pyplot as plt
import pandas as pd

import sys
sys.path.append("../../../../ldpc/src/")
from ldpc import bp_decoder
from ldpc.osd import bposd_decoder

from scipy.optimize import curve_fit
from result import Result, save_new_res, res_to_line

codes = [
    [6,6,3,1,2,3,1,2],
    [15,3,9,1,2,0,2,7],
    [9,6,3,1,2,3,1,2],
    [12,6,3,1,2,3,1,2],
    [12,12,3,2,7,3,1,2],
    [30,6,9,1,2,3,25,26],
    [21,18,3,10,17,5,3,19]
]
code = codes[0]
# code = [18,12,1,11,3,2,15,1]
# code = [28,14,26,6,8,7,9,20]
code = [4,4,3,1,2,3,1,2]

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
# H = np.vstack([Hx, Hz])

bp_dec = bp_decoder(
    Hx,
    error_rate=0.001,
    max_iter=Hx.shape[1],
    bp_method="msl"
)

bposd_dec = bposd_decoder(
    Hx,#the parity check matrix
    error_rate=0.001,
    channel_probs=[None], #assign error_rate to each qubit. This will override "error_rate" input variable
    max_iter=Hx.shape[1], #the maximum number of iterations for BP)
    bp_method="msl",
    ms_scaling_factor=0, #min sum scaling factor. If set to zero the variable scaling factor method is used
    osd_method="osd_cs", #the OSD method. Choose from:  1) "osd_e", "osd_cs", "osd0"
    osd_order=0 #the osd search depth
)

error = np.zeros(Hx.shape[1])
error[0] = 1
mask = np.zeros(Hx.shape[0])
mask[[2,4]] = 1

syndrome = (Hx @ error) % 2

# print(syndrome)
# print(mask)

bp_dec.decode(syndrome, mask)
print(bp_dec.bp_decoding)
print(bp_dec.converge)

# print(bposd_dec.decode(syndrome, mask))