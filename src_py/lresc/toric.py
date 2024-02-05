import numpy as np
import sys
sys.path.append("../../../../ldpc_masked2/src/")
from ldpc.osd import bposd_decoder
from result import Result, save_new_res, res_to_line

def repetition_code(n):
    H = np.zeros((n-1,n))
    for i in range(n-1):
        H[i][i] = H[i][i+1] = 1
    return H

rep_size = 10
H = repetition_code(rep_size)

hx1 = np.kron(H, np.eye(H.shape[1], dtype=bool))
hx2 = np.kron(np.eye(H.shape[0], dtype=bool), H.T)
Hx = np.hstack([hx1, hx2])

hz1 = np.kron(np.eye(H.shape[1], dtype=bool), H)
hz2 = np.kron(H.T, np.eye(H.shape[0], dtype=bool))
Hz = np.hstack([hz1, hz2])

print(f"[[{Hx.shape[1]},{Hx.shape[1]-2*Hx.shape[0]},{rep_size}]]")

p = 0.001
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



num_iters = 1000000000
p_mask = 0
res_file_name = f"results/toric/bp_osd.res"
Ts = np.arange(10, 51, 10).astype(int)
rs = []

for i in range(num_iters):
    for T in Ts:
        mask = (np.random.random(Hx.shape[0]) < p_mask).astype(np.uint8)
        error = np.zeros(Hx.shape[1]).astype(np.uint8)

        for j in range(T):
            error ^= (np.random.random(Hx.shape[1]) < p).astype(np.uint8)
            syndrome = (Hx @ error) % 2

            # bp_dec.decode(syndrome, mask)
            # error ^= bp_dec.bp_decoding.astype(np.uint8)
            bposd_dec_unmasked.decode(syndrome, mask)
            error ^= bposd_dec_unmasked.osdw_decoding.astype(np.uint)

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