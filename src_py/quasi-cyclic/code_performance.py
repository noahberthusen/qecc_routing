import numpy as np
from numpy.linalg import matrix_power, matrix_rank
import itertools
import galois
import pandas as pd
from mec import make_circle
from scipy.sparse import csr_matrix
from scipy.stats import skew

import sys
sys.path.append("../../../../ldpc_masked2/src/")
from ldpc import bp_decoder
from ldpc.osd import bposd_decoder
from result import Result, save_new_res, res_to_line


def test_code(code):
    def cyclic_shift_matrix(l):
        arr = np.eye(l, dtype=int)
        return np.roll(arr, axis=1, shift=1)

    ell = int(code["ell"])
    m = int(code["m"])

    x = np.kron(cyclic_shift_matrix(ell), np.eye(m))
    y = np.kron(np.eye(ell), cyclic_shift_matrix(m))

    A1 = matrix_power(x, int(code["a1"]))
    A2 = matrix_power(y, int(code["a2"]))
    A3 = matrix_power(y, int(code["a3"]))
    A = ( A1 + A2 + A3 ) % 2

    B1 = matrix_power(y, int(code["b1"]))
    B2 = matrix_power(x, int(code["b2"]))
    B3 = matrix_power(x, int(code["b3"]))
    B = ( B1 + B2 + B3 ) % 2

    Hx = np.hstack([A, B]).astype(int)
    # Hz = np.hstack([B.T, A.T]).astype(int)

    def get_mask(code, init):
        emb_m, emb_ell, A_ind, B_ind = code

        lattice = np.empty((2*emb_m, 2*emb_ell), dtype=object)
        lattice[0][0] = f"x{init}"

        # As = [[A1, A2.T], [A2, A3.T], [A1, A3.T]]
        # Bs = [[B1, B2.T], [B2, B3.T], [B1, B3.T]]
        As = [[A1, A2.T], [A2, A1.T], [A2, A3.T], [A3, A2.T], [A1, A3.T], [A3, A1.T]]
        Bs = [[B1, B2.T], [B2, B1.T], [B2, B3.T], [B3, B2.T], [B1, B3.T], [B3, B1.T]]

        def get_nbr(i, j):
            if (i % 2 == 0):
                if (j % 2 == 0):
                    return "x"
                else:
                    return "r"
            else:
                if (j % 2 == 0):
                    return "l"
                else:
                    return "z"

        for i in range(2*emb_m - 1):
            for j in range(2*emb_ell):
                curr_ind = int(lattice[i][j][1:])

                if (i % 2 == 0):
                    tmp_A = As[A_ind][1]
                else:
                    tmp_A = As[A_ind][0]
                if (j % 2 == 0):
                    tmp_B = Bs[B_ind][1]
                else:
                    tmp_B = Bs[B_ind][0]

                lattice[(i+1)%(2*emb_m)][j] = f"{get_nbr((i+1)%(2*emb_m), j)}{np.where(tmp_A @ np.eye(m*ell)[curr_ind])[0][0]}"
                lattice[i][(j+1)%(2*emb_ell)] = f"{get_nbr(i, (j+1)%(2*emb_ell))}{np.where(tmp_B @ np.eye(m*ell)[curr_ind])[0][0]}"

        for i in range(2*emb_m):
            for j in range(2*emb_ell):
                if (lattice[i][j][0] == "z"):
                    lattice[i][j] = f"z{int(lattice[i][j][1:]) + m*ell}"
                elif (lattice[i][j][0] == "r"):
                    lattice[i][j] = f"r{int(lattice[i][j][1:]) + m*ell}"


        qbts = np.array([None for i in range(2*m*ell)])
        for i in range(lattice.shape[0]):
            for j in range(lattice.shape[1]):
                if lattice[i][j][0] == "r" or lattice[i][j][0] == "l":
                    qbts[int(lattice[i][j][1:])] = (i, j)

        x_rs = []
        for i in range(m*ell):
            gen_qbts = qbts[np.where(Hx[i])[0]]
            x_rs.append(make_circle(gen_qbts)[2])
        # for i in range(m*ell):
        #     gen_qbts = qbts[np.where(Hz[i])[0]]
        #     z_rs.append(make_circle(gen_qbts)[2])

        x_mask = np.zeros(Hx.shape[0])
        for i, x in enumerate(x_rs):
            if (x > (max(x_rs)-min(x_rs))/2):
                x_mask[i] = 1

        return x_mask

    mask = get_mask((int(code["emb_m"]),int(code["emb_ell"]),int(code["aa"]),int(code["bb"])), 0)
    Hxm = Hx[mask == 0]

    bposd_dec_masked = bposd_decoder(
        Hx, # the parity check matrix
        Hxm, #the masked parity check matrix
        error_rate=0.005,
        channel_probs=[None], #assign error_rate to each qubit. This will override "error_rate" input variable
        max_iter=Hx.shape[1], #the maximum number of iterations for BP)
        bp_method="msl",
        ms_scaling_factor=0, #min sum scaling factor. If set to zero the variable scaling factor method is used
        osd_method="osd0", #the OSD method. Choose from:  1) "osd_e", "osd_cs", "osd0"
        osd_order=min(Hx.shape[0],60) #the osd search depth
    )

    bposd_dec_unmasked = bposd_decoder(
        Hx, # the parity check matrix
        Hx, #the masked parity check matrix
        error_rate=0.005,
        channel_probs=[None], #assign error_rate to each qubit. This will override "error_rate" input variable
        max_iter=Hx.shape[1], #the maximum number of iterations for BP)
        bp_method="msl",
        ms_scaling_factor=0, #min sum scaling factor. If set to zero the variable scaling factor method is used
        osd_method="osd_cs", #the OSD method. Choose from:  1) "osd_e", "osd_cs", "osd0"
        osd_order=min(Hx.shape[0],60) #the osd search depth
    )


    num_iters = 1001
    p = 0.005
    p_mask = np.count_nonzero(mask)/Hx.shape[0]
    res_file_name = f"results/q{Hx.shape[1]}/combined_results2.res"
    Ts = [10]
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
                bposd_dec_masked.decode(syndrome, mask)
                error ^= bposd_dec_masked.osd0_decoding.astype(np.uint)

            error ^= (np.random.random(Hx.shape[1]) < p).astype(np.uint8)
            syndrome = (Hx @ error) % 2

            # bp_dec.decode(syndrome, np.zeros(Hx.shape[0]))
            # error ^= bp_dec.bp_decoding.astype(np.uint8)
            bposd_dec_unmasked.decode(syndrome, np.zeros(Hx.shape[0]))
            error ^= bposd_dec_unmasked.osdw_decoding.astype(np.uint8)

            res = Result(T, m, ell,
                            int(code["a1"]), int(code["a2"]), int(code["a3"]), int(code["b1"]), int(code["b2"]), int(code["b3"]),
                            int(code["emb_m"]),int(code["emb_ell"]),int(code["aa"]),int(code["bb"]),
                            p, p_mask, 1, int(not np.any(error)))
            rs.append(res)

        if (i % 500 == 0):
            save_new_res(res_file_name, rs)
            rs = []


col_names = ["n","k","ell","m","a1","a2","a3","b1","b2","b3","emb_m","emb_ell","aa","bb","skew","range"]
df = pd.read_csv("codes/40_9.code", names=col_names, header=None)
sorted_df = df.sort_values(by=['range', 'skew'], ascending=False)

for index, row in sorted_df[44:].iterrows():
    print(row)
    test_code(row)