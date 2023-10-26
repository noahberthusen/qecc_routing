import numpy as np
from numpy.linalg import matrix_power, matrix_rank
import itertools
import galois
import pandas as pd
from mec import make_circle
from scipy.sparse import csr_matrix
from scipy.stats import skew
import argparse
import os

full_path = os.path.realpath(__file__)
path, filename = os.path.split(full_path)

# import sys
# sys.path.append("ldpc_masked/src/")
# from ldpc.osd import bposd_decoder
# from result import Result, save_new_res, res_to_line


def test_code(code, res_file_name, num_iters):
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
    Hz = np.hstack([B.T, A.T]).astype(int)

    GF = galois.GF(2)
    k = 2 * (Hz.T.shape[1] - matrix_rank(GF(Hz.T)))
    print(f"{k}")


    # if k > 4 or code['p_log'] > 0.98:
    #     # print(code)
    #     print(f"k = {k}, {code['a1']},{code['a2']},{code['a3']},{code['b1']},{code['b2']},{code['b3']}, {code['p_log']}")
    # else:
        # print(f"k = {2 * (Hz.T.shape[1] - matrix_rank(GF(Hz.T)))}, {code['p_log']}")

    # with open(f"./codes/{ell}_{m}/QZ{ell*m*2}_{res_file_name}.mtx", "w", newline='\n') as f:
    #     H = Hz
    #     f.write("%%MatrixMarket matrix coordinate integer general\n")
    #     f.write("\n")
    #     f.write(f"{H.shape[0]} {H.shape[1]} {np.count_nonzero(H)}")
    #     for i in range(H.shape[0]):
    #         for j in range(H.shape[1]):
    #             if H[i][j]: f.write(f"\n{i+1} {j+1} 1")

    # with open(f"./codes/{ell}_{m}/QX{ell*m*2}_{res_file_name}.mtx", "w", newline='\n') as f:
    #     H = Hx
    #     f.write("%%MatrixMarket matrix coordinate integer general\n")
    #     f.write("\n")
    #     f.write(f"{H.shape[0]} {H.shape[1]} {np.count_nonzero(H)}")
    #     for i in range(H.shape[0]):
    #         for j in range(H.shape[1]):
    #             if H[i][j]: f.write(f"\n{i+1} {j+1} 1")

    def embed_code(code, init):
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
        x_checks = np.array([None for i in range(m*ell)])
        z_checks = np.array([None for i in range(m*ell)])

        for i in range(lattice.shape[0]):
            for j in range(lattice.shape[1]):
                if lattice[i][j][0] == "x":
                    x_checks[int(lattice[i][j][1:])] = (i, j)
                elif lattice[i][j][0] == "z":
                    z_checks[int(lattice[i][j][1:])-(m*ell)] = (i, j)

        # x_rs = []
        # for i in range(m*ell):
        #     gen_qbts = qbts[np.where(Hx[i])[0]]
        #     x_rs.append(make_circle(gen_qbts)[2])
        # #     coord = x_checks[i]
        # #     s = 0
        # #     for qbt in gen_qbts:
        # #         s += (abs(coord[0]-qbt[0]) + abs(coord[1]-qbt[1]))
        # #     # x_rs.append(make_circle(gen_qbts)[2])
        # #     x_rs.append(s)
        # # print(min(x_rs))
        # arr = []
        # for i, x in enumerate(x_rs):
        #     if (x <= (min(x_rs))+np.std(x_rs)):
        #         arr.append(x)
        # print(sum(arr)/sum(x_rs))

        # return lattice

    embed_code((int(code["emb_m"]), int(code["emb_ell"]), int(code["aa"]), int(code["bb"])), 0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-l')
    parser.add_argument('-m')
    parser.add_argument('-t', help="Further test the best codes")
    parser.add_argument('-n', help="Number of iterations")
    args = parser.parse_args()


    col_names = ["t","m","ell","a1","a2","a3","b1","b2","b3","emb_m","emb_ell","aa","bb","p_phys","p_mask","no_test","no_success","p_log"]
    df = pd.read_csv(os.path.join(path, f"./results/{args.l}_{args.m}/refined_results.res"), names=col_names, header=0)
    # df = df.sort_values(by=["p_log"], ascending=False)
    # best_df = df[df["p_log"] > 0.95]

    for index, row in df.iterrows():
        test_code(row, index, int(args.n))