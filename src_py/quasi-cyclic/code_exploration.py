import numpy as np
from numpy.linalg import matrix_power, matrix_rank
import galois
from mec import make_circle
from scipy.sparse import csr_matrix
import os
import pandas as pd
from ldpc import bposd_decoder
import stim
from scipy.sparse import lil_matrix
from result import Result, save_new_res, res_to_line

full_path = os.path.realpath(__file__)
path, filename = os.path.split(full_path)

GF = galois.GF(2)

def find_codes(code_params):
    def cyclic_shift_matrix(l):
        arr = np.eye(l, dtype=int)
        return np.roll(arr, axis=1, shift=1)

    ell = code_params[0]
    m = code_params[1]

    x = np.kron(cyclic_shift_matrix(ell), np.eye(m))
    y = np.kron(np.eye(ell), cyclic_shift_matrix(m))

    A1 = matrix_power(x, code_params[2])
    A2 = matrix_power(y, code_params[3])
    A3 = matrix_power(y, code_params[4])
    A = ( A1 + A2 + A3 ) % 2

    B1 = matrix_power(y, code_params[5])
    B2 = matrix_power(x, code_params[6])
    B3 = matrix_power(x, code_params[7])
    B = ( B1 + B2 + B3 ) % 2

    Hx = np.hstack([A, B]).astype(int)
    Hz = np.hstack([B.T, A.T]).astype(int)
    k = 2 * (Hz.T.shape[1] - matrix_rank(GF(Hz.T)))
    if (k < 12): return
    # if (k != 8): return

    def has_toric_layout():
        # As = [A1 @ A2.T, A2 @ A3.T, A1 @ A3.T]  # A2 @ A3.T cycling up, A3 @ A2.T cycling up, etc.
        # Bs = [B1 @ B2.T, B2 @ B3.T, B1 @ B3.T]
        As = [A1 @ A2.T, A2 @ A1.T, A2 @ A3.T, A3 @ A2.T, A1 @ A3.T, A3 @ A1.T ]
        Bs = [B1 @ B2.T, B2 @ B1.T, B2 @ B3.T, B3 @ B2.T, B1 @ B3.T, B3 @ B1.T]


        def has_toric_layout1():
            def order(arr):
                for i in range(1, m*ell):
                    if not np.any(np.eye(arr.shape[0]) - np.linalg.matrix_power(arr, i)):
                        return i
                return -1

            Aorders = [order(AA) for AA in As]
            Borders = [order(BB) for BB in Bs]

            pot_orders = []
            for i, Ao in enumerate(Aorders):
                for j, Bo in enumerate(Borders):
                    if (Ao*Bo == m*ell):
                        pot_orders.append((Ao,Bo,i,j))
            return pot_orders

        def has_toric_layout2(pot_codes):
            emb_m, emb_ell, A_ind, B_ind = pot_codes

            visited_qbts = set()

            ver = csr_matrix(As[A_ind])
            hor = csr_matrix(Bs[B_ind])

            for i in range(emb_m):
                tmp_qbt = (ver**i)[0].indices[0] if i else 0
                for j in range(emb_ell):
                    visited_qbts |= {(hor**j)[tmp_qbt].indices[0] if j else tmp_qbt}

            return len(visited_qbts) == ell*m

        confirmed_codes = []
        pot_codes = has_toric_layout1()
        for pot_code in pot_codes:
            if has_toric_layout2(pot_code):
                confirmed_codes.append(pot_code)
        return confirmed_codes

    confirmed_codes = has_toric_layout()
    if (len(confirmed_codes) == 0): return

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

        return lattice

    codes_to_test = []
    for c,code in enumerate(confirmed_codes):
        lattice = embed_code(code, 0)

        colors = np.empty(lattice.shape, dtype=object)
        for i in range(lattice.shape[0]):
            for j in range(lattice.shape[1]):
                if lattice[i][j][0] == "x":
                    colors[i][j] = "red"
                elif lattice[i][j][0] == "r":
                    colors[i][j] = "orange"
                elif lattice[i][j][0] == "l":
                    colors[i][j] = "blue"
                else:
                    colors[i][j] = "green"

        qbts = np.array([None for i in range(2*m*ell)])
        for i in range(lattice.shape[0]):
            for j in range(lattice.shape[1]):
                if lattice[i][j][0] == "r" or lattice[i][j][0] == "l":
                    qbts[int(lattice[i][j][1:])] = (i, j)

        x_rs = np.array([])
        z_rs = np.array([])
        for i in range(m*ell):
            gen_qbts = qbts[np.where(Hx[i])[0]]
            x_rs = np.append(x_rs, make_circle(gen_qbts)[2])
        for i in range(m*ell):
            gen_qbts = qbts[np.where(Hz[i])[0]]
            z_rs = np.append(z_rs, make_circle(gen_qbts)[2])

        x_mask = np.zeros(Hx.shape[0])
        for i, x in enumerate(x_rs):
            if (x > (min(x_rs))+np.std(x_rs)):
                x_mask[i] = 1
        p_mask = np.round(np.count_nonzero(x_mask)/(m*ell), 3)
        adv = (1-p_mask) / (sum(x_rs[x_mask==0])/sum(x_rs))
        if (adv < 1.2): break
        if (p_mask > 0.4): break

        codes_to_test.append(f"{2*m*ell},{k}," + ','.join(map(str, code_params)) + ',' + ','.join(map(str, code)) + f",{p_mask},{adv}")

        # with open(os.path.join(path, f"./codes/{code_params[0]}_{code_params[1]}.code"), "a+") as f:
        #     f.write(f"{2*m*ell},{k},")
        #     f.write(','.join(map(str, code_params)))
        #     f.write(',')
        #     f.write(','.join(map(str, code)))
        #     f.write(f",{p_mask},{adv}\n")
    return codes_to_test

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
    Hz = np.hstack([B.T, A.T]).astype(int)

    bposd_dec = bposd_decoder(
        Hx, # the parity check matrix
        error_rate=0.005,
        channel_probs=[None], #assign error_rate to each qubit. This will override "error_rate" input variable
        max_iter=Hx.shape[1], #the maximum number of iterations for BP)
        bp_method="msl",
        ms_scaling_factor=0, #min sum scaling factor. If set to zero the variable scaling factor method is used
        osd_method="osd_cs", #the OSD method. Choose from:  1) "osd_e", "osd_cs", "osd0"
        osd_order=min(Hx.shape[0],60) #the osd search depth
    )

    num_iters = 10001
    p = 0.005
    p_mask = code["p_mask"] #np.count_nonzero(mask)/Hx.shape[0]
    res_file_name = os.path.join(path, f"./results/{int(code['ell'])}_{int(code['m'])}/phenom_results.res")
    Ts = [10]
    rs = []
    fail = 0

    for i in range(num_iters):
        if fail > 100:
            break
        for T in Ts:
            error = np.zeros(Hx.shape[1]).astype(np.uint8)

            for j in range(T):
                error ^= (np.random.random(Hx.shape[1]) < p).astype(np.uint8)
                syndrome = (Hx @ error) % 2

                bposd_dec.decode(syndrome)
                error ^= bposd_dec.osdw_decoding.astype(np.uint)

            error ^= (np.random.random(Hx.shape[1]) < p).astype(np.uint8)
            syndrome = (Hx @ error) % 2

            # bp_dec.decode(syndrome, np.zeros(Hx.shape[0]))
            # error ^= bp_dec.bp_decoding.astype(np.uint8)
            bposd_dec.decode(syndrome)
            error ^= bposd_dec.osdw_decoding.astype(np.uint8)

            if np.any(error): fail += 1
            res = Result(T, int(code["k"]), ell, m,
                            int(code["a1"]), int(code["a2"]), int(code["a3"]), int(code["b1"]), int(code["b2"]), int(code["b3"]),
                            int(code["emb_m"]),int(code["emb_ell"]),int(code["aa"]),int(code["bb"]),code["adv"],
                            p, p_mask, 1, int(not np.any(error)))
            rs.append(res)

    if fail < 100:
        save_new_res(res_file_name, rs)
        return True
    return False

def test_code_circuit(code):
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

    def par2gen(H):
        GF = galois.GF(2)
        gfH = GF(H)
        gfH_rank = np.linalg.matrix_rank(gfH)

        rref_H = gfH.row_reduce()

        swaps = []
        col_H = rref_H.copy()
        for i in range(gfH_rank):
            inds = np.where(col_H[i])[0]
            pivot = inds[0]
            col_H[:,[i,pivot]] = col_H[:,[pivot,i]]
            swaps.append((i,pivot))

        col_H = col_H[:gfH_rank]
        col_G = GF(np.hstack([col_H[:,gfH_rank:].T, np.eye(H.shape[1]-gfH_rank, dtype=int)]))

        G = col_G.copy()
        for swap in swaps[::-1]:
            G[:,[swap[1],swap[0]]] = G[:,[swap[0],swap[1]]]

        if (np.any(G @ rref_H[:gfH_rank].T) or np.any(col_G @ col_H.T)):
            print("FAILED")
            return
        return (np.array(G, dtype=int), np.array(col_G, dtype=int))
    def commute(x, z, n):
        # 0 if commute, 1 if anticommute
        x1 = x[:n]
        x2 = x[n:]
        z1 = z[:n]
        z2 = z[n:]
        return (x1 @ z2 % 2) ^ (x2 @ z1 % 2)
    def SGSOP(Gx, Gz, n):
        # symplectic gram-schmidt orthogonalization procedure
        sym_Gx = np.hstack([Gx, np.zeros(Gx.shape, dtype=int)])
        sym_Gz = np.hstack([np.zeros(Gz.shape, dtype=int), Gz])
        sym_G = np.vstack([sym_Gx, sym_Gz])
        logicals = []
        generators = []

        while(sym_G.shape[0]):
            g1 = sym_G[0]

            commutes = True
            for i in range(1, sym_G.shape[0]-1):
                g2 = sym_G[i]
                if (commute(g1,g2,n)):
                    logicals.append((g1, g2))
                    sym_G = np.delete(sym_G, [0, i], axis=0)

                    for j in range(sym_G.shape[0]):
                        gj = sym_G[j]
                        sym_G[j] = gj ^ (commute(gj,g2,n) * g1) ^ (commute(gj,g1,n) * g2)
                    commutes = False
                    break

            if commutes:
                generators.append(g1)
                sym_G = np.delete(sym_G, 0, axis=0)

        return (logicals, generators)
    def get_logicals(gen_type=False):
        n = Hx.shape[1]
        Gx, col_Gx = par2gen(Hx)
        Gz, col_Gz = par2gen(Hz)
        logicals, generators = SGSOP(Gx, Gz, n)

        logX = np.array([l[1][n:] for l in logicals])
        logZ = np.array([l[0][:n] for l in logicals])

        if gen_type: return logX
        else: return logZ
    def embed_code(code, init):
        emb_m, emb_ell, A_ind, B_ind = code

        lattice = np.empty((2*emb_m, 2*emb_ell), dtype=object)
        lattice[0][0] = f"x{init}"

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
        return lattice

    lattice = embed_code((int(code["emb_m"]),int(code["emb_ell"]),int(code["aa"]),int(code["bb"])), 0)
    all_qbts = {}
    qbts = np.array([None for i in range(2*m*ell)])
    for i in range(lattice.shape[0]):
        for j in range(lattice.shape[1]):
            if lattice[i][j][0] == "r" or lattice[i][j][0] == "l":
                all_qbts[(i,j)] = int(lattice[i][j][1:])
                qbts[int(lattice[i][j][1:])] = (i, j)
    x_checks = np.array([None for i in range(m*ell)])
    z_checks = np.array([None for i in range(m*ell)])

    for i in range(lattice.shape[0]):
        for j in range(lattice.shape[1]):
            if lattice[i][j][0] == "x":
                all_qbts[(i,j)] = int(lattice[i][j][1:]) + 2*m*ell
                x_checks[int(lattice[i][j][1:])] = (i, j)
            elif lattice[i][j][0] == "z":
                all_qbts[(i,j)] = int(lattice[i][j][1:]) + 2*m*ell
                z_checks[int(lattice[i][j][1:])-(m*ell)] = (i, j)

    x_rs = []
    z_rs = []
    for i in range(m*ell):
        gen_qbts = qbts[np.where(Hx[i])[0]]
        x_rs.append(make_circle(gen_qbts)[2])
    for i in range(m*ell):
        gen_qbts = qbts[np.where(Hz[i])[0]]
        z_rs.append(make_circle(gen_qbts)[2])

    lr_x_checks = []
    sr_x_checks = []
    lr_z_checks = []
    sr_z_checks = []

    for i, x_check in enumerate(x_checks):
        gen_qbts = qbts[np.where(Hx[i])[0]]
        if (x_rs[i] > (min(x_rs)+np.std(x_rs))):
            lr_x_checks.append(i)
        else:
            sr_x_checks.append(i)

    for i, z_check in enumerate(z_checks):
        gen_qbts = qbts[np.where(Hz[i])[0]]
        if (z_rs[i] > min(z_rs)+np.std(z_rs)):
            lr_z_checks.append(i)
        else:
            sr_z_checks.append(i)

    def manhattan(qbts):
        p, q = qbts
        return np.abs(p[0]-q[0])+np.abs(p[1]-q[1])
    def measure_x_checks(checks, scale=False):
        c = stim.Circuit()
        c.append("DEPOLARIZE1", [all_qbts[x_checks[x_check]] for x_check in checks], 0.001)
        for x in checks:
            gen_qbts = qbts[np.where(Hx[x])[0]]
            for qbt in gen_qbts:
                path_qbts = [all_qbts[x_checks[x]], all_qbts[qbt]]
                c.append("CNOT", path_qbts)
                if scale:
                    c.append("DEPOLARIZE2", path_qbts, 0.001*manhattan([x_checks[x], qbt]))
                else:
                    c.append("DEPOLARIZE2", path_qbts, 0.001)
        c.append("H", [all_qbts[x_checks[x_check]] for x_check in checks])
        c.append("DEPOLARIZE1", [all_qbts[x_checks[x_check]] for x_check in checks], 0.001)
        return c
    def measure_z_checks(checks, scale=False):
        c = stim.Circuit()
        for z in checks:
            gen_qbts = qbts[np.where(Hz[z])[0]]
            for qbt in gen_qbts:
                path_qbts = [all_qbts[qbt], all_qbts[z_checks[z]]]
                c.append("CNOT", path_qbts)
                if scale:
                    c.append("DEPOLARIZE2", path_qbts, 0.001*manhattan([qbt, z_checks[z]]))
                else:
                    c.append("DEPOLARIZE2", path_qbts, 0.001)
        return c
    def all_checks():
        c = stim.Circuit()
        c += measure_z_checks(sr_z_checks, True)
        c += measure_z_checks(lr_z_checks, False)
        c += measure_x_checks(sr_x_checks, True)
        c += measure_x_checks(lr_x_checks, False)
        return c
    def inter_detectors(type, checks, meas_offset, prev_meas_offset): # false for z
        c = stim.Circuit()
        for i, check in enumerate(checks):
            coord = x_checks[check] if type else z_checks[check]
            c.append("DETECTOR", [stim.target_rec(-meas_offset+i), stim.target_rec(-prev_meas_offset+i)], (coord[0], coord[1], 0))
        return c
    def observables(type):
        c = stim.Circuit()
        for i, logical in enumerate(get_logicals(type)):
            incl_qbts = np.where(logical)[0]
            incl_qbts = [-j-1 for j in incl_qbts]
            c.append("OBSERVABLE_INCLUDE", [stim.target_rec(j) for j in incl_qbts], i)
        return c


    num_rounds = 10
    lr_time = 5
    num_gen_meas = []

    c = stim.Circuit()
    for key, value in all_qbts.items():
        c.append("QUBIT_COORDS", value, (key[0],key[1],0))
    c.append("R", [qbt for qbt in all_qbts.values()])

    c += all_checks().without_noise()
    c.append("MR", [all_qbts[z_check] for z_check in z_checks])
    c.append("MR", [all_qbts[x_check] for x_check in x_checks])
    num_gen_meas.append(2*m*ell)

    def sr_round():
        c = stim.Circuit()

        c += measure_z_checks(sr_z_checks, True)
        c += measure_x_checks(sr_x_checks, True)

        c.append("X_ERROR", [all_qbts[z_checks[z_check]] for z_check in sr_z_checks], 0.001)
        c.append("X_ERROR", [all_qbts[x_checks[x_check]] for x_check in sr_x_checks], 0.001)
        c.append("MR", [all_qbts[z_checks[z_check]] for z_check in sr_z_checks])
        c.append("MR", [all_qbts[x_checks[x_check]] for x_check in sr_x_checks])
        c.append("X_ERROR", [all_qbts[z_checks[z_check]] for z_check in sr_z_checks], 0.001)
        c.append("X_ERROR", [all_qbts[x_checks[x_check]] for x_check in sr_x_checks], 0.001)

        c += inter_detectors(False, sr_z_checks, len(sr_z_checks+sr_x_checks), num_gen_meas[-1]+len(sr_z_checks+sr_x_checks))
        # c += inter_detectors(True, sr_x_checks, l+l2+num_gen_meas[-1][1]+len(sr_z_checks+sr_x_checks))

        num_gen_meas.append(len(sr_z_checks+sr_x_checks))
        return c

    def lr_round():
        c = stim.Circuit()

        c += measure_z_checks(sr_z_checks, True)
        c += measure_z_checks(lr_z_checks, False)
        c += measure_x_checks(sr_x_checks, True)
        c += measure_x_checks(lr_x_checks, False)

        c.append("X_ERROR", [all_qbts[z_checks[z_check]] for z_check in sr_z_checks+lr_z_checks], 0.001)
        c.append("X_ERROR", [all_qbts[x_checks[x_check]] for x_check in sr_x_checks+lr_x_checks], 0.001)
        c.append("MR", [all_qbts[z_checks[z_check]] for z_check in sr_z_checks+lr_z_checks])
        c.append("MR", [all_qbts[x_checks[x_check]] for x_check in sr_x_checks+lr_x_checks])
        c.append("X_ERROR", [all_qbts[z_checks[z_check]] for z_check in sr_z_checks+lr_z_checks], 0.001)
        c.append("X_ERROR", [all_qbts[x_checks[x_check]] for x_check in sr_x_checks+lr_x_checks], 0.001)


        last_lr = len(num_gen_meas) - num_gen_meas[::-1].index(2*m*ell) - 1
        tot_gens = sum(num_gen_meas[last_lr+1:])

        c += inter_detectors(False, sr_z_checks, 2*m*ell, num_gen_meas[-1]+2*m*ell)
        c += inter_detectors(False, lr_z_checks, 2*m*ell-len(sr_z_checks), tot_gens+3*m*ell+len(lr_z_checks))

        num_gen_meas.append(2*m*ell)
        return c


    for i in range(1,num_rounds+1):
        c.append("SHIFT_COORDS", [], (0,0,1))
        c.append("DEPOLARIZE1", [all_qbts[qbt] for qbt in qbts], 0.001)
        if (i%lr_time==0): c += lr_round()
        else: c += sr_round()

    c += lr_round().without_noise()
    c.append("M",[all_qbts[qbt] for qbt in qbts[::-1]])
    c += observables(False)

    dem = c.detector_error_model()
    pcm = lil_matrix((dem.num_detectors, dem.num_errors), dtype=np.uint8)
    lcm = lil_matrix((dem.num_observables, dem.num_errors), dtype=np.uint8)

    channel_probs = [e.args_copy()[0] for e in c.detector_error_model() if e.type=="error"]
    for i, error_event in enumerate(c.explain_detector_error_model_errors()):
        dets = [det.dem_target.val for det in error_event.dem_error_terms if det.dem_target.is_relative_detector_id()]
        obs = [ob.dem_target.val for ob in error_event.dem_error_terms if ob.dem_target.is_logical_observable_id()]
        pcm[[dets],i] = 1
        lcm[[obs],i] = 1


    bposd_dec = bposd_decoder(
        pcm, # the parity check matrix
        channel_probs=channel_probs, #assign error_rate to each qubit. This will override "error_rate" input variable
        max_iter=pcm.shape[1], #the maximum number of iterations for BP)
        bp_method="ms",
        ms_scaling_factor=0, #min sum scaling factor. If set to zero the variable scaling factor method is used
        osd_method="osd_cs", #the OSD method. Choose from:  1) "osd_e", "osd_cs", "osd0"
        osd_order=min(pcm.shape[0],60) #the osd search depth
    )

    num_iters = 5001
    p = 0.001
    p_mask = code["p_mask"] #np.count_nonzero(mask)/Hx.shape[0]
    res_file_name = os.path.join(path, f"./results/{int(code['ell'])}_{int(code['m'])}/circuit_results.res")
    rs = []

    sampler = c.compile_detector_sampler()
    for i in range(num_iters):
        detection_events, observable_flips = sampler.sample(1, separate_observables=True)
        # guessed_errors = bp_dec.decode(detection_events[0])
        guessed_errors = bposd_dec.decode(detection_events[0])
        guessed_obs = (lcm @ guessed_errors) % 2
        success = np.all(observable_flips[0].astype(int) == guessed_obs)

        res = Result(10, int(code["k"]), ell, m,
                            int(code["a1"]), int(code["a2"]), int(code["a3"]), int(code["b1"]), int(code["b2"]), int(code["b3"]),
                            int(code["emb_m"]),int(code["emb_ell"]),int(code["aa"]),int(code["bb"]),code["adv"],
                            p, p_mask, 1, int(success))
        rs.append(res)
    save_new_res(res_file_name, rs)


m = 6
ell = 12
variables = range(1, max(m,ell))

for i in range(100000):
    if (i % 1000 == 0): print(".", end="")
    combo = [np.random.randint(ell+1), np.random.randint(m+1), np.random.randint(m+1), np.random.randint(m+1), np.random.randint(ell+1), np.random.randint(ell+1)]
    if ((combo[1] == combo[2]) or (combo[4] == combo[5])): continue
    res = find_codes([ell,m] + combo)
    if res: print(len(res))
    if res and len(res):
        col_names = ["n","k","ell","m","a1","a2","a3","b1","b2","b3","emb_m","emb_ell","aa","bb","p_mask","adv"]
        data = [[float(i) for i in r.split(',')] for r in res]
        df = pd.DataFrame(data, columns=col_names)
        for index, row in df.iterrows():
            res2 = test_code(row)
            # if res2: test_code_circuit(row)