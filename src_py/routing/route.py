from Grid import Grid
from itertools import chain
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
import random

# class Route:
#     def __init__(self, N, M, k):
#         self.N = N
#         self.M = M
#         self.k = k
#         self.all_points = [(x,y) for x in range(M) for y in range(M)]
#         self.grid = Grid(M, k+1)

#     def randomly_draw_generators(self, beta, gamma):
#         def in_circle(point, circle):
#             x, y = point
#             cx, cy, r = circle
#             return np.sqrt((x - cx)**2 + (y - cy)**2) <= r

#         gens = []

#         N = int(self.M**(2*beta))
#         L = np.sqrt(2)*(self.M**gamma)
#         while (len(gens) < N):
#             cx, cy = random.choice(self.all_points)

#             in_points = [point for point in self.all_points if in_circle(point, (cx, cy, L))]

#             if (len(in_points) >= self.k):
#                 points = random.sample(in_points, self.k)
#                 gens.append(points)

#         return gens

#     def configuration_model():
#         m = 10
#         n = m**2

#         r = np.sqrt(2)*((m/2)**0.4)
#         deg_v = 4 # w_c. Every bit is in this many checks
#         deg_c = 5 # w_r. Every check has this many bits in it
#         num_checks = (n*deg_v)//deg_c
#         k = n - num_checks

#         vs = [deg_v for _ in range(n)]
#         qbts = [(x,y) for x in range(m) for y in range(m)]
#         pot_qbts = np.ones((num_checks, n))
#         ops = [[] for i in range(num_checks)]

#         while (np.count_nonzero(vs)):
#             if (np.count_nonzero(pot_qbts)):
#                 c_ind = np.random.choice(np.where(pot_qbts.any(axis=1))[0])
#             else:
#                 # print("Failed")
#                 break

#             # choose a v that is within the specified radius (from list of potential qbts)
#             v_ind = np.random.choice(np.nonzero(pot_qbts[c_ind])[0])
#             ops[c_ind].append(qbts[v_ind])

#             if (len(ops[c_ind]) == deg_c):
#                 pot_qbts[c_ind, :] = 0
#             else:
#                 pot_qbts[c_ind][v_ind] = 0

#             # update potential qbts
#             for pot_ind, pot in enumerate(pot_qbts[c_ind]):
#                 if (pot and (make_circle(ops[c_ind] + [qbts[pot_ind]])[2] > r)):
#                     pot_qbts[c_ind][pot_ind] = 0

#             vs[v_ind] -= 1
#             if (not vs[v_ind]):
#                 pot_qbts[:, v_ind] = 0

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
    if (k == 0): return

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

            zero = np.zeros(m*ell)
            zero[0] = 1

            for i in range(emb_m):
                tmp_qbt = (ver**i) @ zero if i else zero
                for j in range(emb_ell):
                    visited_qbts |= {np.where((hor**j) @ tmp_qbt)[0][0] if j else np.where(tmp_qbt)[0][0]}

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

    def manhattan(qbts):
        p, q = qbts
        return np.abs(p[0]-q[0])+np.abs(p[1]-q[1])

    res = []
    for c, code in enumerate(confirmed_codes):
        lattice = embed_code(code, 0)

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

        gens = []
        tot_paths = 0
        for i in range(m*ell):
            nlqbts = []
            coord = x_checks[i]
            gen_qbts = qbts[np.where(Hx[i])[0]]
            for qbt in gen_qbts:
                if (abs(coord[0]-qbt[0])+ abs(coord[1]-qbt[1]) > 1):
                    # nlqbts.append(qbt)
                    gens.append([[qbt], coord])
                    tot_paths += (manhattan([qbt, coord]) + 1)
            # gens.append([nlqbts, coord])

        grid = Grid(lattice.shape[0], lattice.shape[1], 1)
        rounds = grid.greedy_route_set(gens)

        res.append((tot_paths, lattice.shape[0]*lattice.shape[1], rounds))


        # for i in range(m*ell):
        #     nlqbts = []
        #     coord = z_checks[i]
        #     gen_qbts = qbts[np.where(Hz[i])[0]]
        #     for qbt in gen_qbts:
        #         if (abs(coord[0]-qbt[0])+ abs(coord[1]-qbt[1]) > 1):
        #             nlqbts.append(qbt)


        # print(f"{2*m*ell},{k}," + ','.join(map(str, code_params)) + ',' + ','.join(map(str, code)))
    return res


# m = 20
# ell = 20

# overall_res = []

# variables = range(1, max(m,ell))

# resses = []
# for k in range(1000):
#     if (k % 100 == 0): print(".", end="")
#     combo = [np.random.randint(ell+1), np.random.randint(m+1), np.random.randint(m+1), np.random.randint(m+1), np.random.randint(ell+1), np.random.randint(ell+1)]
#     if ((combo[1] == combo[2]) or (combo[4] == combo[5])): continue
#     res = find_codes([ell,m] + combo)
#     if res: resses += res

# if resses:
#     routing_times = [r[2] for r in resses]
#     theory_times = [r[0]/r[1] for r in resses]


#     overall_res.append([m, ell, np.mean(theory_times), np.mean(routing_times)])

# print(overall_res)

