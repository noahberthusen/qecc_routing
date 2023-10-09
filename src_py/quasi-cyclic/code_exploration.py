import numpy as np
from numpy.linalg import matrix_power, matrix_rank
import itertools
import galois
from mec import make_circle
from scipy.sparse import csr_matrix
from scipy.stats import skew

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
    if (k < 8): return

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

        rs = []
        for i in range(m*ell):
            gen_qbts = qbts[np.where(Hx[i])[0]]
            rs.append(make_circle(gen_qbts)[2])
        for i in range(m*ell):
            gen_qbts = qbts[np.where(Hz[i])[0]]
            rs.append(make_circle(gen_qbts)[2])

        with open(f"codes/{code_params[0]}_{code_params[1]}.code", "a") as f:
            f.write(f"{2*m*ell},{k},")
            f.write(','.join(map(str, code_params)))
            f.write(',')
            f.write(','.join(map(str, code)))
            f.write(f",{skew(rs)},{max(rs)-min(rs)}\n")

m = 10
ell = 36
variables = range(1, max(m,ell))

for i in range(100000):
    combo = [np.random.randint(ell), np.random.randint(m), np.random.randint(m), np.random.randint(m), np.random.randint(ell), np.random.randint(ell)]
    if ((combo[1] == combo[2]) or (combo[4] == combo[5])): continue
    find_codes([ell,m] + combo)