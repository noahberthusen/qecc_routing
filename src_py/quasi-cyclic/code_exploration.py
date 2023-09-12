import numpy as np
from numpy.linalg import matrix_power, matrix_rank
import itertools
import galois
from mec import make_circle
from scipy.sparse import csr_matrix

GF = galois.GF(2)

def generate_quasi_cyclic_code(ell, m, As, Bs):
    def cyclic_shift_matrix(l):
        arr = np.eye(l, dtype=int)
        return np.roll(arr, axis=1, shift=1)

    x = np.kron(cyclic_shift_matrix(ell), np.eye(m))
    y = np.kron(np.eye(ell), cyclic_shift_matrix(m))

    A1 = matrix_power(x, As[0])
    A2 = matrix_power(y, As[1])
    A3 = matrix_power(y, As[2])
    A = ( A1 + A2 + A3 ) % 2

    B1 = matrix_power(y, Bs[0])
    B2 = matrix_power(x, Bs[1])
    B3 = matrix_power(x, Bs[2])
    B = ( B1 + B2 + B3 ) % 2

    Hx = np.hstack([A, B]).astype(int)
    Hz = np.hstack([B.T, A.T]).astype(int)
    # H = np.vstack([Hx, Hz])


    def has_toric_layout():
        As = [A1 @ A2.T, A2 @ A3.T, A1 @ A3.T]  # A2 @ A3.T cycling up, A3 @ A2.T cycling up, etc.
        Bs = [B1 @ B2.T, B2 @ B3.T, B1 @ B3.T]

        def has_toric_layout1():
            def order(arr):
                for i in range(1, m*ell):
                    if not np.any(np.eye(arr.shape[0]) - np.linalg.matrix_power(arr, i)):
                        return i
                return -1

            Aorders = [order(As[0]), order(As[1]), order(As[2])]
            Borders = [order(Bs[0]), order(Bs[1]), order(Bs[2])]

            pot_orders = []
            for i, Ao in enumerate(Aorders):
                for j, Bo in enumerate(Borders):
                    if (Ao*Bo == m*ell):
                        pot_orders.append((Ao,Bo,i,j))
            return pot_orders

        def has_toric_layout2(pot_order):
            emb_m, emb_ell, A_ind, B_ind = pot_order

            visited_qbts = set()

            ver = csr_matrix(As[A_ind])
            hor = csr_matrix(Bs[B_ind])

            for i in range(emb_m):
                tmp_qbt = (ver**i)[0].indices[0] if i else 0
                for j in range(emb_ell):
                    visited_qbts |= {(hor**j)[tmp_qbt].indices[0] if j else tmp_qbt}

            return len(visited_qbts) == ell*m

        pot_orders = has_toric_layout1()
        for pot_order in pot_orders:
            if has_toric_layout2(pot_order):
                return True
        return False


    k = 2 * (Hz.T.shape[1] - matrix_rank(GF(Hz.T)))
    print(2*ell*m, k, has_toric_layout2())
    if (not k):
        return


    # lattice = np.empty((2*m, 2*ell), dtype=object)
    # lattice[0][0] = "x0"

    # def get_nbr(i, j):
    #     if (i % 2 == 0):
    #         if (j % 2 == 0):
    #             return "x"
    #         else:
    #             return "r"
    #     else:
    #         if (j % 2 == 0):
    #             return "l"
    #         else:
    #             return "z"

    # for i in range(2*m - 1):
    #     for j in range(2*ell):
    #         curr_ind = int(lattice[i][j][1:])

    #         if (i % 2 == 0):
    #             tmp_A = A3.T
    #         else:
    #             tmp_A = A2
    #         if (j % 2 == 0):
    #             tmp_B = B3.T
    #         else:
    #             tmp_B = B2

    #         lattice[(i+1)%(2*m)][j] = f"{get_nbr((i+1)%(2*m), j)}{np.where(tmp_A @ np.eye(m*ell)[int(lattice[i][j][1:])])[0][0]}"
    #         lattice[i][(j+1)%(2*ell)] = f"{get_nbr(i, (j+1)%(2*ell))}{np.where(tmp_B @ np.eye(m*ell)[int(lattice[i][j][1:])])[0][0]}"

    # for i in range(2*m):
    #     for j in range(2*ell):
    #         if (lattice[i][j][0] == "z"):
    #             lattice[i][j] = f"z{int(lattice[i][j][1:]) + m*ell}"
    #         elif (lattice[i][j][0] == "r"):
    #             lattice[i][j] = f"r{int(lattice[i][j][1:]) + m*ell}"

    # qbts = np.array([None for i in range(2*m*ell)])
    # for i in range(2*m):
    #     for j in range(2*ell):
    #         if lattice[i][j][0] == "r" or lattice[i][j][0] == "l":
    #             qbts[int(lattice[i][j][1:])] = (i, j)

    # rs = []
    # for i in range(m*ell):
    #     gen_qbts = qbts[np.where(Hx[i])[0]]
    #     rs.append(make_circle(gen_qbts)[2])
    # for i in range(m*ell):
    #     gen_qbts = qbts[np.where(Hz[i])[0]]
    #     rs.append(make_circle(gen_qbts)[2])

    # scaled_rs = [(r - min(rs)) / (max(rs) - min(rs)) for r in rs]

    # with open('results.txt', 'a') as file:
    #     file.write(f'{2*ell*m},{k},{ell},{m},{A[0]},{A[1]},{A[2]},{B[0]},{B[1]},{B[2]},{sum(scaled_rs)}')


ell = 6
m = 6

for combination in itertools.product(range(m), repeat=3):
    generate_quasi_cyclic_code(ell, m, combination, combination)