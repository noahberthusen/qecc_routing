from Grid import Grid
from itertools import chain
from mec import make_circle

import random
import numpy as np

# class Route:
#     def __init__(self, N, M, k):
#         self.N = N
#         self.M = M
#         self.k = k
#         self.all_points = [(x,y) for x in range(M) for y in range(M)]
#         self.grid = Grid(M, k+1)

    # def randomly_draw_generators(self, beta, gamma):
    #     def in_circle(point, circle):
    #         x, y = point
    #         cx, cy, r = circle
    #         return np.sqrt((x - cx)**2 + (y - cy)**2) <= r

    #     gens = []

    #     N = int(self.M**(2*beta))
    #     L = np.sqrt(2)*(self.M**gamma)
    #     while (len(gens) < N):
    #         cx, cy = random.choice(self.all_points)

    #         in_points = [point for point in self.all_points if in_circle(point, (cx, cy, L))]

    #         if (len(in_points) >= self.k):
    #             points = random.sample(in_points, self.k)
    #             gens.append(points)

    #     return gens

    # def configuration_model():
    #     m = 10
    #     n = m**2

    #     r = np.sqrt(2)*((m/2)**0.4)
    #     deg_v = 4 # w_c. Every bit is in this many checks
    #     deg_c = 5 # w_r. Every check has this many bits in it
    #     num_checks = (n*deg_v)//deg_c
    #     k = n - num_checks

    #     vs = [deg_v for _ in range(n)]
    #     qbts = [(x,y) for x in range(m) for y in range(m)]
    #     pot_qbts = np.ones((num_checks, n))
    #     ops = [[] for i in range(num_checks)]

    #     while (np.count_nonzero(vs)):
    #         if (np.count_nonzero(pot_qbts)):
    #             c_ind = np.random.choice(np.where(pot_qbts.any(axis=1))[0])
    #         else:
    #             # print("Failed")
    #             break

    #         # choose a v that is within the specified radius (from list of potential qbts)
    #         v_ind = np.random.choice(np.nonzero(pot_qbts[c_ind])[0])
    #         ops[c_ind].append(qbts[v_ind])

    #         if (len(ops[c_ind]) == deg_c):
    #             pot_qbts[c_ind, :] = 0
    #         else:
    #             pot_qbts[c_ind][v_ind] = 0

    #         # update potential qbts
    #         for pot_ind, pot in enumerate(pot_qbts[c_ind]):
    #             if (pot and (make_circle(ops[c_ind] + [qbts[pot_ind]])[2] > r)):
    #                 pot_qbts[c_ind][pot_ind] = 0

    #         vs[v_ind] -= 1
    #         if (not vs[v_ind]):
    #             pot_qbts[:, v_ind] = 0


def print_gens(N, M, gens, label):
    lattice = [["." for i in range(N)] for j in range(M)]

    for qbts, gen in gens:
        for qbt in qbts:
            lattice[qbt[0]][qbt[1]] = "q"
        lattice[gen[0]][gen[1]] = label
    for j in range(N):
        for i in range(M):
    # for i in range(M-1, -1, -1):
        # for j in range(N):
            print(lattice[i][j], end="")
        print()


# def independent_gens(gens):
#     ind_sets = []
#     gen_sets = [set(gen) for gen in gens]

#     while gen_sets:
#         ind_set = [gen_sets[0]]
#         gen_sets.pop(0)

#         for i, gen_set in reversed(list(enumerate(gen_sets))):
#             add = True
#             for gen in ind_set:
#                 if (gen & gen_set):
#                     add = False

#             if add:
#                 ind_set.append(gen_set)
#                 gen_sets.pop(i)

#         if (ind_set):
#             ind_sets.append(ind_set)

#     return ind_sets

# gens = randomly_draw_generators(1, 0)

long_gens = [
[[(5, 0), (3, 24), (0, 29), (0, 27)], (0, 0)],[[(5, 24), (2, 29), (2, 27)], (2, 0)],[[(1, 24), (4, 29), (4, 27)], (4, 0)],[[(5, 2), (3, 26), (0, 29)], (0, 2)],[[(5, 26), (2, 29)], (2, 2)],[[(1, 26), (4, 29)], (4, 2)],[[(5, 4), (3, 28), (0, 1)], (0, 4)],[[(5, 28), (2, 1)], (2, 4)],[[(1, 28), (4, 1)], (4, 4)],[[(3, 0), (3, 2), (0, 5)], (3, 29)],[[(5, 0), (5, 2), (0, 29), (2, 5)], (5, 29)],[[(1, 0), (1, 2), (4, 5)], (1, 29)],[[(3, 28), (0, 1)], (3, 25)],[[(5, 28), (2, 1), (0, 25)], (5, 25)],[[(1, 28), (4, 1)], (1, 25)],[[(3, 0), (0, 3)], (3, 27)],[[(5, 0), (2, 3), (0, 27)], (5, 27)],[[(1, 0), (4, 3)], (1, 27)]
]


short_gens = [
[[(3, 0), (5, 6), (0, 3)], (0, 6)],[[(5, 0), (2, 3)], (2, 6)],[[(1, 0), (4, 3)], (4, 6)],[[(3, 2), (5, 8), (0, 5)], (0, 8)],[[(5, 2), (2, 5)], (2, 8)],[[(1, 2), (4, 5)], (4, 8)],[[(3, 4), (5, 10), (0, 7)], (0, 10)],[[(5, 4), (2, 7)], (2, 10)],[[(1, 4), (4, 7)], (4, 10)],[[(3, 6), (5, 12), (0, 9)], (0, 12)],[[(5, 6), (2, 9)], (2, 12)],[[(1, 6), (4, 9)], (4, 12)],[[(3, 8), (5, 14), (0, 11)], (0, 14)],[[(5, 8), (2, 11)], (2, 14)],[[(1, 8), (4, 11)], (4, 14)],[[(3, 10), (5, 16), (0, 13)], (0, 16)],[[(5, 10), (2, 13)], (2, 16)],[[(1, 10), (4, 13)], (4, 16)],[[(3, 12), (5, 18), (0, 15)], (0, 18)],[[(5, 12), (2, 15)], (2, 18)],[[(1, 12), (4, 15)], (4, 18)],[[(3, 14), (5, 20), (0, 17)], (0, 20)],[[(5, 14), (2, 17)], (2, 20)],[[(1, 14), (4, 17)], (4, 20)],[[(3, 16), (5, 22), (0, 19)], (0, 22)],[[(5, 16), (2, 19)], (2, 22)],[[(1, 16), (4, 19)], (4, 22)],[[(3, 18), (5, 24), (0, 21)], (0, 24)],[[(5, 18), (2, 21)], (2, 24)],[[(1, 18), (4, 21)], (4, 24)],[[(3, 20), (5, 26), (0, 23)], (0, 26)],[[(5, 20), (2, 23)], (2, 26)],[[(1, 20), (4, 23)], (4, 26)],[[(3, 22), (5, 28), (0, 25)], (0, 28)],[[(5, 22), (2, 25)], (2, 28)],[[(1, 22), (4, 25)], (4, 28)],[[(3, 4), (0, 7)], (3, 1)],[[(5, 4), (0, 1), (2, 7)], (5, 1)],[[(1, 4), (4, 7)], (1, 1)],[[(3, 6), (0, 9)], (3, 3)],[[(5, 6), (0, 3), (2, 9)], (5, 3)],[[(1, 6), (4, 9)], (1, 3)],[[(3, 8), (0, 11)], (3, 5)],[[(5, 8), (0, 5), (2, 11)], (5, 5)],[[(1, 8), (4, 11)], (1, 5)],[[(3, 10), (0, 13)], (3, 7)],[[(5, 10), (0, 7), (2, 13)], (5, 7)],[[(1, 10), (4, 13)], (1, 7)],[[(3, 12), (0, 15)], (3, 9)],[[(5, 12), (0, 9), (2, 15)], (5, 9)],[[(1, 12), (4, 15)], (1, 9)],[[(3, 14), (0, 17)], (3, 11)],[[(5, 14), (0, 11), (2, 17)], (5, 11)],[[(1, 14), (4, 17)], (1, 11)],[[(3, 16), (0, 19)], (3, 13)],[[(5, 16), (0, 13), (2, 19)], (5, 13)],[[(1, 16), (4, 19)], (1, 13)],[[(3, 18), (0, 21)], (3, 15)],[[(5, 18), (0, 15), (2, 21)], (5, 15)],[[(1, 18), (4, 21)], (1, 15)],[[(3, 20), (0, 23)], (3, 17)],[[(5, 20), (0, 17), (2, 23)], (5, 17)],[[(1, 20), (4, 23)], (1, 17)],[[(3, 22), (0, 25)], (3, 19)],[[(5, 22), (0, 19), (2, 25)], (5, 19)],[[(1, 22), (4, 25)], (1, 19)],[[(3, 24), (0, 27)], (3, 21)],[[(5, 24), (0, 21), (2, 27)], (5, 21)],[[(1, 24), (4, 27)], (1, 21)],[[(3, 26), (0, 29)], (3, 23)],[[(5, 26), (2, 29), (0, 23)], (5, 23)],[[(1, 26), (4, 29)], (1, 23)],
]


short_gens_split = []
long_gens_split = []
for gen in short_gens:
    for qbt in gen[0]:
        short_gens_split.append([[qbt], gen[1]])
for gen in long_gens:
    for qbt in gen[0]:
        long_gens_split.append([[qbt], gen[1]])

r = Grid(30, 6, 2)
# r.set_available_ancillas(needed_ancillas)

# rounds = r.greedy_route_set(short_gens_split)
# rounds = r.greedy_route_set(long_gens_split)


# r.print_grid()
# print(rounds)

# print_gens(30, 6, short_gens[:len(short_gens)//2], "x")

gens = short_gens
for gen in sorted(gens[:len(gens)//2], key=lambda x: x[1]):
    print(gen)
    print_gens(30, 6, [gen], "x")
    print()

# for gen in sorted(gens[len(gens)//2:], key=lambda x: x[1]):
#     print(gen)
#     print_gens(30, 6, [gen], "z")
#     print()
# c = r.route_generator([(3, 6), (2, 11)], prior_dest=(2, 14))
# print(c)
# print(gens)
# ind_gens = independent_gens(gens)
# print(len(ind_gens))

# for ind_gen in ind_gens:
    # print(ind_gen)
    # flat = sum(list(chain(*ind_gen)), ())
    # print(len(flat)//2 / M**2)
    # print()