from Grid import Grid
from itertools import chain
from mec import make_circle

import random
import numpy as np

class Route:
    def __init__(self, M, k):
        self.M = M
        self.k = k
        self.all_points = [(x,y) for x in range(M) for y in range(M)]
        self.grid = Grid(M, k+1)

    def randomly_draw_generators(self, beta, gamma):
        def in_circle(point, circle):
            x, y = point
            cx, cy, r = circle
            return np.sqrt((x - cx)**2 + (y - cy)**2) <= r

        gens = []
        
        N = int(self.M**(2*beta))
        L = np.sqrt(2)*(self.M**gamma)
        while (len(gens) < N):
            cx, cy = random.choice(self.all_points)

            in_points = [point for point in self.all_points if in_circle(point, (cx, cy, L))]

            if (len(in_points) >= self.k):
                points = random.sample(in_points, self.k)
                gens.append(points)

        return gens

    def configuration_model():
        m = 10
        n = m**2

        r = np.sqrt(2)*((m/2)**0.4)
        deg_v = 4 # w_c. Every bit is in this many checks
        deg_c = 5 # w_r. Every check has this many bits in it
        num_checks = (n*deg_v)//deg_c
        k = n - num_checks

        vs = [deg_v for _ in range(n)]
        qbts = [(x,y) for x in range(m) for y in range(m)]
        pot_qbts = np.ones((num_checks, n))
        ops = [[] for i in range(num_checks)]

        while (np.count_nonzero(vs)):
            if (np.count_nonzero(pot_qbts)):
                c_ind = np.random.choice(np.where(pot_qbts.any(axis=1))[0])
            else:
                # print("Failed")
                break

            # choose a v that is within the specified radius (from list of potential qbts)
            v_ind = np.random.choice(np.nonzero(pot_qbts[c_ind])[0])
            ops[c_ind].append(qbts[v_ind])
            
            if (len(ops[c_ind]) == deg_c):
                pot_qbts[c_ind, :] = 0
            else:
                pot_qbts[c_ind][v_ind] = 0

            # update potential qbts
            for pot_ind, pot in enumerate(pot_qbts[c_ind]):
                if (pot and (make_circle(ops[c_ind] + [qbts[pot_ind]])[2] > r)):
                    pot_qbts[c_ind][pot_ind] = 0

            vs[v_ind] -= 1
            if (not vs[v_ind]):
                pot_qbts[:, v_ind] = 0


def independent_gens(gens):
    ind_sets = []
    gen_sets = [set(gen) for gen in gens]

    while gen_sets:
        ind_set = [gen_sets[0]]
        gen_sets.pop(0)

        for i, gen_set in reversed(list(enumerate(gen_sets))):
            add = True
            for gen in ind_set:
                if (gen & gen_set):
                    add = False

            if add:
                ind_set.append(gen_set)
                gen_sets.pop(i)

        if (ind_set):
            ind_sets.append(ind_set)

    return ind_sets

# gens = randomly_draw_generators(1, 0)
r = Route(10, 5)
# gens = [[(8, 1),(1 ,2),(7, 4),(8, 5),(7, 6)]]
r.grid.greedy_route_set(gens)
# print(gens)
# ind_gens = independent_gens(gens)
# print(len(ind_gens))

# for ind_gen in ind_gens:
    # print(ind_gen)
    # flat = sum(list(chain(*ind_gen)), ())
    # print(len(flat)//2 / M**2)
    # print()