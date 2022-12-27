from Grid import Grid
from itertools import chain

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
# print(gens)
# ind_gens = independent_gens(gens)
# print(len(ind_gens))

# for ind_gen in ind_gens:
    # print(ind_gen)
    # flat = sum(list(chain(*ind_gen)), ())
    # print(len(flat)//2 / M**2)
    # print()
# print(grid.greedy_route_set(gens))