from itertools import islice
import random
import numpy as np
import copy

class Generator:
    def __init__(self, qbts, dest, i):
        # have to route all qubits. if dest in generator, then it's essentially a free route
        self.qbts = qbts
        self.key = i
        self.qbts_to_route = qbts.copy()
        self.routed = []
        self.dest = dest # could be in generator or not


class Grid:
    # class that contains a grid of qubit (objects?), site objects perhaps containing data qubit and some number of ancillas
    # will be in charge of actual routing and generating sets to route

    def __init__(self, M, N, k):
        # construct a MxN grid of sites each with k ancilla qubits
        self.reset(M, N, k)

    def reset(self, M, N, k):
        self.N = N
        self.M = M
        self.k = k
        self.generators = []
        self.full_chains = {}
        self.bell_pairs = {}
        self.grid = [[k for _ in range(self.N)] for _ in range(self.M)]
        self.dests = [[False for _ in range(self.N)] for _ in range(self.M)]
        self.in_progress = [[0 for _ in range(self.N)] for _ in range(self.M)]


    def find_chain(self, grid, qSite1, qSite2):
        # this function checks to make sure that the intermediate sites have enough ancilla qubits
        # free to handle the Bell pair chain. Also checks endpoints of chain for one available ancilla
        # returns empty chain if impossible

        visited = [[False for x in range(self.N)] for y in range(self.M)]
        parent = [[None for x in range(self.N)] for y in range(self.M)]
        queue = [qSite1]
        visited[qSite1[0]][qSite1[1]] = True

        if ((grid[qSite1[0]][qSite1[1]] <= 0) or
            (grid[qSite2[0]][qSite2[1]] <= 0)):
            return []

        while queue:
            x, y = queue.pop(0)
            if ((x, y) == qSite2):

                curr_site = qSite2
                chain = [curr_site]

                while (not (curr_site == qSite1)):
                    curr_site = parent[curr_site[0]][curr_site[1]]
                    chain.append(curr_site)
                return chain[::-1]

            if (((x, y) != qSite1) and (grid[x][y] < 1)):
                continue
            pot_nbrs = [(x, y+1), (x, y-1), (x+1, y), (x-1, y)]
            for nbr in pot_nbrs:
                new_x, new_y = nbr
                if ((0 <= new_x < self.M) and (0 <= new_y < self.N)):
                    if (not visited[new_x][new_y]):
                        parent[new_x][new_y] = (x, y)
                        queue.append(nbr)
                        visited[new_x][new_y] = True

        return []


    def add_chain(self, chain, gen):
        # can assume that the path parameter is valid, i.e. there are enough ancillas available
        for i, pair in enumerate(chain):
            x, y = pair
            self.grid[x][y] -= 1
        self.full_chains[gen].append(chain)


    def perform_bell_measurements(self):
        # goes through each of the full chains and creates long range bell pairs
        for _, gen in enumerate(self.generators):
            for _, chain in enumerate(self.full_chains[gen.key]):
                for _, pair in enumerate(chain[1:len(chain)-1]):
                    x, y = pair
                    self.grid[x][y] += 1
                self.bell_pairs[gen.key].append([chain[0], chain[-1]])

            self.full_chains[gen.key] = []


    def perform_syndrome_measurements(self):
        # performs syndrome measurement on any fully routed generators and frees their bell pairs (and extra readout ancilla)
        for _, gen in enumerate(self.generators):
            if (not gen.qbts_to_route):
                # free bell pairs
                for _, pair in enumerate(self.bell_pairs[gen.key]):
                    for _, qbt in enumerate(pair):
                        x, y = qbt
                        self.grid[x][y] += 1

                self.bell_pairs[gen.key] = []


    @staticmethod
    def tmp_add_chain(grid, chain):
        # performs add chain on a provided ancilla grid
        for i, pair in enumerate(chain):
            x, y = pair
            if ((i == 0) or (i == len(chain)-1)):
                grid[x][y] -= 1
            else:
                grid[x][y] -= 1


    def print_grid(self):
        for i in range(self.M):
            for j in range(self.N):
                print(self.grid[i][j], end='')
            print()
        print()

    def set_available_ancillas(self, locs):
        for loc in locs:
            x,y = loc[0]
            if loc[1] > self.k: self.grid[x][y] = loc[1]

    @staticmethod
    def tmp_print_grid(grid):
        for i in range(len(grid)):
            for j in range(len(grid[0])):
                print(grid[i][j], end='')
            print()
        print()


    def route_generator(self, gen, prior_dest=None):
        # attempts to route a generator, if so returns chains and dest
        # if not returns most routings it was able to do and dest
        # able to provide pre-existing destination site, could be in generator or not...

        out = []
        possible_dests = gen
        if (prior_dest):
            possible_dests = [prior_dest]

        for _, dest in enumerate(possible_dests):
            if ((not prior_dest) and self.in_progress[dest[0]][dest[1]]):
                continue

            tmp_grid = copy.deepcopy(self.grid)
            chains = []
            routed_qbts = []
            tot_len = 0

            # dest is the meeting site
            # if (not prior_dest):
            #     tmp_grid[dest[0]][dest[1]] -= 1

            for _, site in enumerate(gen):
                if ((dest != site) and (not self.dests[site[0]][site[1]])):
                    chain = self.find_chain(tmp_grid, dest, site)
                    if chain:
                        self.tmp_add_chain(tmp_grid, chain)
                        chains.append(chain)
                        routed_qbts.append(site)
                        tot_len += len(chain)
                    else:
                        continue
            out.append({"dest":dest, "chains":chains, "routed_qbts":routed_qbts, "num":len(routed_qbts), "len":tot_len})

        if (not prior_dest):
            out.sort(key=lambda x: (x["num"], -x["len"]), reverse=True)

        return out[0] if out else out