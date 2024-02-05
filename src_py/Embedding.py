import numpy as np
import random
import itertools
from scipy import sparse
import osqp
from read_code import read_code

class Embedding:

    def __init__(self):
        self.n = 0
        self.grid = None
        self.tanner = None
        self.connect = None

    

    def find_nearest_available(self, x, y):
        visited = np.zeros((self.n,self.n))

        queue = [(x, y)]
        for _ in range(self.n**2):
            x, y = queue.pop(0)
            
            if (not self.grid[x][y]):
                return (x, y)
            
            pot_nbrs = [(x, y+1), (x, y-1), (x+1, y), (x-1, y)]
            for nbr in pot_nbrs:
                new_x, new_y = nbr
                if ((0 <= new_x < self.n) and (0 <= new_y < self.n)):
                    if (not visited[new_x][new_y]):
                        queue.append(nbr)
                        visited[new_x][new_y] = True

    def points_with_manhattan_distance(self, p, d):
        result = []
        for i in range(-d, d+1):
            for j in range(-d, d+1):
                new_x, new_y = (p[0]+i, p[1]+j)
                if ((self.manhattan(p, (new_x, new_y)) <= d)
                    and (0 <= new_x < self.n) and (0 <= new_y < self.n) and (not self.grid[new_x][new_y])):
                    result.append((new_x, new_y))
        return result
    
    def manhattan(self, p1, p2):
        return abs(p1[0] - p2[0]) + abs(p1[1] - p2[1])

    def tot_edge_len(self, G, v, p):
        return sum([self.manhattan(p, (G.nodes[nbr]['x'], G.nodes[nbr]['y'])) for nbr in G.neighbors(v)])
    
    def compact_osqp(self, G, dir, gamma):
        def S(dir, gamma):
            edges = []
            for i in range(len(self.grid[0])):
                if (dir): # horizontal
                    nodes = np.where(self.grid[:,i])[0]
                    gamma_p = min(gamma, len(self.grid[0])//len(nodes)) if len(nodes) else 0
                    for j in range(len(nodes)-1):
                        edges.append((self.grid[nodes[j]][i]-1, self.grid[nodes[j+1]][i]-1, gamma_p))
                else: # vertical
                    nodes = np.where(self.grid[i])[0]
                    gamma_p = min(gamma, len(self.grid[0])//len(nodes)) if len(nodes) else 0
                    for j in range(len(nodes)-1):
                        edges.append((self.grid[i][nodes[j]]-1, self.grid[i][nodes[j+1]]-1, gamma_p))
            return edges
        
        vis = S(dir, gamma)

        A1 = np.zeros((len(vis), G.number_of_nodes()))
        l1 = np.zeros(len(vis))

        for i, e in enumerate(vis):
            A1[i][e[1]] = 1
            A1[i][e[0]] = -1
            l1[i] = e[2]

        A2 = np.eye(G.number_of_nodes())
        l2 = np.zeros(G.number_of_nodes())

        A = np.vstack([A1, A2])
        l = np.concatenate([l1, l2])
        u = np.full(len(vis) + G.number_of_nodes(), int(5*np.sqrt(n))-1)
        
        q = np.zeros(G.number_of_nodes())
        A = sparse.csc_matrix(A)
        l = l
        u = u

        prob = osqp.OSQP()
        prob.setup(P, q, A, l, u, alpha=1.0, verbose=False)
        res = prob.solve()
        
        for i, p in enumerate((res.x).astype(int)):
            v = G.nodes[i]
            if (dir):
                v['x'] = p
            else:
                v['y'] = p
        


    def optimize_embedding(self, G, iters, T, compact_iters):
        k = (0.2/T)**(1/iters)

        for i in range(iters//2):
            for j in range(G.number_of_nodes()):
                if (list(G.neighbors(j))):
                    nbrMedX = np.median([G.nodes[v]['x'] for v in G.neighbors(j)])
                    nbrMedY = np.median([G.nodes[v]['y'] for v in G.neighbors(j)])

                    v = G.nodes[j]
                    x = int(min(max(nbrMedX + random.uniform(-T, T), 0), int(5*np.sqrt(self.n))-1))
                    y = int(min(max(nbrMedY + random.uniform(-T, T), 0), int(5*np.sqrt(self.n))-1))

                    new_x, new_y = self.find_nearest_available(x, y)
                    d = self.manhattan((x,y), (new_x, new_y))
                    pot_locs = self.points_with_manhattan_distance((new_x, new_y), d+1)

                    dist = np.inf
                    for pot_loc in pot_locs:
                        tmp_dist = self.tot_edge_len(G, j, (pot_loc[0], pot_loc[1]))
                        if (tmp_dist < dist):
                            dist = tmp_dist
                            new_x, new_y = pot_loc

                    if (self.tot_edge_len(G, j, (new_x, new_y)) < self.tot_edge_len(G, j, (v['x'], v['y']))):
                        grid[v['x']][v['y']] = 0
                        v['x'] = new_x
                        v['y'] = new_y
                        grid[new_x][new_y] = j+1

            if (i % compact_iters == 0):
                self.compact_osqp(G, compactDir, 3)
                compactDir = not compactDir

                grid = np.array([[0 for i in range(int(5*np.sqrt(self.n)))] for j in range(int(5*np.sqrt(self.n)))])

                for j in range(G.number_of_nodes()):
                    grid[G.nodes[j]['x']][G.nodes[j]['y']] = j+1

            T *= k

        self.compact_osqp(G, True, 3)
        grid = np.array([[0 for i in range(int(5*np.sqrt(self.n)))] for j in range(int(5*np.sqrt(self.n)))])

        for j in range(G.number_of_nodes()):
            grid[G.nodes[j]['x']][G.nodes[j]['y']] = j+1
        self.compact_osqp(G, False, 3)
        grid = np.array([[0 for i in range(int(5*np.sqrt(self.n)))] for j in range(int(5*np.sqrt(self.n)))])

        for j in range(G.number_of_nodes()):
            grid[G.nodes[j]['x']][G.nodes[j]['y']] = j+1

        for i in range(iters//2+1):
            print(i, iters//2+1)
            for j in range(G.number_of_nodes()):
                if (list(G.neighbors(j))):
                    nbrMedX = np.median([G.nodes[v]['x'] for v in G.neighbors(j)])
                    nbrMedY = np.median([G.nodes[v]['y'] for v in G.neighbors(j)])

                    v = G.nodes[j]
                    x = int(min(max(nbrMedX + random.uniform(-T, T), 0), int(5*np.sqrt(self.n))-1))
                    y = int(min(max(nbrMedY + random.uniform(-T, T), 0), int(5*np.sqrt(self.n))-1))

                    new_x, new_y = self.find_nearest_available(x, y)
                    d = self.manhattan((x,y), (new_x, new_y))
                    pot_locs = self.points_with_manhattan_distance((new_x, new_y), d+1)

                    dist = np.inf
                    for pot_loc in pot_locs:
                        tmp_dist = self.tot_edge_len(G, j, (pot_loc[0], pot_loc[1]))
                        if (tmp_dist < dist):
                            dist = tmp_dist
                            new_x, new_y = pot_loc

                    if (self.tot_edge_len(G, j, (new_x, new_y)) < self.tot_edge_len(G, j, (v['x'], v['y']))):
                        grid[v['x']][v['y']] = 0
                        v['x'] = new_x
                        v['y'] = new_y
                        grid[new_x][new_y] = j+1

            if (i % compact_iters == 0):
                self.compact_osqp(G, compactDir, max(1, 1+(2*(iters//2-i-30)/(0.5*iters))))
                compactDir = not compactDir

                grid = np.array([[0 for i in range(int(5*np.sqrt(self.n)))] for j in range(int(5*np.sqrt(self.n)))])

                for j in range(G.number_of_nodes()):
                    grid[G.nodes[j]['x']][G.nodes[j]['y']] = j+1

            T *= k