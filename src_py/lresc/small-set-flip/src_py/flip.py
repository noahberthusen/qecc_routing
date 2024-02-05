# from configuration_model import*
# from magic_graph import*
# import matplotlib
# import matplotlib.pyplot as plt
import numpy as np
import random
#Doesn't maximize for best possible flip which is something that the quantum code does

class FLIP():
        def __init__(self,n,m,bit_nbhd,check_nbhd):
                self.n = n
                self.m = m
                self.bit_nbhd = bit_nbhd
                self.check_nbhd = check_nbhd
                self.error = np.zeros(n)
                self.no_unsat = np.zeros(n)
                self.valid = np.zeros(n)
                self.synd = np.zeros(m)

        def reset(self):
                for v in np.arange(self.n):
                        self.error[v] = 0
                        self.no_unsat[v] = 0
                        self.valid[v] = 0

                for c in np.arange(self.m):
                        self.synd[c] = 0

        def generate_error(self,p):
                for v in np.arange(self.n):
                        r = np.random.rand()
                        if (r < p):
                                self.error[v] = 1

        def generate_syndrome(self):
                for c in np.arange(self.m):
                        for v in self.check_nbhd[c]:
                                self.synd[c] = (self.synd[c] + self.error[v])%2

        def generate_no_unsat(self):
                for v in np.arange(self.n):
                        self.no_unsat[v] = 0
                        for c in self.bit_nbhd[v]:
                                self.no_unsat[v] += self.synd[c]

        def generate_valid(self):
                for v in np.arange(self.n):
                        self.valid[v] = len(self.bit_nbhd[v]) - 2*self.no_unsat[v]
                        


        def flip(self,v):
                self.error[v] = (self.error[v] + 1)%2
                for c in self.bit_nbhd[v]:
                        self.synd[c] = (self.synd[c] + 1)%2
                        
        # Returns 1 if corrected, 2 if non zero syndrome and 0 if logical error
        def decode(self):
                self.generate_syndrome()
                self.generate_no_unsat()
                self.generate_valid()
                while (np.min(self.valid) < 0):
                        v = self.find_best_var_node()
                        self.flip(v)
                        self.generate_no_unsat()
                        self.generate_valid()

                if sum(self.synd) != 0:
                        return 2
                elif sum(self.error) == 0:
                        return 1
                else:
                        return 0

class FLIP_MAX(FLIP):
        def find_best_var_node(self):
                return(np.argmin(self.valid))

class FLIP_FAST(FLIP):
        def find_best_var_node(self):
                L = [v for v in range(len(self.valid)) if self.valid[v] < 0]
                return L[random.randint(0, len(L) - 1)]
                
# Run the bit-flip algorithm on a classical code
# Output: returns 1 if corrected, 2 if non zero syndrome and 0 if logical error
def run_flip_ccode(ccode,p, algo):
        if algo == 0:
                graph = FLIP_MAX(ccode.n,ccode.m,ccode.bit_nbhd,ccode.check_nbhd)
        elif algo == -1 :
                graph = FLIP_FAST(ccode.n,ccode.m,ccode.bit_nbhd,ccode.check_nbhd)
        graph.reset()
        graph.generate_error(p)
        return graph.decode()

# ###### Tests #######
# n = 240
# m = 200
# dv = 5
# dc = 6

# name = "_".join((str(n),str(m),str(dv),str(dc)))

# bit_nbhd = []
# check_nbhd  = []
# f = open('ccodes/' + name + '.code','r')
# id = f.readline()
# tag = f.readline()
# for i in range(n):
# 	bit_nbhd.append(list(map(int,f.readline().strip(',\n').split(','))))
# tag = f.readline()
# for i in range(m):
# 	check_nbhd.append(list(map(int,f.readline().strip(',\n').split(','))))

# patience = 100000000

# no_runs = 400
# P = [0.001 + 0.005*j for j in np.arange(10)]

# fail_array = []
# #bit_nbhd, check_nbhd = configuration_model(n,m,dv,dc,patience)
# #bit_nbhd, check_nbhd = magic_graph(n,m,dv)

# graph = FLIP(n, m, bit_nbhd, check_nbhd)
# for p in P:
# 	f = 0
# 	graph.reset()
# 	for j in np.arange(no_runs):
# 		graph.generate_error(p)
# 		f += graph.decode()
# 	f = f/no_runs
# 	fail_array.append(f)
# plt.scatter(P,fail_array)

# plt.show()
