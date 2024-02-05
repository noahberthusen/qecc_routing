cimport cython


import random

print("----------------- Using pyx file ----------------")

############ Score generator ############
# The goal is to find for a given generator the subset of this generator with the best score.
# The score of a subset F for a given syndrome s is (|s| - |s xor synd(F)|)/|F|
# We use the fact that we have an hypergraphe product. The syndrome which correspond to a generator is a rectangle of 0,1. We can flip some row and some columns and the goal is the find the flips with the best score.
# This problem seems hard, it is very close to the Max-cut problem for a bipartite graph (each row and each column is interpreted as a vertex of the bipartite graph)
# Given a flip for the columns, it is possible to find the best flip for the lines. So we do an exhaustive search on the collumn.
# The 'gray code' is a way to go through all the elements of {0,1}^k by flipping one bit at each step.

def compute_gray_code(dv):
    """
    returns the gray code
    If you begin with [0 for i in range(dv)] and flip the bits res[0], res[1], ..., res[dv-2] then you go through {0,1}^dv
    """
    res = []
    for i in range(dv):
        res = res + [i] + res
    return res

############ Tests ############
# print(compute_gray_code(4))
###############################

@cython.binding(True)
cdef hor_subset_score(hor_synd_diff, hor_wweight, ver_synd_diff, dv, dc):
    """
    'hor' means horizontal
    'ver' means vertical
    Input:
      'hor_synd_diff' is |s| - |s xor synd(F)| for the current F which is a flip of the columns
      'hor_wweight' is the weighted weight of the horizontal flips = dc * hor_flips
      ver_synd_diff[i] is the difference of syndrome size when we flip the line 'i'
    Output:
      'ver_flips' the optimal set of lines to flip for this given flips of columns
      'synd_diff' the syndrome difference for this flips
    When hor_weight = 0, i.e F = 0 then len(ver_flips) > 0 even if the flipping ver_flips increases the syndrome weight
    wweight = weighted weight = dv * len(ver_flips) + dc * len (hor_flips)
    """
    cdef int synd_diff = 0
    cdef int s = 0
    cdef int i = 0
    cdef int weight = 0
    cdef int len_ver_flips = 0

    synd_diff = hor_synd_diff
    # print(synd_diff, hor_synd_diff, ver_synd_diff)
    ver_flips = []
    sorted_ver_synd_diff = [(ver_synd_diff[i],i) for  i in range(len(ver_synd_diff))]
    sorted_ver_synd_diff.sort(reverse=True)

    wweight = hor_wweight
    for (s,i) in sorted_ver_synd_diff:
        if s*wweight >= dv * synd_diff:
            synd_diff = synd_diff + s
            ver_flips.append(i)
            wweight = wweight + dv

    return (synd_diff,ver_flips)

@cython.binding(True)
cdef score_gen(synd_gen, synd_gen_mask):
    """
    Input:
      'synd_gen' is a 0,1 matrix which reprensents the syndrome of the current generator
      'gray_code' is the output of 'compute_gray_code'
    Output:
      'best_flips' = (ver_flips,hor_flips) are two lists of lines and columns which are optimal for the generator
      'best_synd_diff' is the syndrome difference for these flips
      'best_wweight' = dv * len(ver_flips) + dc * len(hor_flips)
    We go through all the possible flips of columns and use the function 'hor_subset_score'
    At the end, best_weight > 0 even it is better to flip nothing
    """
    cdef int i = 0
    cdef int j = 0
    cdef int dc = 0
    cdef int dv = 0

    cdef int hor_wweight = 0
    cdef int hor_synd_diff = 0

    cdef int best_wweight = 0
    cdef int best_synd_diff = 0

    dc = len(synd_gen)
    dv = len(synd_gen[0])
    hor_wweight = 0
    hor_flips_array = [False for j in range(dv)]
    hor_synd_diff = 0
    ver_synd_diff = [0 for i in range(dc)]
    for i in range(dc):
        for j in range(dv):
            if (synd_gen_mask[i][j]):
                ver_synd_diff[i] = ver_synd_diff[i] + 2*synd_gen[i][j] - 1

    gray_code = compute_gray_code(dv)
    (best_synd_diff,ver_flips) = hor_subset_score(hor_synd_diff, hor_wweight, ver_synd_diff, dv, dc)
    best_wweight = len(ver_flips) * dv
    best_flips = (ver_flips, [])
    for j in gray_code:
        if hor_flips_array[j]:
            hor_wweight = hor_wweight - dc
            hor_flips_array[j] = False
            for i in range(dc):
                if (synd_gen_mask[i][j]):
                    ver_synd_diff[i] = ver_synd_diff[i] + 4*synd_gen[i][j] - 2
                    hor_synd_diff = hor_synd_diff - 2*synd_gen[i][j] + 1
        else:
            hor_wweight = hor_wweight + dc
            hor_flips_array[j] = True
            for i in range(dc):
                if (synd_gen_mask[i][j]):
                    ver_synd_diff[i] = ver_synd_diff[i] - 4*synd_gen[i][j] + 2
                    hor_synd_diff = hor_synd_diff + 2*synd_gen[i][j] - 1

        (synd_diff,ver_flips) = hor_subset_score(hor_synd_diff, hor_wweight, ver_synd_diff, dv, dc)
        wweight = hor_wweight + dv * len(ver_flips)
        if synd_diff*best_wweight > best_synd_diff*wweight:
            best_synd_diff = synd_diff
            best_wweight = wweight
            best_flips = (ver_flips, [j for j in range(dv) if hor_flips_array[j]])

    return (best_synd_diff,best_wweight,best_flips)


############ Lookup table ############
class Lookup_table:
    """
    'synd_matrix' is a dv*dc matrix which represents the syndrome of the current error
    self.lookup_table contains the score of the generators
    """
    def __init__(self, ccode, synd_matrix, mask):
        self.ccode = ccode
        self.synd_matrix = synd_matrix
        self.mask = mask

        # The weight should not be necessary but we print it
        self.synd_weight = sum([sum(ll and mm for (ll, mm) in zip(l, m)) for (l, m) in zip(self.synd_matrix, self.mask)])
        self.masked_synd_weight = sum([sum(l) for l in self.synd_matrix])
        # self.gray_code = compute_gray_code(dv)

        cdef int c1 = 0
        cdef int v2 = 0

        self.round = 0
        self.last_update = [[self.round for v2 in range(self.ccode.n)] for c1 in range(self.ccode.m)]

        self.lookup_table = [[(0,0,([],[])) for v2 in range(self.ccode.n)] for c1 in range(self.ccode.m)]
        for c1 in range(self.ccode.m):
            for v2 in range(self.ccode.n):
                self.update_score_generator((c1,v2))

    @cython.binding(True)
    def update_score_generator(self,gen):
        """
        compute the score of the generator 'gen'
        """
        cdef int c1 = 0
        cdef int v2 = 0
        cdef int i = 0
        cdef int j = 0
        (c1,v2) = gen
        ver = self.ccode.check_nbhd[c1]
        hor = self.ccode.bit_nbhd[v2]

        synd_gen = [[self.synd_matrix[ver[i]][hor[j]] for j in range(len(hor))] for i in range(len(ver))]
        synd_gen_mask = [[self.mask[ver[i]][hor[j]] for j in range(len(hor))] for i in range(len(ver))]
        self.lookup_table[c1][v2] = score_gen(synd_gen, synd_gen_mask)

    def find_best_gen(self):
        """
        returns the best generator to flip for the current syndrome
        """
        best_gen = None
        cdef int synd_diff = 0
        cdef int best_synd_diff = 0
        cdef int wweight = 0
        cdef int best_wweight = 1
        cdef int c1 = 0
        cdef int v2 = 0
        for c1 in range(self.ccode.m):
            for v2 in range(self.ccode.n):
                (synd_diff,wweight,_) = self.lookup_table[c1][v2]
                if (best_synd_diff*wweight < synd_diff*best_wweight):
                    best_gen = (c1,v2)
                    best_synd_diff = synd_diff
                    best_wweight = wweight

        return best_gen


    # Not efficient but used few times
    def compute_synd(self, gen, flips):
        """
        Compute the syndrome when we flip the row and columns of 'flips' of the syndrome corresponding to the generator 'gen'
        """
        cdef int c1 = 0
        cdef int v2 = 0
        cdef int i = 0
        cdef int j = 0
        (c1,v2) = gen
        ver = self.ccode.check_nbhd[c1]
        hor = self.ccode.bit_nbhd[v2]
        (ver_flips, hor_flips) = flips

        synd_matrix = [[False for j in range(len(hor))] for i in range(len(ver))]
        for i in ver_flips:
            for j in range(len(hor)):
                synd_matrix[i][j] = not synd_matrix[i][j]
        for j in hor_flips:
            for i in range(len(ver)):
                synd_matrix[i][j] = not synd_matrix[i][j]

        synd = []
        for i in range(len(ver)):
            for j in range(len(hor)):
                if synd_matrix[i][j]:
                    v1 = self.ccode.check_nbhd[c1][i]
                    c2 = self.ccode.bit_nbhd[v2][j]
                    synd.append((v1,c2))
        return synd

    def compute_qbits(self, gen, flips):
        """
        compute the qbits which correspond to the rows and columns of 'flips' for the generator 'gen'
        """
        cdef int c1 = 0
        cdef int v2 = 0
        cdef int i = 0
        cdef int j = 0
        (c1,v2) = gen
        (ver_flips, hor_flips) = flips

        vv_qbits = []
        for i in ver_flips:
            vv_qbits.append((self.ccode.check_nbhd[c1][i], v2))
        cc_qbits = []
        for j in hor_flips:
            cc_qbits.append((c1, self.ccode.bit_nbhd[v2][j]))

        return (vv_qbits,cc_qbits)


    def update(self, gen):
        """
        update the lookup table under the assumption that we flip the best subset of the generator 'gen'
        """
        cdef int c1 = 0
        cdef int c2 = 0
        cdef int v1 = 0
        cdef int v2 = 0

        self.round = self.round + 1
        (c1,v2) = gen
        (synd_diff,_,flips) = self.lookup_table[c1][v2]
        synd = self.compute_synd(gen, flips)

        self.synd_weight = self.synd_weight - synd_diff

        for (v1,c2) in synd:
            self.synd_matrix[v1][c2] = not self.synd_matrix[v1][c2]

        for (v1,c2) in synd:
            for v2 in self.ccode.check_nbhd[c2]:
                for c1 in self.ccode.bit_nbhd[v1]:
                    if self.last_update[c1][v2] != self.round:
                        self.update_score_generator((c1,v2))
                        self.last_update[c1][v2] = self.round

        return self.compute_qbits(gen, flips)


############ Decoder ############
def xerror_to_list(ccode, xerror_array):
    """
    Create two lists of qbits given two 0,1 matrices which represent these qbits
    """
    cdef int v1 = 0
    cdef int v2 = 0
    cdef int c1 = 0
    cdef int c2 = 0

    (vv_xerror_array, cc_xerror_array) = xerror_array
    return (
        [(v1,v2) for v1 in range(ccode.n) for v2 in range(ccode.n) if vv_xerror_array[v1][v2]],
        [(c1,c2) for c1 in range(ccode.m) for c2 in range(ccode.m) if cc_xerror_array[c1][c2]])

def xerror_to_array(ccode, xerror):
    """
    Create two 0,1 matrices of qbits given two lists which represent these qbits
    """
    cdef int v1 = 0
    cdef int v2 = 0
    cdef int c1 = 0
    cdef int c2 = 0

    (vv_xerror, cc_xerror) = xerror
    vv_xerror_array = [[False for v2 in range(ccode.n)] for v1 in range(ccode.n)]
    cc_xerror_array = [[False for c2 in range(ccode.m)] for c1 in range(ccode.m)]
    for (v1, v2) in vv_xerror:
        vv_xerror_array[v1][v2] = not vv_xerror_array[v1][v2]
    for (c1, c2) in cc_xerror:
        cc_xerror_array[c1][c2] = not cc_xerror_array[c1][c2]
    return (vv_xerror_array, cc_xerror_array)


@cython.binding(True)
def decoder(ccode, synd_matrix, mask):
    """
    Run the decoder for a given syndrome
    """
    cdef int v1 = 0
    cdef int v2 = 0
    cdef int c1 = 0
    cdef int c2 = 0

    vv_guessed_xerror = [[False for v2 in range(ccode.n)] for v1 in range(ccode.n)]
    cc_guessed_xerror = [[False for c2 in range(ccode.m)] for c1 in range(ccode.m)]
    lookup_table = Lookup_table(ccode,synd_matrix, mask)

    # print("synd weight:", lookup_table.synd_weight, " masked synd weight:", lookup_table.masked_synd_weight)
    gen = lookup_table.find_best_gen()
    while gen != None:
        (vv_qbits,cc_qbits) = lookup_table.update(gen)
        for (v1,v2) in vv_qbits:
            vv_guessed_xerror[v1][v2] = not vv_guessed_xerror[v1][v2]
        for (c1,c2) in cc_qbits:
            cc_guessed_xerror[c1][c2] = not cc_guessed_xerror[c1][c2]
        # print("synd weight:", lookup_table.synd_weight, " masked synd weight:", lookup_table.masked_synd_weight)
        gen = lookup_table.find_best_gen()
    # print("synd weight:", lookup_table.synd_weight, " masked synd weight:", lookup_table.masked_synd_weight)

    return (lookup_table.synd_weight,xerror_to_list(ccode, (vv_guessed_xerror, cc_guessed_xerror)))

############ Measurements ############
def compute_synd_matrix(ccode, xerror):
    """
    returns the syndrome of the list of errors xerror
    """
    cdef int v1 = 0
    cdef int v2 = 0
    cdef int c1 = 0
    cdef int c2 = 0

    (vv_xerror, cc_xerror) = xerror
    synd_matrix = [[False for c2 in range(ccode.m)] for v1 in range(ccode.n)]

    for (v1,v2) in vv_xerror:
        for c2 in ccode.bit_nbhd[v2]:
            synd_matrix[v1][c2] = not synd_matrix[v1][c2]
    for (c1,c2) in cc_xerror:
        for v1 in ccode.check_nbhd[c1]:
            synd_matrix[v1][c2] = not synd_matrix[v1][c2]

    return synd_matrix

def compute_syndrome_weight(ccode, xerror, mask):
    cdef int synd_weight = 0
    synd = compute_synd_matrix(ccode, xerror)

    for c1 in range(ccode.m):
        for v2 in range(ccode.n):
            synd_weight += synd[v2][c1]*mask[v2][c1]
    return synd_weight

############ Gaussian elim ############

def first_line_true(M, rank, j):
    cdef int i = 0
    for i in range(rank, len(M)):
        if M[i][j]:
            return i

# rref = reduced row echelon form
def compute_rref(M0):
    """
    Compute the reduced row echelon form of the matrix M0
    """
    #print_erase = Print_erase(0)

    cdef int n = 0
    cdef int m = 0
    cdef int rank = 0
    cdef int i = 0
    cdef int j = 0

    M = M0.copy()
    n = len(M)
    m = len(M[0])

    for j in range(m):
        #print_erase.print("Gaussian elim: " + str(j) + '/' + str(m))
        i0 = first_line_true(M, rank, j)
        if i0 != None:
            # TODO: do not exchange the rows to optimize
            M[i0], M[rank] = M[rank], M[i0]
            for i in range(rank+1,n):
                if M[i][j]:
                    for j2 in range(j,m):
                        M[i][j2] = M[i][j2] ^ M[rank][j2]
            rank = rank + 1

    return M


def first_true(L, j_min):
    cdef j = 0
    for j in range(j_min, len(L)):
        if L[j]:
            return j

def is_spanned(M, v):
    """
    Return True when the the vector v is spanned by the rows of M
    BE CARREFUL ASSUMPTION: M is rref
    lgu = last generator used

    """
    cdef int m = 0
    cdef int j = 0
    m = len(v)
    if (len(M[0]) != m):
        raise NameError('Bad dimensions')

    cdef int lgu = -1
    cdef int pivot_lgu = -1
    for j in range(m):
        if v[j]:
            while(pivot_lgu < j):
                if lgu == len(M) - 1:
                    return False
                lgu = lgu + 1
                pivot_lgu = first_true(M[lgu], pivot_lgu + 1)
                if pivot_lgu == None:
                    return False
            if (pivot_lgu > j):
                return False
            else:
                for j in range(pivot_lgu, m):
                    v[j] = v[j] ^ M[lgu][j]
    return True


# To do tests
def is_full_rank(gen_matrix_rref):
    """
    Return True when the matrix in rref 'gen_matrix_rref' has rank equal to number of rows
    """
    return sum(gen_matrix_rref[len(gen_matrix_rref) - 1]) != 0




############ Logical error test 1st version ############
# The second version below is more efficient

class Logical1:
    def __init__(self, ccode):
        self.no_qbits = ccode.n * ccode.n + ccode.m * ccode.m
        self.no_checks = ccode.n * ccode.m
        self.no_gen = ccode.n * ccode.m
        self.ccode = ccode
        self.gen_matrix_rref = self.compute_gen_matrix_rref(self.ccode)

    def vc_to_index(self,v1,c2):
        return self.ccode.m * v1 + c2

    def cv_to_index(self,c1,v2):
        return self.ccode.n * c1 + v2

    def vv_to_index(self,v1,v2):
        return v1 * self.ccode.n + v2

    def cc_to_index(self,c1,c2):
        return self.ccode.n * self.ccode.n + c1 * self.ccode.m + c2

    def compute_gen_matrix_rref(self,ccode):
        cdef int c1 = 0
        cdef int c2 = 0
        cdef int v1 = 0
        cdef int v2 = 0

        cdef int qubit = 0
        cdef int gen = 0

        gen_matrix = [[False for qbit in range(self.no_qbits)] for gen in range(self.no_gen)]

        for c1 in range(ccode.m):
            for v2 in range(ccode.n):
                gen = self.cv_to_index(c1,v2)

                for v1 in ccode.check_nbhd[c1]:
                    vv_qbit = self.vv_to_index(v1,v2)
                    gen_matrix[gen][vv_qbit] = not gen_matrix[gen][vv_qbit]
                for c2 in ccode.bit_nbhd[v2]:
                    cc_qbit = self.cc_to_index(c1,c2)
                    gen_matrix[gen][cc_qbit] = not gen_matrix[gen][cc_qbit]

        return compute_rref(gen_matrix)


    def test(self, xerror):
        cdef int c1 = 0
        cdef int c2 = 0
        cdef int v1 = 0
        cdef int v2 = 0

        cdef int qbit = 0

        (vv_xerror, cc_xerror) = xerror
        xerror_array = [False for qbit in range(self.no_qbits)]

        for (v1, v2) in vv_xerror:
            xerror_array[self.vv_to_index(v1,v2)] = not xerror_array[self.vv_to_index(v1,v2)]
        for (c1, c2) in cc_xerror:
            xerror_array[self.cc_to_index(c1,c2)] = not xerror_array[self.cc_to_index(c1,c2)]

        return is_spanned(self.gen_matrix_rref, xerror_array)


############ Logical error test 2nd version ############

# To compute whether an error is trivial or is a logical error we do not compute the rref of the quantum code but only the rref of the classical code
# rref = reduced row echelon form
# cpc = classical parity check
# i = identity
# cpci = classical parity check and identity = (H|I)
# ctpc = classical transpose parity check matrix

def compute_cpci_rref(ccode):
    cpci = [[False for bit in range(ccode.n)] + [i == check for i in range(ccode.m)] for check in range(ccode.m)]
    for bit in range(ccode.n):
        for check in ccode.bit_nbhd[bit]:
            cpci[check][bit] = not cpci[check][bit]
    cpci_rref = compute_rref(cpci)
    cpc_rref = [cpci_rref[check][:ccode.n] for check in range(ccode.m)]
    i_rref = [cpci_rref[check][ccode.n:] for check in range(ccode.m)]

    # We don't need ctpc_rref in the case where the classical code is full rank:
    ctpc_rref = None
    # Case where the parity check matrix of the classical code is not full rank:
    if sum(cpc_rref[ccode.m - 1]) == 0:
        ctpc = [[False for check in range(ccode.m)] for bit in range(ccode.n)]
        for bit in range(ccode.n):
            for check in ccode.bit_nbhd[bit]:
                ctpc[bit][check] = not ctpc[bit][check]
        ctpc_rref = compute_rref(ctpc)

    return (cpc_rref, i_rref, ctpc_rref)


class Logical2:
    def __init__(self, ccode):
        self.no_qbits = ccode.n * ccode.n + ccode.m * ccode.m
        self.no_checks = ccode.n * ccode.m
        self.no_gen = ccode.n * ccode.m
        self.ccode = ccode
        (self.cpc_rref, self.i_rref, self.ctpc_rref) = compute_cpci_rref(self.ccode)

    # lgu = last generator used
    def test(self, xerror):
        (vv_xerror, cc_xerror) = xerror
        (vv_xerror_array, cc_xerror_array) = xerror_to_array(self.ccode, (vv_xerror, cc_xerror))

        for v2 in range(self.ccode.n):
            lgu = -1
            pivot_lgu = -1
            for v1 in range(self.ccode.n):
                if vv_xerror_array[v1][v2]:
                    while pivot_lgu != None and pivot_lgu < v1 and lgu < self.ccode.m - 1:
                        lgu = lgu + 1
                        pivot_lgu = first_true(self.cpc_rref[lgu], pivot_lgu + 1)
                    if pivot_lgu != v1:
                        return False
                    for v10 in range(pivot_lgu, self.ccode.n):
                        vv_xerror_array[v10][v2] = vv_xerror_array[v10][v2] ^ self.cpc_rref[lgu][v10]
                    for c2 in self.ccode.bit_nbhd[v2]:
                        for c1 in range(self.ccode.m):
                            cc_xerror_array[c1][c2] = cc_xerror_array[c1][c2] ^ self.i_rref[lgu][c1]

        while lgu < self.ccode.m - 1 and pivot_lgu != None:
            lgu = lgu + 1
            pivot_lgu = first_true(self.cpc_rref[lgu], pivot_lgu + 1)
        # Case where the parity check matrix of the classical code is full rank
        if pivot_lgu != None:
            return sum([sum(cc_xerror_array[c1]) for c1 in range(self.ccode.m)]) == 0
        lgu_c1_ini = lgu
        pivot_lgu_c1_ini = first_true(self.i_rref[lgu_c1_ini], 0)

        lgu_v2 = -1
        pivot_lgu_v2 = -1
        for c2 in range(self.ccode.m):
            lgu_c1 = lgu_c1_ini
            pivot_lgu_c1 = pivot_lgu_c1_ini
            for c1 in range(self.ccode.m):
                if cc_xerror_array[c1][c2]:
                    while pivot_lgu_c1 != None and pivot_lgu_c1 < c1 and lgu_c1 < self.ccode.m - 1:
                        lgu_c1 = lgu_c1 + 1
                        pivot_lgu_c1 = first_true(self.i_rref[lgu_c1], pivot_lgu_c1 + 1)
                    if pivot_lgu_c1 != c1:
                        return False
                    while pivot_lgu_v2 != None and pivot_lgu_v2 < c2 and lgu_v2 < self.ccode.n - 1:
                        lgu_v2 = lgu_v2 + 1
                        pivot_lgu_v2 = first_true(self.ctpc_rref[lgu_v2], pivot_lgu_v2 + 1)
                    if pivot_lgu_v2 != c2:
                        return False
                    for c10 in range(pivot_lgu_c1, self.ccode.m):
                        if self.i_rref[lgu_c1][c10]:
                            for c20 in range(pivot_lgu_v2, self.ccode.m):
                                cc_xerror_array[c10][c20] = cc_xerror_array[c10][c20] ^ self.ctpc_rref[lgu_v2][c20]

        return True

################# Classical code #####################

class Classical_code:
    """
    - n is the number of bits
    - m is the number of checknodes
    - bit_nbhd is a size n list of lists. bit_nbhd[no_bit] is the list
    of checknodes which involve the bit number no_bit.
    - check_nbhd is a size m list of lists. check_nbhd[no_check] is the list
    of bits which are involved in the check number no_check.
    """
    def __init__(self, n, m, bit_nbhd, check_nbhd, id):
        self.n = n
        self.m = m
        self.bit_nbhd = bit_nbhd
        self.check_nbhd = check_nbhd
        self.id = id

############ Random error ############
def random_error(ccode, p):
    """
    Return a random iid error of proba 'p'
    """
    vv_xerror = [(v1,v2) for v1 in range(ccode.n) for v2 in range(ccode.n) if p > random.uniform(0,1)]
    cc_xerror = [(c1,c2) for c1 in range(ccode.m) for c2 in range(ccode.m) if p > random.uniform(0,1)]

    return (vv_xerror, cc_xerror)

def random_sized_error(ccode, s):
    """
    Return a random iid error of proba 'p'
    """
    arr = [i for i in range(ccode.n**2 + ccode.m**2)]
    errors = random.sample(arr, s)
    vv_xerror = [divmod(e, ccode.n) for e in errors if e < ccode.n**2]
    cc_xerror = [divmod(e - ccode.n**2, ccode.m) for e in errors if e > ccode.n**2]

    return (vv_xerror, cc_xerror)

def random_mask(ccode, maskp):
    """
    Return a random mask for the code
    """
    mask = [[False if maskp > random.uniform(0,1) else True for c2 in range(ccode.m)] for v1 in range(ccode.n)]
    return mask


# Output: 1 if corrected, 2 if non zero syndrome and 0 if logical error
def run_algo_qcode(ccode, xerror, mask, logical2):
    synd_matrix = compute_synd_matrix(ccode, xerror)
    (synd_weight,guessed_xerror) = decoder(ccode, synd_matrix, mask)

    if synd_weight != 0:
        return 2, guessed_xerror
    else:
        return int(logical2.test((xerror[0] + guessed_xerror[0], xerror[1] + guessed_xerror[1]))), guessed_xerror
