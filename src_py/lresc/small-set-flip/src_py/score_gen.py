# TODO:
# - Change the rule to decide when to flip qbits (weighted using degrees dv and dc; heuristics;...)
# - Improve the efficiency of the function 'score_gen' (NP-complete? Very close to Max-Cut)
# - More tests for the function 'score_gen'
# - More tests for Logical2
# - Improve the function 'find_best_gen'
# - Rewrite the function 'compute_synd' more efficiently

# - Parallelize
# - Profiling

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

def flip_array_to_flip_list(flip_array):
    return [j for j in range(len(flip_array)) if flip_array[j]]

############ Tests ############
# print(compute_gray_code(4))
###############################

###############################
# Let \Delta := |synd| - |synd xor sigma(F)| the decreasing in syndrome weight when we flip the qubits F
# Let |F| be the usual cardinality of F
# Let ||F|| := dv * |F \cap V^2| + dc * |F \cap C^2|
# - Algo 1: we optimize \Delta / ||F||
# - Algo 2: we optimize \Delta / |F|
# - Algo 3: we optimize \Delta
###############################

def flip_condition(algo, synd_diff, s, denom, dv):
    """
    Input:
      'algo' is the algorithm we are using
      'synd_diff' is \Delta for some flip F of rows and collumns
      's' is the decreasing in the syndrome size if we flip a new row \notin F
      'denom' is the denominator in the objective function hen we flip F
      'dv' is the left degree
    Output:
      True if flipping the new row increases the score
      False otherwise
    """
    if algo == 1:
        return s * denom >= dv * synd_diff
    elif algo == 2:
        return s * denom >= synd_diff
    elif algo == 3:
        return s >= 0


def update_denom(algo, denom, dv):
    """
    Input:
      'denom' is the denominator in the objective function for a given F
      'dv' is the left degree
    Output:
      The denominator in the objective function if a new row \notin F
    """
    if algo == 1:
        return denom + dv
    elif algo == 2:
        return denom + 1
    elif algo == 3:
        return 1


def initialize_denom(algo, hor_weight, dc):
    """
    Input:
      'hor_weight' is the number of collumns that we flip
      'dc' is the right degree
    Output:
      The denominator in the objective function if we flip only columns
    """
    if algo == 1:
        return dc * hor_weight
    elif algo == 2:
        return hor_weight
    elif algo == 3:
        return 1
    

def hor_subset_score(algo, hor_synd_diff, hor_weight, ver_synd_diff,dv ,dc):
    """
    'hor' means horizontal
    'ver' means vertical
    Input:
      'hor_synd_diff' is |s| - |s xor synd(F)| for the current F which is a flip of the columns
      'hor_weight' is the size of F
      ver_synd_diff[i] is the difference of syndrome size when we flip the line 'i'
    Output:
      'ver_flips' the optimal set of lines to flip for this given flips of columns
      'synd_diff' the syndrome difference for this flips
    When hor_weight = 0, i.e F = 0 then len(ver_flips) > 0 even if the flipping ver_flips increases the syndrome weight
    """
    synd_diff = hor_synd_diff
    # print(synd_diff, hor_synd_diff, ver_synd_diff)
    ver_flips = []
    sorted_ver_synd_diff = [(ver_synd_diff[i],i) for  i in range(len(ver_synd_diff))]
    sorted_ver_synd_diff.sort(reverse=True)
    
    denom = initialize_denom(algo, hor_weight, dc)
    for (s,i) in sorted_ver_synd_diff:
        if flip_condition(algo, synd_diff, s, denom, dv):
            synd_diff = synd_diff + s
            ver_flips.append(i)
            denom = update_denom(algo, denom, dv)

    return (synd_diff, denom, ver_flips)




def score_gen(algo, synd_gen, gray_code):
    """
    Input:
      'synd_gen' is a 0,1 matrix which reprensents the syndrome of the current generator
      'gray_code' is the output of 'compute_gray_code'
    Output:
      'best_flips' = (ver_flips,hor_flips) are two lists of lines and columns which are optimal for the generator
      'best_synd_diff' is the syndrome difference for these flips
      'best_denom' = len(ver_flips) + len(hor_flips)
    We go through all the possible flips of columns and use the function 'hor_subset_score'
    At the end, we flip at least one qubit even if it is better to flip nothing
    """
    dc = len(synd_gen)
    dv = len(synd_gen[0])
    hor_weight = 0
    hor_flips_array = [False for j in range(dv)]
    hor_synd_diff = 0
    ver_synd_diff = [0 for i in range(dc)]
    for i in range(dc):
        for j in range(dv):
            ver_synd_diff[i] = ver_synd_diff[i] + 2*synd_gen[i][j] - 1

    (best_synd_diff, best_denom, ver_flips) = hor_subset_score(algo, hor_synd_diff, hor_weight, ver_synd_diff, dv, dc)
    best_flips = (ver_flips, [])
    for j in gray_code:
        if hor_flips_array[j]:
            hor_weight = hor_weight - 1
            hor_flips_array[j] = False
            for i in range(dc):
                ver_synd_diff[i] = ver_synd_diff[i] + 4*synd_gen[i][j] - 2
                hor_synd_diff = hor_synd_diff - 2*synd_gen[i][j] + 1
        else:
            hor_weight = hor_weight + 1
            hor_flips_array[j] = True
            for i in range(dc):
                ver_synd_diff[i] = ver_synd_diff[i] - 4*synd_gen[i][j] + 2
                hor_synd_diff = hor_synd_diff + 2*synd_gen[i][j] - 1
            
        (synd_diff, denom, ver_flips) = hor_subset_score(algo, hor_synd_diff, hor_weight, ver_synd_diff, dv, dc)
        if synd_diff*best_denom > best_synd_diff*denom:
            best_synd_diff = synd_diff
            best_denom = denom
            best_flips = (ver_flips, flip_array_to_flip_list(hor_flips_array))
            
    return (best_synd_diff,best_denom,best_flips)
            



######################################################################################
# Using real relaxation
######################################################################################
import numpy as np

def bool_to_int_matrix(bool_matrix):
    return [[-2 * int(b) + 1 for b in L] for L in bool_matrix]


def compute_flip_real_approx(synd_gen):
    int_matrix = np.array(bool_to_int_matrix(synd_gen))
    (u,_,vh) = np.linalg.svd(int_matrix, full_matrices = False)
    ver_flip_list = [i for i in np.arange(u.shape[0]) if u[i,0] <= 0]
    hor_flip_list = [j for j in np.arange(vh.shape[1]) if vh[0,j] <= 0]
    return (ver_flip_list,hor_flip_list)

def compute_weight_synd_gen(synd_gen):
    return sum([sum(synd_gen[i]) for i in range(len(synd_gen))])

def compute_synd_diff(synd_gen, flips):
    (ver_flip_list,hor_flip_list) = flips
    new_synd_gen = [[synd_gen[i][j] ^ (i in ver_flip_list) ^ (j in hor_flip_list)  for j in range(len(synd_gen[0]))] for i in range(len(synd_gen))]
    return compute_weight_synd_gen(synd_gen) - compute_weight_synd_gen(new_synd_gen)

def compute_score_real_approx(synd_gen):
    flips = compute_flip_real_approx(synd_gen)
    return (compute_synd_diff(synd_gen, flips),1,flips)

######################################################################################
# To compare the results
######################################################################################

# Do the complement if 0 notin u
def complement(s,dv,dc):
    (a,b,(u,v)) = s
    if not 0 in u:
        return (a,b,([i for i in range(dc) if not i in u],[j for j in range(dv) if not j in v]))
    else:
        return s

# sort the flips
def sort_scores(s):
    (a,b,(u,v)) = s
    u.sort()
    v.sort()
    return (a,b,(u,v))

def print_all_scores(synd_gen):
    dv = len(synd_gen[0])
    dc = len(synd_gen)
    gray_code = compute_gray_code(dv)
    
    s1 = score_gen(1,synd_gen, gray_code)
    s1 = sort_scores(s1)
    print(s1)

    s2 = score_gen(2,synd_gen, gray_code)
    s2 = sort_scores(s2)
    print(s2)

    s3 = score_gen(3,synd_gen, gray_code)
    s3 = sort_scores(s3)
    s3 = complement(s3,dv,dc)
    print(s3)
    
    s_real = compute_score_real_approx(synd_gen)
    s_real = complement(s_real,dv,dc)
    s_real = sort_scores(s_real)
    print(s_real)
    
    print()


######################################################################################
# To do tests
######################################################################################
import random

def random_synd(dv,dc,p):
    return [[random.uniform(0,1) > p for _ in range(dv)] for _ in range(dc)]

def all_tests(dv,dc):
    big_gray_code = compute_gray_code(dv * dc)
    gray_code = compute_gray_code(dv)
    synd_gen = [[False for j in range(dv)] for i in range(dc)]
    no_fails = 0
    for f in big_gray_code:
        synd_gen[f % dc][f // dc] = not synd_gen[f % dc][f // dc]
        (diff3,_,_) = score_gen(3,synd_gen, gray_code)
        (diff_real,_,_) = compute_score_real_approx(synd_gen)
        if diff3 != diff_real:
            print(diff3, diff_real ,end='||', flush = True)
            no_fails = no_fails + 1
        # print_all_scores(synd_gen)
    print("\nno_fails = ", no_fails , "/", len(big_gray_code))


############ Tests ############
all_tests(8,4)
###############################

############ Tests ############
synd_gen = np.array(random_synd(6,10,0.25))
print(synd_gen)
print_all_scores(synd_gen)
###############################
    
abort
        


############ Tests ############
synd_gen = [
    [True,True,True,False],
    [True,True,True,False],
    [True,True,True,False],
    [True,True,True,False],
    [False,False,False,True],
    [False,False,False,True],
    [False,False,False,True],
    [False,False,False,False]
]
print_all_scores(synd_gen)
###############################

############ Tests ############
synd_gen = [
    [True,True,False,True],
    [False,False,True,False],
    [False,False,True,False],
    [False,False,True,False],
    [False,True,False,False],
    [False,False,True,False],
    [False,False,True,False],
    [False,False,True,False]
]
print_all_scores(synd_gen)
###############################

############ Tests ############
synd_gen = [
    [False,False,False,True],
    [False,True,False,True],
    [False,False,False,True],
    [False,False,False,True],
    [True,False,True,True],
    [False,False,False,True],
    [False,False,False,True],
    [False,False,False,True]
]
print_all_scores(synd_gen)
###############################

